#!/usr/bin/env python3
"""Exercise a built module through the real OpenCRAVAT BaseAnnotator CLI.

Run with Python from an environment containing open-cravat. Two synthetic CRX
records check an exact BRAF MANE pair and rejection of a different transcript
version. Diagnostic files stay under this directory's test-output/ directory.
"""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
MODULE = 'refseq_transcripts'
REFSEQ = 'NM_004333.6'
TRANSCRIPT = 'ENST00000646891.2'
WRONG_VERSION = 'ENST00000646891.99999'
EXPECTED_TYPES = {
    'uid': 'int', 'primary_refseq': 'string', 'primary_status': 'string',
    'primary_enst': 'string', 'all_refseq': 'string', 'matched_count': 'int',
    'unmatched_count': 'int', 'all_transcripts': 'string', 'gencode_release': 'string',
    'mane_release': 'string',
}


def require(condition, message):
    # Explicit checks also operate when Python is launched with -O.
    if not condition:
        raise ValueError(message)


def verify_reference_pair():
    """The fixed fixture was verified against the official MANE v1.5 summary."""
    summary = HERE / '.cache' / 'MANE.GRCh38.v1.5.summary.txt.gz'
    if not summary.is_file():
        return
    with gzip.open(summary, 'rt', encoding='utf-8') as handle:
        pairs = [(row['RefSeq_nuc'], row['Ensembl_nuc'])
                 for row in csv.DictReader(handle, delimiter='\t')
                 if row['symbol'] == 'BRAF' and row['MANE_status'] == 'MANE Select']
    require(pairs == [(REFSEQ, TRANSCRIPT)], 'Cached MANE v1.5 disagrees with the BRAF fixture.')


def write_input(path):
    columns = [('uid', 'int'), ('chrom', 'string'), ('pos', 'int'),
               ('ref_base', 'string'), ('alt_base', 'string'), ('hugo', 'string'),
               ('transcript', 'string'), ('all_mappings', 'string')]
    with path.open('w', encoding='utf-8') as handle:
        for index, (name, kind) in enumerate(columns):
            handle.write('#column=' + json.dumps(
                {'index': index, 'name': name, 'title': name, 'type': kind}) + '\n')
        for uid, transcript in ((1, TRANSCRIPT), (2, WRONG_VERSION)):
            # OpenCRAVAT all_mappings arrays: protein, protein change, SO,
            # transcript, and coding DNA change. Consequences are input only.
            mappings = {'BRAF': [['ENSP00000493543.1', 'p.Val600Glu', 'MIS',
                                  transcript, 'c.1799T>A']]}
            values = [str(uid), 'chr7', '140753336', 'A', 'T', 'BRAF', transcript,
                      json.dumps(mappings, separators=(',', ':'))]
            handle.write('\t'.join(values) + '\n')


def read_output(path):
    require(path.is_file(), f'OpenCRAVAT did not create {path.name}.')
    columns, records = {}, []
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#column='):
                column = json.loads(line.split('=', 1)[1])
                index = column['index']
                require(isinstance(index, int) and index not in columns,
                        'Output has invalid or duplicate column indices.')
                columns[index] = column
            elif line.strip() and not line.startswith('#'):
                records.append(line.rstrip('\r\n').split('\t'))
    require(sorted(columns) == list(range(len(EXPECTED_TYPES))),
            'Output schema has missing or unexpected columns.')
    require({c['name']: c['type'] for c in columns.values()} == EXPECTED_TYPES,
            'Output column names or types differ from the module schema.')
    require(len(records) == 2, f'Expected two annotation rows, found {len(records)}.')
    require(all(len(row) == len(columns) for row in records), 'Output rows do not match their schema.')
    return [{columns[i]['name']: value for i, value in enumerate(row)} for row in records]


def check_annotation(records, manifest):
    rows = {row['uid']: row for row in records}
    require(set(rows) == {'1', '2'}, 'Variant UIDs were not preserved.')
    for row in records:
        require(row['gencode_release'] == str(manifest['gencode_release']) and
                row['mane_release'] == str(manifest['mane_release']),
                'Output source releases differ from the built module manifest.')
    row = rows['1']
    require(row['primary_enst'] == TRANSCRIPT, 'Original primary transcript was not preserved.')
    require(row['primary_status'] == 'mane_select', 'BRAF exact mapping is not marked MANE Select.')
    accessions = row['primary_refseq'].split(';')
    require(row['primary_refseq'] == REFSEQ, f'Expected only the MANE counterpart {REFSEQ}.')
    require(all(re.fullmatch(r'(?:NM|NR|XM|XR)_\d+\.\d+', a) for a in accessions),
            'Primary transcript output contains an unversioned or non-RefSeq identifier.')
    require(REFSEQ in row['all_refseq'].split(';'), 'BRAF RefSeq is missing from all transcript mappings.')
    require(row['matched_count'] == '1' and row['unmatched_count'] == '0',
            'Exact-match transcript counts are wrong.')
    details = json.loads(row['all_transcripts'])
    require(isinstance(details, list) and len(details) == 1, 'Expected one unique transcript in detailed output.')
    selected = details[0]
    require(selected['gene'] == 'BRAF' and selected['gencode_transcript'] == TRANSCRIPT and
            selected['is_primary'] is True and selected['status'] == 'mane_select',
            'Detailed BRAF transcript identity or status is wrong.')
    require({'refseq': REFSEQ, 'mapping_source': 'MANE Select'} in selected['refseq_mappings'],
            'The exact BRAF MANE Select pair is missing from transcript details.')
    wrong = rows['2']
    require(wrong['primary_enst'] == WRONG_VERSION, 'Unmatched primary transcript was not preserved.')
    require(wrong['primary_status'] == 'no_exact_version_match', 'Different transcript version was not rejected.')
    require(not wrong['primary_refseq'] and not wrong['all_refseq'], 'Different version incorrectly received a RefSeq ID.')
    require(wrong['matched_count'] == '0' and wrong['unmatched_count'] == '1',
            'Version-mismatch transcript counts are wrong.')
    require(json.loads(wrong['all_transcripts']) == [{
        'gene': 'BRAF', 'gencode_transcript': WRONG_VERSION, 'is_primary': True,
        'status': 'no_exact_version_match', 'refseq_mappings': []}],
        'Version-mismatch details do not report an unmapped transcript.')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--module', type=Path, default=HERE / 'build' / MODULE,
                        help='Built refseq_transcripts module directory.')
    parser.add_argument('--timeout', type=int, default=120,
                        help='Maximum annotator subprocess duration in seconds (default: 120).')
    args = parser.parse_args(argv)
    if args.timeout <= 0:
        parser.error('--timeout must be positive.')
    run_dir = None
    try:
        module = args.module.resolve()
        for relative in (f'{MODULE}.py', f'{MODULE}.yml', 'refseq_transcripts_engine.py',
                         'data/mappings.sqlite', 'data/manifest.json'):
            require((module / relative).is_file(), f'Missing built module file: {module / relative}')
        manifest = json.loads((module / 'data' / 'manifest.json').read_text())
        require(manifest.get('mane_release'), 'This smoke test requires a module built with MANE data.')
        verify_reference_pair()
        output_root = HERE / 'test-output'
        output_root.mkdir(exist_ok=True)
        run_dir = Path(tempfile.mkdtemp(prefix='oc-smoke-', dir=output_root))
        modules = run_dir / 'modules'
        annotators = modules / 'annotators'
        annotators.mkdir(parents=True)
        (annotators / MODULE).symlink_to(module, target_is_directory=True)
        conf_dir = run_dir / 'conf'
        conf_dir.mkdir()
        (conf_dir / 'cravat.yml').write_text('{}\n')
        system_conf = conf_dir / 'cravat-system.yml'
        configuration = {'modules_dir': str(modules), 'conf_dir': str(conf_dir), 'save_metrics': False,
                         'store_url': 'https://store.opencravat.org', 'base_modules': []}
        for key, directory in (('jobs_dir', 'jobs'), ('log_dir', 'logs'), ('metrics_dir', 'metrics')):
            path = run_dir / directory
            path.mkdir()
            configuration[key] = str(path)
        system_conf.write_text(json.dumps(configuration))
        source = run_dir / 'braf.crx'
        write_input(source)
        env = {**os.environ, 'OPENCRAVAT_MD': str(modules),
               'OPENCRAVAT_SYSTEM_CONF_PATH': str(system_conf), 'PYTHONDONTWRITEBYTECODE': '1'}
        command = [sys.executable, str(module / f'{MODULE}.py'), str(source),
                   '-d', str(run_dir), '-n', 'braf']
        with (run_dir / 'stdout.txt').open('w') as stdout, (run_dir / 'stderr.txt').open('w') as stderr:
            result = subprocess.run(command, env=env, cwd=run_dir, stdout=stdout, stderr=stderr,
                                    timeout=args.timeout, check=False)
        require(result.returncode == 0,
                f'OpenCRAVAT exited with code {result.returncode}; use Python from the OpenCRAVAT environment.')
        for error_file in run_dir.glob('*.err'):
            require(not error_file.read_text().strip(), f'OpenCRAVAT logged errors in {error_file.name}.')
        for diagnostic in [*run_dir.glob('*.log'), run_dir / 'stdout.txt', run_dir / 'stderr.txt']:
            require(not re.search(r'Traceback \(most recent call last\)|\bERROR\b|\b\w*Exception\b',
                                  diagnostic.read_text()),
                    f'OpenCRAVAT reported an error in {diagnostic.name}.')
        output = run_dir / f'braf.{MODULE}.var'
        check_annotation(read_output(output), manifest)
        print(f'PASS: genuine OpenCRAVAT mapped {TRANSCRIPT} to {REFSEQ} (MANE Select) and rejected a different version.')
        print(f'Output and diagnostic logs: {run_dir}')
    except (OSError, ValueError, KeyError, TypeError, subprocess.TimeoutExpired) as exc:
        parser.exit(1, f'Smoke test failed: {exc}\n' + (f'Diagnostic files: {run_dir}\n' if run_dir else ''))


if __name__ == '__main__':
    main()
