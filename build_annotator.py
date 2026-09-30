#!/usr/bin/env python3
"""Build a RefSeq transcript cross-reference annotator for OpenCRAVAT (stdlib only)."""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import re
import shutil
import sqlite3
import sys
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE = 'refseq_transcripts'
BUILDER_VERSION = '1.1.0'
RNA_ID = re.compile(r'^(?:NM|NR|XM|XR)_\d+\.\d+$')
TX_ID = re.compile(r'^(ENSTR?\d+)\.\d+(_PAR_Y)?$')


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8')


def digest(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def download(url, destination):
    """Cache complete downloads; interrupted transfers never replace a cached file."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return destination
    temporary = destination.with_suffix(destination.suffix + '.part')
    print(f'Downloading {url}', file=sys.stderr)
    try:
        with urllib.request.urlopen(url, timeout=120) as src, open(temporary, 'wb') as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def source_files(args):
    if args.refseq or args.mane:
        if not (args.refseq and (args.mane or args.without_mane)):
            raise ValueError('Local input needs --refseq and either --mane or --without-mane.')
        paths = {'refseq': args.refseq}
        if args.mane:
            paths['mane'] = args.mane
        return paths, {k: {'path': str(p.resolve()), 'sha256': digest(p)} for k, p in paths.items()}
    base = f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{args.release}/'
    name = f'gencode.v{args.release}.metadata.RefSeq.gz'
    md5file = download(base + 'MD5SUMS', args.cache / f'gencode-v{args.release}' / 'MD5SUMS')
    expected = {}
    for line in md5file.read_text().splitlines():
        fields = line.split()
        if len(fields) == 2:
            expected[fields[1].lstrip('*').removeprefix('./')] = fields[0]
    path = download(base + name, args.cache / f'gencode-v{args.release}' / name)
    if name not in expected or digest(path, 'md5') != expected[name]:
        raise ValueError(f'GENCODE MD5 verification failed: {path}; remove that cached file and retry.')
    paths = {'refseq': path}
    sources = {'refseq': {'url': base + name, 'sha256': digest(path), 'publisher_md5': expected[name]}}
    if not args.without_mane:
        name = f'MANE.GRCh38.v{args.mane_release}.summary.txt.gz'
        url = f'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_{args.mane_release}/{name}'
        paths['mane'] = download(url, args.cache / name)
        sources['mane'] = {'url': url, 'sha256': digest(paths['mane'])}
    return paths, sources


def build_database(refseq, mane, destination, release, mane_release='1.5'):
    """Keep versioned transcript pairs from both sources, without fuzzy joins."""
    conn = sqlite3.connect(destination)
    try:
        conn.executescript('''
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE mapping (
                transcript TEXT NOT NULL, transcript_stem TEXT NOT NULL,
                refseq TEXT NOT NULL, source TEXT NOT NULL,
                PRIMARY KEY(transcript, refseq, source));
        ''')
        counts = {'refseq_non_rna': 0, 'refseq_rows': 0, 'mane_rows': 0}

        def insert(tx, accession, source):
            match = TX_ID.fullmatch(tx)
            if not match:
                raise ValueError(f'Expected versioned Ensembl transcript, got {tx!r}.')
            stem = match.group(1) + (match.group(2) or '')
            conn.execute('INSERT OR IGNORE INTO mapping VALUES (?,?,?,?)', (tx, stem, accession, source))

        with open_text(refseq) as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\r\n').split('\t')
                if len(fields) not in {2, 3}:
                    raise ValueError(f'RefSeq metadata line {number}: expected 2 or 3 columns.')
                tx, accession = fields[:2]
                if not RNA_ID.fullmatch(accession):
                    counts['refseq_non_rna'] += 1
                    continue
                insert(tx, accession, 'GENCODE_RefSeq_xref')
                counts['refseq_rows'] += 1
        if counts['refseq_rows'] == 0:
            raise ValueError('No versioned RefSeq RNA mappings found in GENCODE metadata.')
        if mane:
            with open_text(mane) as handle:
                rows = csv.DictReader(handle, delimiter='\t')
                if not {'Ensembl_nuc', 'RefSeq_nuc', 'MANE_status'}.issubset(rows.fieldnames or []):
                    raise ValueError('MANE summary is missing required named columns.')
                for row in rows:
                    tx, accession, status = row['Ensembl_nuc'], row['RefSeq_nuc'], row['MANE_status']
                    if not RNA_ID.fullmatch(accession) or status not in {'MANE Select', 'MANE Plus Clinical'}:
                        raise ValueError(f'Invalid MANE pair: {tx}, {accession}, {status}')
                    insert(tx, accession, status)
                    counts['mane_rows'] += 1
            if not counts['mane_rows']:
                raise ValueError('MANE summary contains no transcript pairs.')
        conn.execute('CREATE INDEX mapping_stem ON mapping(transcript_stem)')
        counts['mapping_rows'] = conn.execute('SELECT COUNT(*) FROM mapping').fetchone()[0]
        counts['ensembl_transcripts'] = conn.execute('SELECT COUNT(DISTINCT transcript) FROM mapping').fetchone()[0]
        metadata = {'schema_version': '1', 'gencode_release': str(release),
                    'mane_release': str(mane_release) if mane else None,
                    'builder_version': BUILDER_VERSION, 'counts': counts}
        conn.executemany('INSERT INTO metadata VALUES (?,?)', [(k, json.dumps(v)) for k, v in metadata.items()])
        conn.commit()
        return metadata
    finally:
        conn.close()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--release', type=int, default=47,
                        help='GENCODE xref release; select to match your mapper (default: 47).')
    parser.add_argument('--mane-release', default='1.5', help='Pinned MANE release (default: 1.5).')
    parser.add_argument('--without-mane', action='store_true', help='Use only GENCODE RefSeq cross-references.')
    parser.add_argument('--refseq', type=Path, help='Local metadata.RefSeq[.gz] from the chosen GENCODE release.')
    parser.add_argument('--mane', type=Path, help='Local MANE GRCh38 summary.txt[.gz].')
    parser.add_argument('--cache', type=Path, default=HERE / '.cache')
    parser.add_argument('--output', type=Path, default=HERE / 'build', help='Parent for the generated module.')
    args = parser.parse_args(argv)
    if args.release < 20:
        parser.error('Use a GENCODE release for GRCh38 (20 or later).')
    if not re.fullmatch(r'\d+\.\d+', args.mane_release):
        parser.error('--mane-release must have the form 1.5.')
    if args.mane and args.without_mane:
        parser.error('--mane and --without-mane are mutually exclusive.')
    target = args.output.resolve() / MODULE
    if target.exists():
        parser.error(f'Output already exists: {target}. Use a new --output directory for a rebuild.')
    try:
        paths, sources = source_files(args)
        args.output.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.refseq-build-', dir=args.output) as stage:
            module_dir = Path(stage) / MODULE
            shutil.copytree(HERE / 'templates', module_dir,
                            ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            shutil.copy2(HERE / 'annotation.py', module_dir / 'refseq_transcripts_engine.py')
            data_dir = module_dir / 'data'
            data_dir.mkdir()
            metadata = build_database(paths['refseq'], paths.get('mane'),
                                      data_dir / 'mappings.sqlite', args.release, args.mane_release)
            manifest = {**metadata, 'sources': sources,
                        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        'database_sha256': digest(data_dir / 'mappings.sqlite')}
            (data_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
            config = module_dir / f'{MODULE}.yml'
            config.write_text(config.read_text().replace('@RELEASE@', str(args.release)).replace(
                '@MANE@', args.mane_release if not args.without_mane else 'none'))
            module_dir.rename(target)
        print(json.dumps({'module': str(target), 'counts': metadata['counts']}, indent=2))
        print('Build complete. See README.md for installation and hg19 input commands.')
    except (OSError, ValueError, sqlite3.Error) as exc:
        parser.exit(1, f'Build failed: {exc}\n')


if __name__ == '__main__':
    main()
