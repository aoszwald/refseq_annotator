"""Offline fixtures for exact ENST-to-RefSeq annotation of mapper output."""
import contextlib
import copy
import gzip
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT))
from annotation import AnnotationDatabase  # noqa: E402
from build_annotator import build_database  # noqa: E402


PRIMARY = 'ENST00000000001.1'
SECONDARY = 'ENST00000000002.3'
OLDER_VERSION = 'ENST00000000003.1'
MANE_ONLY = 'ENST00000000004.2'
PLUS_ONLY = 'ENST00000000005.1'
UNKNOWN = 'ENST00000000006.1'
AMBIGUOUS = 'ENST00000000007.1'
DUPLICATE_PAIR = 'ENST00000000008.1'


def write_inputs(directory):
    """Mapping fixtures include version disagreements and source-only pairs."""
    refseq = directory / 'gencode.v47.metadata.RefSeq'
    mane = directory / 'MANE.GRCh38.v1.5.summary.txt'
    refseq.write_text(
        f'{PRIMARY}\tNM_000001.1\n'
        f'{PRIMARY}\tNM_000001.2\n'
        f'{PRIMARY}\tNM_000002.2\n'
        f'{PRIMARY}\tNM_000002.2\n'
        f'{PRIMARY}\tNP_000001.1\n'
        f'{SECONDARY}\tNR_000003.1\tadditional_metadata\n'
        f'{OLDER_VERSION}\tXM_000005.1\n',
        encoding='utf-8',
    )
    mane.write_text(
        '#NCBI_GeneID\tEnsembl_nuc\tRefSeq_nuc\tMANE_status\n'
        f'1\t{PRIMARY}\tNM_000001.1\tMANE Select\n'
        f'7\t{AMBIGUOUS}\tNM_000008.1\tMANE Select\n'
        f'7\t{AMBIGUOUS}\tNM_000009.1\tMANE Plus Clinical\n'
        f'8\t{DUPLICATE_PAIR}\tNM_000010.1\tMANE Select\n'
        f'8\t{DUPLICATE_PAIR}\tNM_000010.1\tMANE Plus Clinical\n'
        f'4\t{MANE_ONLY}\tNM_000006.2\tMANE Select\n'
        f'5\t{PLUS_ONLY}\tNM_000007.1\tMANE Plus Clinical\n',
        encoding='utf-8',
    )
    return refseq, mane


def mapper_row(transcript, protein='ENSP00000000001.1', achange='p.A1V',
               so='MIS', cchange='c.2C>T'):
    return [protein, achange, so, transcript, cchange]


class AnnotationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix='refseq-test-')
        cls.directory = Path(cls.temporary.name)
        cls.inputs = write_inputs(cls.directory)
        cls.database_path = cls.directory / 'annotations.sqlite'
        with contextlib.redirect_stderr(io.StringIO()):
            cls.metadata = build_database(*cls.inputs, cls.database_path, release=47)
        cls.database = AnnotationDatabase(cls.database_path)

    @classmethod
    def tearDownClass(cls):
        cls.database.close()
        cls.temporary.cleanup()

    def annotate(self, primary=PRIMARY, gene='GENE1', all_mappings=None, **extra):
        data = {'transcript': primary, 'hugo': gene, **extra}
        if all_mappings is not None:
            data['all_mappings'] = all_mappings
        result = self.database.annotate(data)
        return result, {row['gencode_transcript']: row for row in json.loads(result['all_transcripts'])}

    def test_primary_selects_only_mane_and_retains_other_accessions_in_details(self):
        result, transcripts = self.annotate()
        self.assertEqual(result['primary_enst'], PRIMARY)
        self.assertEqual(result['primary_status'], 'mane_select')
        self.assertEqual(result['primary_refseq'], 'NM_000001.1')
        self.assertEqual(set(result['all_refseq'].split(';')),
                         {'NM_000001.1', 'NM_000001.2', 'NM_000002.2'})
        self.assertEqual(transcripts[PRIMARY]['refseq_mappings'], [
            {'refseq': 'NM_000001.1', 'mapping_source': 'MANE Select'},
            {'refseq': 'NM_000001.2', 'mapping_source': 'GENCODE_RefSeq_xref'},
            {'refseq': 'NM_000002.2', 'mapping_source': 'GENCODE_RefSeq_xref'},
        ])
        self.assertTrue(transcripts[PRIMARY]['is_primary'])
        self.assertEqual(transcripts[PRIMARY]['gene'], 'GENE1')

    def test_secondary_mappings_do_not_replace_primary_transcript(self):
        result, transcripts = self.annotate(
            primary=SECONDARY, gene='GENE2',
            all_mappings={'GENE1': [mapper_row(PRIMARY)], 'GENE2': [mapper_row(SECONDARY)]},
        )
        self.assertEqual(result['primary_enst'], SECONDARY)
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(result['primary_status'], 'no_mane_match')
        self.assertEqual(transcripts[SECONDARY]['status'], 'xref_only')
        self.assertTrue(transcripts[SECONDARY]['is_primary'])
        self.assertFalse(transcripts[PRIMARY]['is_primary'])
        self.assertIn('NM_000001.1', result['all_refseq'].split(';'))
        self.assertEqual(result['matched_count'], 2)

    def test_same_gene_mane_transcript_cannot_replace_selected_enst(self):
        result, _ = self.annotate(primary=SECONDARY, all_mappings={
            'GENE1': [mapper_row(PRIMARY), mapper_row(SECONDARY)]})
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(result['primary_status'], 'no_mane_match')

    def test_conflicting_mane_counterparts_leave_primary_empty(self):
        result, transcripts = self.annotate(primary=AMBIGUOUS)
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(result['primary_status'], 'ambiguous_mane_match')
        self.assertEqual(len(transcripts[AMBIGUOUS]['refseq_mappings']), 2)
        self.assertEqual(set(result['all_refseq'].split(';')), {'NM_000008.1', 'NM_000009.1'})

    def test_same_refseq_pair_with_two_mane_tags_is_not_ambiguous(self):
        result, _ = self.annotate(primary=DUPLICATE_PAIR)
        self.assertEqual(result['primary_refseq'], 'NM_000010.1')
        self.assertEqual(result['primary_status'], 'mane_select')

    def test_mane_source_can_add_a_transcript_absent_from_gencode_cross_references(self):
        result, transcripts = self.annotate(primary=MANE_ONLY)
        self.assertEqual(result['primary_status'], 'mane_select')
        self.assertEqual(result['primary_refseq'], 'NM_000006.2')
        self.assertEqual(transcripts[MANE_ONLY]['refseq_mappings'], [
            {'refseq': 'NM_000006.2', 'mapping_source': 'MANE Select'},
        ])

    def test_plus_clinical_has_distinct_status(self):
        result, _ = self.annotate(primary=PLUS_ONLY)
        self.assertEqual(result['primary_status'], 'mane_plus_clinical')
        self.assertEqual(result['primary_refseq'], 'NM_000007.1')

    def test_mismatched_version_is_not_silently_mapped_by_stable_identifier(self):
        result, transcripts = self.annotate(primary='ENST00000000003.2')
        self.assertEqual(result['primary_status'], 'no_exact_version_match')
        self.assertIsNone(result['primary_refseq'])
        self.assertIsNone(result['all_refseq'])
        self.assertEqual(transcripts['ENST00000000003.2']['refseq_mappings'], [])
        self.assertEqual(result['matched_count'], 0)
        self.assertEqual(result['unmatched_count'], 1)

    def test_mane_version_mismatch_is_not_silently_mapped(self):
        result, _ = self.annotate(primary='ENST00000000004.1')
        self.assertEqual(result['primary_status'], 'no_exact_version_match')
        self.assertIsNone(result['primary_refseq'])

    def test_unversioned_input_is_not_given_a_current_version(self):
        result, transcripts = self.annotate(primary='ENST00000000001')
        self.assertEqual(result['primary_status'], 'unversioned_input')
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(transcripts['ENST00000000001']['refseq_mappings'], [])

    def test_unknown_transcript_remains_visible(self):
        result, transcripts = self.annotate(primary=UNKNOWN)
        self.assertEqual(result['primary_status'], 'unmapped')
        self.assertEqual(transcripts[UNKNOWN]['status'], 'unmapped')
        self.assertEqual(transcripts[UNKNOWN]['gene'], 'GENE1')
        self.assertEqual(result['unmatched_count'], 1)

    def test_refseq_input_is_not_reinterpreted_as_ensembl(self):
        result, transcripts = self.annotate(primary='NM_000001.1')
        self.assertEqual(result['primary_status'], 'not_ensembl')
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(transcripts['NM_000001.1']['refseq_mappings'], [])

    def test_missing_primary_can_still_annotate_secondary_transcripts(self):
        result, transcripts = self.annotate(primary='', all_mappings={'GENE2': [mapper_row(SECONDARY)]})
        self.assertEqual(result['primary_status'], 'no_transcript')
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(result['all_refseq'], 'NR_000003.1')
        self.assertEqual(set(transcripts), {SECONDARY})
        self.assertFalse(transcripts[SECONDARY]['is_primary'])

    def test_no_transcripts_has_empty_results(self):
        result, transcripts = self.annotate(primary=None, gene=None)
        self.assertEqual(result['primary_status'], 'no_transcript')
        self.assertEqual(result['matched_count'], 0)
        self.assertEqual(result['unmatched_count'], 0)
        self.assertEqual(transcripts, {})
        self.assertIsNone(result['all_refseq'])

    def test_primary_is_included_even_when_absent_from_all_mappings(self):
        result, transcripts = self.annotate(all_mappings={'GENE2': [mapper_row(SECONDARY)]})
        self.assertEqual(set(transcripts), {PRIMARY, SECONDARY})
        self.assertEqual(result['matched_count'], 2)

    def test_json_and_dictionary_all_mappings_have_the_same_result(self):
        mappings = {'GENE1': [mapper_row(PRIMARY)], 'GENE2': [mapper_row(SECONDARY), mapper_row(UNKNOWN)]}
        self.assertEqual(self.annotate(all_mappings=mappings), self.annotate(all_mappings=json.dumps(mappings)))

    def test_repeated_transcript_consequence_rows_do_not_inflate_counts(self):
        result, transcripts = self.annotate(all_mappings={
            'GENE1': [mapper_row(PRIMARY), mapper_row(PRIMARY, so='SPL', achange=None)],
            'GENE2': [mapper_row(SECONDARY), mapper_row(UNKNOWN), mapper_row(UNKNOWN)],
        })
        self.assertEqual(result['matched_count'], 2)
        self.assertEqual(result['unmatched_count'], 1)
        self.assertEqual(set(transcripts), {PRIMARY, SECONDARY, UNKNOWN})

    def test_ensembl_hgvs_and_consequences_are_neither_mutated_nor_relabelled(self):
        data = {
            'transcript': PRIMARY, 'hugo': 'GENE1', 'cchange': 'c.2C>T', 'achange': 'p.A1V',
            'all_mappings': {'GENE1': [mapper_row(PRIMARY)]},
        }
        original = copy.deepcopy(data)
        result = self.database.annotate(data)
        self.assertEqual(data, original)
        self.assertNotIn('cchange', result)
        self.assertNotIn('achange', result)
        self.assertNotIn('so', result)
        for row in json.loads(result['all_transcripts']):
            self.assertEqual(set(row), {'gene', 'gencode_transcript', 'is_primary', 'status', 'refseq_mappings'})

    def test_malformed_mapper_output_is_rejected_instead_of_silently_ignored(self):
        for mappings in ('not_json', ['unexpected_array'], {'GENE1': 'not_a_list'},
                         {'GENE1': [['too', 'few', 'fields']]},
                         {'GENE1': [[None, None, None, 123]]}):
            with self.subTest(mappings=mappings):
                with self.assertRaises(ValueError):
                    self.annotate(all_mappings=mappings)


class BuilderTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='refseq-builder-test-')
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.inputs = write_inputs(self.directory)

    def build(self, inputs=None):
        with contextlib.redirect_stderr(io.StringIO()):
            return build_database(*(inputs or self.inputs), self.directory / 'annotations.sqlite', release=47)

    def cli(self, *extra):
        refseq, mane = self.inputs
        return subprocess.run(
            [sys.executable, str(PROJECT / 'build_annotator.py'), '--refseq', str(refseq),
             '--mane', str(mane), '--release', '47',
             '--output', str(self.directory / 'build'), *extra],
            capture_output=True, text=True, check=False,
        )

    def test_builder_reads_gzip_and_plain_inputs(self):
        compressed = []
        for original in self.inputs:
            path = original.with_suffix(original.suffix + '.gz')
            with gzip.open(path, 'wt', encoding='utf-8') as handle:
                handle.write(original.read_text())
            compressed.append(path)
        metadata = self.build(compressed)
        self.assertEqual(metadata['gencode_release'], '47')
        self.assertEqual(metadata['mane_release'], '1.5')
        database = AnnotationDatabase(self.directory / 'annotations.sqlite')
        self.addCleanup(database.close)
        result = database.annotate({'transcript': PRIMARY, 'hugo': 'GENE1'})
        self.assertEqual(result['primary_status'], 'mane_select')

    def test_mane_is_optional(self):
        metadata = self.build((self.inputs[0], None))
        self.assertIsNone(metadata['mane_release'])
        database = AnnotationDatabase(self.directory / 'annotations.sqlite')
        self.addCleanup(database.close)
        result = database.annotate({'transcript': PRIMARY, 'hugo': 'GENE1'})
        self.assertEqual(result['primary_status'], 'no_mane_match')
        self.assertIsNone(result['primary_refseq'])
        self.assertEqual(result['all_refseq'], 'NM_000001.1;NM_000001.2;NM_000002.2')

    def test_dataset_with_no_valid_rna_mapping_is_rejected(self):
        self.inputs[0].write_text(f'{PRIMARY}\tNP_000001.1\n')
        with self.assertRaises(ValueError):
            self.build((self.inputs[0], None))

    def test_malformed_refseq_table_is_rejected(self):
        self.inputs[0].write_text('not_a_tab_separated_mapping\n')
        with self.assertRaises(ValueError):
            self.build()

    def test_mane_missing_required_columns_is_rejected(self):
        self.inputs[1].write_text('RefSeq_nuc\tMANE_status\nNM_000001.1\tMANE Select\n')
        with self.assertRaises(ValueError):
            self.build()

    def test_failed_cli_build_does_not_publish_partial_module(self):
        self.inputs[0].write_text('invalid_record\n')
        result = self.cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.directory / 'build' / 'refseq_transcripts').exists())
        self.assertEqual(list((self.directory / 'build').iterdir()), [])

    def test_successful_cli_module_is_complete_and_cannot_be_overwritten(self):
        result = self.cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        module = self.directory / 'build' / 'refseq_transcripts'
        for name in ('refseq_transcripts.py', 'refseq_transcripts_engine.py', 'refseq_transcripts.yml',
                     'data/mappings.sqlite', 'data/manifest.json'):
            self.assertTrue((module / name).is_file(), name)
        config = (module / 'refseq_transcripts.yml').read_text()
        self.assertNotIn('@RELEASE@', config)
        self.assertNotIn('@MANE@', config)
        manifest = json.loads((module / 'data' / 'manifest.json').read_text())
        self.assertEqual(set(manifest['sources']), {'refseq', 'mane'})
        self.assertEqual(len(manifest['database_sha256']), 64)
        before = {p.relative_to(module): p.read_bytes() for p in module.rglob('*') if p.is_file()}
        repeated = self.cli()
        self.assertNotEqual(repeated.returncode, 0)
        self.assertIn('already exists', repeated.stderr)
        after = {p.relative_to(module): p.read_bytes() for p in module.rglob('*') if p.is_file()}
        self.assertEqual(before, after)


if __name__ == '__main__':
    unittest.main()
