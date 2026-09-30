"""Resolve OpenCRAVAT's existing ENST transcripts to versioned NCBI RNA IDs."""
import json
import re
import sqlite3
from functools import lru_cache
from pathlib import Path

SOURCE_ORDER = {'MANE Select': 0, 'MANE Plus Clinical': 1, 'GENCODE_RefSeq_xref': 2}
SOURCE_STATUS = {'MANE Select': 'mane_select', 'MANE Plus Clinical': 'mane_plus_clinical',
                 'GENCODE_RefSeq_xref': 'xref_only'}
TX_ID = re.compile(r'^(ENSTR?\d+)\.\d+(_PAR_Y)?$')


class AnnotationDatabase:
    def __init__(self, path):
        self.conn = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)
        self.conn.row_factory = sqlite3.Row
        try:
            self.metadata = {r['key']: json.loads(r['value']) for r in self.conn.execute('SELECT * FROM metadata')}
            if self.metadata.get('schema_version') != '1':
                raise ValueError('Unsupported RefSeq mapping database schema.')
        except Exception:
            self.conn.close()
            raise
        self.lookup = lru_cache(maxsize=8192)(self._lookup)

    def close(self):
        self.lookup.cache_clear()
        self.conn.close()

    def _lookup(self, transcript):
        if not transcript:
            return 'no_transcript', []
        match = TX_ID.fullmatch(transcript)
        if not match:
            if re.fullmatch(r'ENSTR?\d+(?:_PAR_Y)?', transcript):
                return 'unversioned_input', []
            return 'not_ensembl', []
        pairs = {}
        for row in self.conn.execute('SELECT refseq,source FROM mapping WHERE transcript=?', (transcript,)):
            accession, source = row['refseq'], row['source']
            if accession not in pairs or SOURCE_ORDER[source] < SOURCE_ORDER[pairs[accession]]:
                pairs[accession] = source
        if not pairs:
            stem = match.group(1) + (match.group(2) or '')
            exists = self.conn.execute('SELECT 1 FROM mapping WHERE transcript_stem=? LIMIT 1', (stem,)).fetchone()
            return ('no_exact_version_match' if exists else 'unmapped'), []
        mappings = [{'refseq': accession, 'mapping_source': source}
                    for accession, source in sorted(pairs.items(), key=lambda x: (SOURCE_ORDER[x[1]], x[0]))]
        return SOURCE_STATUS[mappings[0]['mapping_source']], mappings

    def annotate(self, input_data):
        primary = input_data.get('transcript') or ''
        if not isinstance(primary, str):
            raise ValueError('Expected a transcript identifier string.')
        raw = input_data.get('all_mappings') or {}
        raw = json.loads(raw) if isinstance(raw, str) else raw
        if not isinstance(raw, dict):
            raise ValueError('Expected all_mappings to be a gene-to-mappings JSON object.')
        transcripts = set()
        for gene, rows in raw.items():
            if not isinstance(rows, list):
                raise ValueError(f'Expected transcript mapping list for {gene}.')
            for row in rows:
                if not isinstance(row, (list, tuple)) or len(row) < 4 or not isinstance(row[3], str):
                    raise ValueError(f'Invalid OpenCRAVAT transcript mapping for {gene}.')
                if row[3]:
                    transcripts.add((gene, row[3]))
        if primary and not any(tx == primary for _, tx in transcripts):
            transcripts.add((input_data.get('hugo') or '', primary))
        rows = []
        for gene, tx in sorted(transcripts, key=lambda x: (x[1] != primary, x[0], x[1])):
            status, mappings = self.lookup(tx)
            rows.append({'gene': gene, 'gencode_transcript': tx, 'is_primary': tx == primary,
                         'status': status, 'refseq_mappings': mappings})
        primary_status, primary_mappings = self.lookup(primary)
        mane_matches = [m for m in primary_mappings if m['mapping_source'] in
                        {'MANE Select', 'MANE Plus Clinical'}]
        primary_refseq = None
        if len(mane_matches) == 1:
            primary_refseq = mane_matches[0]['refseq']
            primary_status = SOURCE_STATUS[mane_matches[0]['mapping_source']]
        elif len(mane_matches) > 1:
            primary_status = 'ambiguous_mane_match'
        elif primary_mappings:
            primary_status = 'no_mane_match'
        all_refseq = sorted({m['refseq'] for row in rows for m in row['refseq_mappings']})
        matched = {row['gencode_transcript'] for row in rows if row['refseq_mappings']}
        unmatched = {row['gencode_transcript'] for row in rows if not row['refseq_mappings']}
        return {'primary_enst': primary or None,
                'primary_refseq': primary_refseq,
                'primary_status': primary_status,
                'all_refseq': ';'.join(all_refseq) or None,
                'matched_count': len(matched), 'unmatched_count': len(unmatched),
                'all_transcripts': json.dumps(rows, separators=(',', ':')),
                'gencode_release': str(self.metadata['gencode_release']),
                'mane_release': self.metadata['mane_release']}
