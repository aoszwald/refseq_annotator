# NCBI RefSeq transcript cross-references

Adds versioned RefSeq RNA accessions for the primary and all other Ensembl
transcripts already reported by the OpenCRAVAT mapper. Joins require the full
Ensembl transcript version. MANE pairs and ordinary GENCODE RefSeq associations
are distinguished; missing and different-version matches remain explicit.

`primary_refseq` contains only the unique MANE counterpart of the exact selected
ENST version. It is blank with `no_mane_match` for association-only matches, or
`ambiguous_mane_match` for multiple distinct MANE counterparts. All associations
remain in `all_refseq` and `all_transcripts`; their counts include ordinary xrefs.

This adds identifier cross-references. It does not replace the mapper, update
its gene models, or calculate RefSeq HGVS/protein consequences. Ordinary xrefs
are not equivalent transcript models. Existing mapper-derived changes must not
simply be relabeled as RefSeq changes.

Sources, versions, counts and SHA256 hashes are in `data/manifest.json`.
