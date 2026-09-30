# RefSeq transcript annotator for OpenCRAVAT

`build_annotator.py` builds an OpenCRAVAT annotator named **refseq_transcripts**.
It reads the Ensembl transcripts already reported by your mapper and adds their
**versioned NCBI RefSeq RNA accessions**. It handles the selected primary
transcript and all other affected transcripts separately. This project is
standalone; building it does not modify OpenCRAVAT configuration.

This implements the transcript cross-reference workaround. It is **not a RefSeq
mapper**: it does not find additional RefSeq-only transcripts, change the primary
transcript, or recompute coding/protein consequences or HGVS on RefSeq models.

## Why these sources

The defaults are GENCODE **47**'s RefSeq metadata and **MANE 1.5**. These
were tested with OpenCRAVAT **2.17.0** and its `hg38` mapper **1.47.2**, which
uses GENCODE 47. Select the GENCODE release appropriate to your installed mapper;
the builder does not detect it automatically. Check the mapper metadata with
`oc module info hg38` and its release notes. These are pinned, reproducible
inputs; they are not live NCBI
lookups. The GENCODE release number is the cross-reference source, not a RefSeq
annotation release number.

Each accession pair is labeled as one of:

- **MANE Select / MANE Plus Clinical:** matched, full-versioned Ensembl and RefSeq
  transcripts. MANE defines matching sequence, exon coordinates and coding
  region for the pair on GRCh38.
- **GENCODE_RefSeq_xref:** an association supplied by GENCODE's cross-reference
  pipeline. This alone does **not** establish identical transcript sequence,
  exon structure or cDNA coordinates.

The annotator keeps all associated RNA accessions in the detailed output,
preserving versions. It never chooses one arbitrary RefSeq identifier from an
ambiguous association. For an
identical pair present in both sources, MANE takes precedence. Other RefSeq
associations remain visible with their own source labels. MANE pairs are included
even when their ENST version is absent from the older GENCODE metadata, but can
only match the exact ENST version supplied by the mapper.

Existing Ensembl-derived cDNA/protein changes are not copied into the new
columns or relabeled as RefSeq changes. For ordinary cross-references that would
be unreliable. Even a MANE match here is an identifier match, not independent
HGVS validation.

## Requirements and working directory

- Python **3.9+** for the builder and unit tests; no third-party Python packages
  are needed for these steps. Local verification used Python 3.10.
- An existing OpenCRAVAT installation with its base modules for annotation.
  Runtime integration has been tested with **OpenCRAVAT 2.17.0**. Other versions,
  including 3.x, have not been verified here.
- Network access to EMBL-EBI and NCBI for downloading references, or the local
  source files described below. Annotation itself makes no network requests.
- The shell examples target Linux/macOS with Bash. The smoke test uses directory
  symlinks; Windows support has not been verified.

Run every command below from this project's root, the directory containing
`build_annotator.py`. In the original parent repository, first run
`cd refseq_annotator`; in a detached clone, use the clone root. Activate the
OpenCRAVAT environment before commands using `oc` or `python smoke_test.py`.

## Build

```bash
python3 build_annotator.py --release 47 --mane-release 1.5
```

This downloads about 2 MB of compressed metadata into `.cache/` here and creates:

```text
build/refseq_transcripts/
  refseq_transcripts.py
  refseq_transcripts_engine.py
  refseq_transcripts.yml
  refseq_transcripts.md
  data/
    mappings.sqlite
    manifest.json
```

GENCODE downloads are checked against the publisher's MD5SUMS. The manifest
records URLs or local paths, source SHA256 hashes, database SHA256, releases,
build time and row counts. MANE's SHA256 is recorded for reproducibility, not
compared against a separate publisher digest. No patient data are downloaded or
used by the builder. Cache, generated modules and test artifacts are gitignored
within this directory. Allow roughly 50 MB of free disk space for one build.

The build refuses to overwrite an existing module and only publishes the output
after completing successfully. Rebuild into a fresh parent directory:

```bash
python3 build_annotator.py \
  --release 47 --mane-release 1.5 --output build/rebuilt
```

GENCODE **50** is the current release checked on 2026-09-28. After a mapper
upgrade, choose its corresponding release explicitly, for example `--release 50`.
Changing this metadata release does **not** update the mapper. Using a newer
cross-reference file with an older mapper may increase version-mismatch misses.
The annotator intentionally does not strip versions to hide these misses.

Offline build (gzip or plain text):

```bash
python3 build_annotator.py --release 47 --mane-release 1.5 \
  --refseq /path/gencode.v47.metadata.RefSeq.gz \
  --mane /path/MANE.GRCh38.v1.5.summary.txt.gz \
  --output build/offline
```

Use `--without-mane` instead of `--mane ...` to omit MANE; this leaves
`primary_refseq` empty while retaining ordinary cross-references in the other
outputs. The declared release labels for local files are your responsibility:
the TSV headers do not contain
release identifiers. Source hashes are still recorded. Only versioned RNA
accessions (`NM_`, `NR_`, `XM_`, `XR_`) are retained; protein accessions are not
transcript identifiers.

## Install and run

Activate your OpenCRAVAT environment and find its modules directory:

```bash
oc config md
```

Copy the generated module into `annotators/` below that path. Replace
`/path/to/oc/modules` with the printed path; the guard avoids overwriting an
installed copy:

```bash
modules_dir=/path/to/oc/modules
test ! -e "$modules_dir/annotators/refseq_transcripts" && \
  cp -R build/refseq_transcripts "$modules_dir/annotators/"
oc module ls -t annotator
```

That guard is for the **first installation only**: when the destination already
exists, it skips the copy. To update an installed module after rebuilding, back
it up and copy the new contents into the existing directory:

```bash
modules_dir=/path/to/oc/modules
backup_dir=".cache/installed-backups/refseq_transcripts-$(date +%Y%m%dT%H%M%S)"
mkdir -p "$(dirname "$backup_dir")"
cp -a "$modules_dir/annotators/refseq_transcripts" "$backup_dir" && \
  cp -R build/refseq_transcripts/. "$modules_dir/annotators/refseq_transcripts/"
oc module ls refseq_transcripts
```

The `/.` copies the directory contents. Restart an already running OpenCRAVAT
GUI to reload the module, and rerun annotation to update existing result files.

For a TSO500 VCF aligned to **GRCh37/hg19**:

```bash
oc run /path/to/tso500.vcf -l hg19 -a refseq_transcripts \
  -d test-output/tso500
```

Append `refseq_transcripts` to your usual annotator list to run the others as
well. Use `-l hg38` only for input already in GRCh38. OpenCRAVAT handles the
hg19-to-hg38 conversion before mapping. This module reads `.crx` transcript
identifiers after mapping and performs no coordinate conversion.

The OpenCRAVAT CLI also supports adding an annotator to an existing result:

```bash
oc run /path/to/existing.sqlite -a refseq_transcripts
```

That updates the specified result database; use a copy if you want to keep it
unchanged. Recreate any desired report separately with `oc report`. The builder
and tests here do not modify your existing result databases or installed modules.

## Output

OpenCRAVAT's mapper runs first and writes the selected `transcript` and
`all_mappings` to its `.crx` output. This annotator receives those fields for each
variant during the annotator stage. It follows the mapper's actual selection;
it does not choose the ENST from its own reference table.

All columns below have the prefix `refseq_transcripts__` in the result database.

| Column | Meaning |
| --- | --- |
| `primary_enst` | The existing mapper's selected Ensembl transcript, unchanged. |
| `primary_refseq` | One unique MANE RefSeq counterpart of the exact selected ENST version, or blank. |
| `primary_status` | Outcome of selecting the primary transcript's unique MANE counterpart; see below. |
| `all_refseq` | Unique RefSeq accessions for all affected mapper transcripts. |
| `matched_count` / `unmatched_count` | Counts of distinct affected ENST identifiers with/without exact-version matches. |
| `all_transcripts` | JSON rows with gene, original transcript, primary flag, status, and every RefSeq/source pair. |
| `gencode_release` / `mane_release` | Source metadata releases. |

A blank `primary_refseq` stays blank even when another transcript of the same
gene has a match. No gene-symbol fallback or silent primary-transcript change
occurs. The status is `mane_select` or `mane_plus_clinical` when there is exactly
one distinct MANE counterpart. Ordinary associations alone give `no_mane_match`,
even if only one association exists. Multiple distinct MANE counterparts give
`ambiguous_mane_match`; neither status fills `primary_refseq`. Other empty-result
statuses remain `no_exact_version_match`, `unmapped`, `unversioned_input`,
`not_ensembl`, and `no_transcript`.

`all_refseq`, the detailed pairs and transcript counts continue to include
ordinary associations. Detailed row statuses describe their strongest evidence
(including `xref_only`), independently of the stricter `primary_status`.

Version **1.1.0** introduces this single-counterpart policy. Rebuild and replace
the installed module to use it, then rerun the annotation on existing results;
previously generated result columns do not change automatically.
After building, the module is available at `build/refseq_transcripts/`. If you
use `--output`, adjust the install and smoke-test paths accordingly.

Example exact pair in MANE 1.5:

```text
BRAF: ENST00000646891.2 -> NM_004333.6 (MANE Select)
```

The lookup does not infer a mapping for `ENST00000646891.99999` or an unversioned
`ENST00000646891`. A missing match means this reference snapshot has no eligible
pair, not necessarily that NCBI has no corresponding transcript.

## Verification

Version 1.1.0 passed **26 regression tests** and the genuine OpenCRAVAT smoke
test on 2026-09-30. The smoke test now requires the primary BRAF output to equal
`NM_004333.6` exactly, while continuing to reject an incorrect ENST version.

Synthetic regression tests, without OpenCRAVAT or network access:

```bash
python3 -m unittest discover -s tests -v
```

After building, run the generated module through your installed OpenCRAVAT's
actual BaseAnnotator using an isolated modules directory:

```bash
python smoke_test.py
```

This uses synthetic `.crx` records, tests a real BRAF MANE match and a deliberate
version mismatch, and checks the output and error logs. It uses a process-local
`OPENCRAVAT_MD` override and stores diagnostic files under `test-output/` here.
It does not change global OpenCRAVAT configuration. The tests exercise the
identifier lookup; they are not validation of a new consequence mapper.

An earlier end-to-end check with OpenCRAVAT 2.17.0 converted synthetic hg19
BRAF `7:140453136 A>T` to hg38 `chr7:140753336 A>T` and associated the selected
`ENST00000646891.2` with `NM_004333.6`. This is a limited integration check, not
validation of variant consequences or coverage across a clinical panel.

## Exporting to a separate repository

Copy the source files only. Existing local directories such as `.cache/`,
`build/`, `build_old_v1.0.0/`, `test-output/`, and `__pycache__/` contain generated
or machine-specific material. Some old build copies were tracked in the parent
repository; `.gitignore` does not remove already tracked files. A source-only
archive avoids carrying them or the parent's Git history into a public repository:

```bash
tar --exclude='__pycache__' --exclude='*.pyc' -czf /tmp/refseq-annotator-source.tar.gz \
  README.md .gitignore annotation.py build_annotator.py smoke_test.py templates tests
```

Extract this archive into a new repository and run the unit tests there before
publishing. No built database is included; users generate it with the builder.
When a standalone license is added, include its file in the export as well.

### License decision before publication

This directory currently has **no standalone software license**. The original
parent repository contains a proprietary license. The maintainer needs to choose
and add the appropriate license for this separate project; this review does not
assign one. Reference data from GENCODE and NCBI have their own terms and citation
requirements, separate from any license chosen for the code.

## What an actual RefSeq mapper would require

A genuine mapper would start with RefSeq genomic alignments, exon/CDS models,
and transcript/protein sequences, then calculate variant effects in those
models. It would need defined primary-transcript selection, handling of
transcript/genome sequence differences, indels and splice boundaries, and
HGVS normalization. Existing OpenCRAVAT annotators that consume ENST identifiers
would also need compatibility review. That is a substantially larger task than
this lookup module. Merely replacing `base__transcript` with an NM accession
would not implement it.

## Documentation and sources checked

- [OpenCRAVAT annotator reference](https://docs.opencravat.org/en/latest/Annotator-Reference.html): module layout, `.crx` input, setup/annotate/cleanup.
- [OpenCRAVAT file formats](https://docs.opencravat.org/en/latest/File-Formats.html): `transcript` and the `all_mappings` JSON layout.
- [OpenCRAVAT command line](https://docs.opencravat.org/en/latest/2.-Command-line-usage.html): module installation, assembly choice, adding annotators to existing results.
- [GENCODE release 47](https://www.gencodegenes.org/human/release_47.html) and [release 47 metadata files](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/).
- [Current GENCODE release](https://www.gencodegenes.org/human/): metadata describes RefSeq associations from the Ensembl xref pipeline.
- [NCBI MANE](https://www.ncbi.nlm.nih.gov/refseq/MANE/): what a matched transcript pair guarantees.
- [MANE format README](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/README.txt) and [pinned release 1.5](https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/).
