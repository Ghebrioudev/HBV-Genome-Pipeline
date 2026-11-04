# Task 2 — Data preprocessing

Task 2 reads the standardized output of Task 1 and produces an HBV-focused,
validated dataset for later genomic analysis. It uses only the Python standard
library; no packages need to be installed.

## Run

From the repository root:

```bash
python3 preprocess.py
```

The default input directory is `data-collection/data/processed/`. Paths can be
overridden when needed:

```bash
python3 preprocess.py \
  --input-dir data-collection/data/processed \
  --output-dir data/processed \
  --reports-dir results/preprocessing
```

The script verifies the SHA-256 hashes recorded in the Task 1 `manifest.json`
before processing. A checksum mismatch, missing input or mismatch between the
manifest record count and the FASTA record count causes a non-zero exit. The Task
1 files are read-only and are never rewritten.

## Input contract

The required Task 1 files are:

| File | Role |
|---|---|
| `sequences.fasta` | Nucleotide sequences and FASTA headers; the first whitespace-delimited header token is the versioned accession |
| `metadata.csv` | NCBI metadata, including `accession`, `tax_id`, `lineage`, `length` and `completeness` |
| `manifest.json` | Task 1 record count, source details and SHA-256 hashes |

The metadata CSV and FASTA may contain the same accessions in different orders.
Task 2 therefore joins by exact `accession.version`; it never assumes that a
metadata row corresponds to the FASTA record at the same row number.

## Outputs

### Cleaned dataset — `data/processed/`

- `hbv_sequences.fasta` — accepted HBV nucleotide sequences, preserving FASTA
  headers and uppercasing sequence symbols; sequence lines are wrapped at 80 bases.
- `hbv_metadata.csv` — corresponding metadata in the same order as the FASTA,
  with added sequence length, SHA-256, ambiguity, duplicate-group and validation
  columns.
- `manifest.json` — input paths and hashes, output paths and hashes, record
  counts, processing version and configuration. The manifest does not hash
  itself (that would be self-referential).

### Audit reports — `results/preprocessing/`

- `cleaning_report.json` — record counts, sequence-length statistics, ambiguity
  statistics, completeness counts and the rules applied.
- `rejected_records.csv` — invalid sequences, conflicting duplicate accession
  records, missing-accession records and metadata rows without a FASTA record.
- `excluded_non_hbv.csv` — records whose NCBI lineage places them in
  Hepadnaviridae but outside HBV.
- `unresolved_taxonomy.csv` — records without sufficient or internally
  consistent NCBI taxonomy information, including FASTA records with no metadata.
- `review_records.csv` — records containing U, alignment gaps, or a mismatch
  between metadata length and observed FASTA length.
- `duplicate_accessions.csv` — duplicate accession occurrences and the action
  taken (identical copy removed or conflicting records quarantined).
- `identical_sequences.csv` — one row per accession in each group of identical
  nucleotide strings. Different accessions are all preserved in the cleaned set.

Empty reports contain their header row, so the audit outputs are always present.

## Scientific rules

### Taxonomic scope

A record is included when its NCBI `tax_id` is 10407, or its NCBI lineage
contains the exact `Hepatitis B virus` taxon (which also recognises a descendant
whose lineage contains that ancestor). A lineage containing `Hepadnaviridae`
but not `Hepatitis B virus` is excluded from the primary HBV dataset. Missing,
incomplete or contradictory lineage information is reported for review rather
than guessed from organism-name text alone.

### Sequence alphabet and ambiguity

The accepted DNA alphabet is `ACGT` plus the IUPAC ambiguity symbols
`RYSWKMBDHVN`. Ambiguity symbols are retained and counted; no base is imputed.
A mutation analysis should not make a definitive base call at an ambiguous
position without an explicit downstream policy.

`U` and `-` are not silently converted or deleted. Records containing U
(uracil) or alignment gaps are placed in `review_records.csv` and are not passed
to the cleaned FASTA until reviewed. Other unsupported symbols cause rejection
with the offending symbols recorded in `rejected_records.csv`.

### Partial sequences

`COMPLETE`, `PARTIAL` and missing/unknown completeness classifications are
preserved when the sequence passes validation. No minimum sequence length is
applied. Task 3 must choose whether an individual sequence is suitable for a
whole-genome or region-specific comparison; absence of sequence coverage must
not be interpreted as a deletion.

### Duplicate records

Duplicate accession detection uses the full versioned accession string. If the
same accession.version occurs more than once with identical normalized sequence
content, the first record is retained and subsequent copies are audited. If one
accession.version has conflicting sequence contents, all records for that
accession.version are quarantined. Different accession versions are never
merged automatically.

Identical sequence detection is separate. SHA-256 is used to find candidate
groups and sequence strings are compared exactly before reporting a group.
Records with different accessions remain in the cleaned dataset even when their
sequence strings are identical. Reverse complements are not treated as identical.

### Alignment

Task 2 does not perform sequence alignment. Alignment and mutation calling are
left to Task 3, which must account for sequence coverage and ambiguity.

## Current input snapshot included in this archive

The included `data-collection/data/raw/retrieval.json` records an NCBI request
for taxon 10407, and the Task 1 processed manifest reports 140,855 records.
The `data-collection/dataset_documentation.md` and `SOURCES.md` still describe an
earlier family-wide (taxon 10404) snapshot. For the actual files in this archive,
`retrieval.json` and `data-collection/data/processed/manifest.json` are the
machine-readable provenance for the current snapshot. Task 2 checks the
processed FASTA and metadata hashes against that manifest.

### Results for the bundled snapshot

The production run on this archive retained all **140,855** input records and
**129,365,332** nucleotide positions. There were no rejected records, non-HBV
exclusions, unresolved taxonomy records, review records, or duplicate
accession.version groups. FASTA and metadata accessions matched as sets, but
were not in the same source order; the cleaned output is joined by accession.

- **Completeness:** 127,395 `PARTIAL`; 13,460 `COMPLETE`. Both classes are
  retained. No minimum sequence length is imposed.
- **Length:** minimum 1 nt, median 606 nt, mean 918.4291 nt, maximum 10,123 nt.
  In particular, the shortest records should not be assumed suitable for
  whole-genome comparisons; Task 3 must choose suitable genomic coverage.
- **Ambiguity:** 508,032 valid IUPAC ambiguity positions across 18,268 records
  (0.392711% of retained positions). These symbols remain in the FASTA.
- **Exact sequence identity:** 11,912 groups contain 43,200 records. Different
  accessions in these groups are all retained; the 31,288 additional group
  members are not duplicate accession records.
- **Accession duplicates, invalid symbols, U and gaps:** none were observed in
  this snapshot. The corresponding reports are still generated for auditability.

These are measurements of the dataset bundled in this archive, not fixed
expectations for future NCBI downloads.

## Tests

Run Task 2 tests:

```bash
python3 -m unittest discover -s tests/preprocessing -v
```

Run the original Task 1 suite:

```bash
bash data-collection/tests/run_tests.sh
```

Run both suites using the repository-level test runner:

```bash
bash tests/run_tests.sh
```

All tests use `unittest` and temporary fixtures. They do not access NCBI or
modify the real Task 1 dataset.
