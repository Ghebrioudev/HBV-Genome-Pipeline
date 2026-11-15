# Task 2 — Data preprocessing

Task 2 reads the standardised output of Task 1 and produces a validated,
HBV-focused dataset for downstream genomic analysis. It uses only the Python
standard library.

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

Before processing, the script verifies the SHA-256 hashes recorded in the Task 1
`manifest.json`. A checksum mismatch, missing input or mismatch between the
manifest record count and FASTA count causes a non-zero exit. Task 1 files are
read-only and are never rewritten.

## Input contract

| File | Role |
|---|---|
| `sequences.fasta` | Nucleotide sequences; the first whitespace-delimited FASTA header token is the versioned accession |
| `metadata.csv` | NCBI metadata, including accession, taxonomic information, length and completeness |
| `manifest.json` | Task 1 record count, source details and SHA-256 hashes |

The metadata and FASTA may contain the same accessions in different orders.
Task 2 therefore joins by exact `accession.version`; it does not assume matching
row order.

## Outputs

### Cleaned dataset — `data/processed/`

- `hbv_sequences.fasta` — accepted HBV sequences, preserving FASTA headers and
  uppercasing symbols; sequence lines are wrapped at 80 bases.
- `hbv_metadata.csv` — corresponding metadata in FASTA order, with additional
  sequence length, SHA-256, ambiguity, duplicate-group and validation columns.
- `manifest.json` — input/output checksums, record counts, processing version
  and configuration. The manifest does not hash itself.

### Audit reports — `results/preprocessing/`

- `cleaning_report.json` — record counts, length and ambiguity statistics,
  completeness counts and applied rules.
- `rejected_records.csv` — invalid sequences, conflicting duplicate accessions,
  missing-accession records and metadata rows without a FASTA record.
- `excluded_non_hbv.csv` — records classified within Hepadnaviridae but outside
  HBV.
- `unresolved_taxonomy.csv` — records with insufficient or inconsistent
  taxonomy evidence.
- `review_records.csv` — records containing U, alignment gaps or a discrepancy
  between metadata length and observed FASTA length.
- `duplicate_accessions.csv` — duplicate versioned accession occurrences and
  the action taken.
- `identical_sequences.csv` — accessions grouped by exactly identical
  nucleotide strings. Different accessions remain in the cleaned dataset.

Empty audit reports contain their header row so the report set is consistent.

## Scientific rules

### Taxonomic scope

A record is included when `tax_id` is 10407 or its NCBI lineage contains the
exact *Hepatitis B virus* taxon. A record whose lineage identifies another
Hepadnaviridae member is excluded. Missing, incomplete or contradictory lineage
information is reported for review rather than inferred from organism-name
text alone.

### Sequence alphabet and ambiguity

The accepted DNA alphabet is `ACGT` plus the IUPAC ambiguity symbols
`RYSWKMBDHVN`. Valid ambiguity symbols are preserved and counted; no base is
imputed. `U` and `-` are not silently converted or deleted: records containing
them are placed in `review_records.csv`. Other unsupported symbols cause
rejection and are listed in `rejected_records.csv`.

### Partial sequences

`COMPLETE`, `PARTIAL` and missing/unknown completeness classifications are
preserved when the record passes validation. No minimum sequence length is
applied in Task 2. Task 3 applies its own explicit eligibility rules. Unobserved
terminal sequence must not be interpreted as a deletion.

### Duplicate records

Duplicate accession detection uses the full versioned accession. Identical
copies of one accession.version retain the first record; conflicting sequences
for the same accession.version are quarantined. Different accession versions
are not merged. Separately, identical sequence candidates are grouped by
SHA-256 and then confirmed by exact sequence-string comparison. Different
accessions with identical sequence content remain in the cleaned dataset.
Reverse complements are not treated as identical.

### Alignment

Task 2 does not align sequences. Alignment and mutation calling are performed
in Task 3, which accounts for sequence coverage and ambiguity.

## Bundled snapshot and measured results

The included Task 1 snapshot was retrieved on **14 November 2026** from NCBI
using the family-level taxon 10404 (*Hepadnaviridae*). Its manifest contains
**142,288 records** and 132,005,421 nucleotide positions. The Task 2 run
retained **140,889 HBV records** (129,414,793 nucleotide positions) and
excluded **1,399 non-HBV hepadnaviruses**. No records were unresolved, rejected
or marked for manual review; no duplicate accession.version groups were found.

- **Completeness:** 127,420 `PARTIAL` and 13,469 `COMPLETE`; both are retained.
- **Length:** minimum 1 nt, median 606 nt, mean 918.5585 nt and maximum
  10,123 nt.
- **Ambiguity:** 508,170 valid IUPAC ambiguity positions across 18,285 records
  (0.392668% of retained nucleotide positions). These symbols remain unchanged.
- **Identical sequences:** 11,913 groups contain 43,203 records. The 31,290
  additional group members are retained because they have different accessions.
- **Accession duplicates:** no duplicate accession.version groups were found.

These figures are measured from the files included in this archive. They are
snapshot-specific and may change after a new NCBI retrieval. The detailed
provenance and checksums are in
[`data-collection/dataset_documentation.md`](data-collection/dataset_documentation.md).

## Tests

```bash
python3 -m unittest discover -s tests/preprocessing -p 'test_*.py' -v
bash data-collection/tests/run_tests.sh
bash tests/run_tests.sh
```

The tests use `unittest` and temporary fixtures. They do not contact NCBI or
modify the bundled production dataset.
