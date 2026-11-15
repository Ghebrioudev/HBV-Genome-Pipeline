# Dataset documentation — HBV nucleotide snapshot (Task 1)

This document describes the exact NCBI dataset snapshot bundled with the
project. The machine-readable provenance is recorded in:

- `data/raw/retrieval.json` — request parameters, retrieval time, archive size
  and archive SHA-256.
- `data/checksums/ncbi_validation.json` — validation measurements for the
  extracted NCBI package.
- `data/processed/manifest.json` — Task 1 processed record count and file
  checksums.
- `data/checksums/dataset.sha256` and `data/checksums/ncbi_extracted.sha256` —
  SHA-256 checksums for the raw package, extracted files and normalised outputs.

## 1. Dataset scope and retrieval

The raw NCBI request was made on **14 November 2026 at 12:36:39 UTC** through
the NCBI Datasets Virus v2 API. It requested nucleotide genome records for
*Hepadnaviridae* (taxon **10404**), the family that includes HBV and related
hepadnaviruses.

| Field | Value |
|---|---|
| Service | NCBI Datasets Virus v2 |
| Endpoint | `POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download` |
| Taxon requested | 10404 — *Hepadnaviridae* |
| Sequence type | Genome nucleotide sequences |
| Retrieval time (UTC) | 2026-11-14 12:36:39 |
| Input records | 142,288 |
| Total nucleotide positions | 132,005,421 |
| Raw ZIP size | 32,898,096 bytes |
| Extracted package size | 425,368,669 bytes (405.66 MiB) |
| ZIP SHA-256 | `968ce5165734b61df7164451d57de429c0ac22317997b3d18e88ff9d42bd324b` |

The downloaded package was validated successfully (`ok: true`) with no recorded
warnings. Its `genomic.fna` file contains 142,288 records and 132,005,421
nucleotide positions. The accompanying `data_report.jsonl` contains the NCBI
record metadata used to create the normalised metadata CSV.

## 2. HBV filtering and Task 2 dataset

The NCBI request is family-level, not HBV-only. Task 2 uses the NCBI taxonomic
identifier and lineage fields to include records assigned to taxon 10407 or
whose lineage contains the exact *Hepatitis B virus* taxon. It excluded **1,399
non-HBV hepadnavirus records** and retained **140,889 HBV records**.

| Task 2 measurement | Value |
|---|---:|
| Input FASTA and metadata records | 142,288 |
| HBV records retained | 140,889 |
| Non-HBV records excluded | 1,399 |
| Retained nucleotide positions | 129,414,793 |
| `COMPLETE` records retained | 13,469 |
| `PARTIAL` records retained | 127,420 |
| Sequence length, minimum / median / mean / maximum | 1 / 606 / 918.5585 / 10,123 nt |

Task 2 found no unresolved taxonomy records, rejected records, records requiring
manual review, or duplicate accession.version groups in this snapshot. It
preserved partial sequences and valid IUPAC ambiguity symbols.

### Ambiguity and identical sequences

There are **508,170** valid IUPAC ambiguity positions across **18,285** retained
records, representing **0.392668%** of all retained nucleotide positions. The
symbols are preserved rather than replaced with inferred bases.

Exact sequence comparison identified **11,913 groups** containing **43,203
records**. The **31,290 additional group members** are retained because they
have different accessions. Identical sequence strings are not treated as
identical accession records, and reverse complements are not grouped together.

## 3. Files and formats

| File | Description |
|---|---|
| `data/raw/ncbi_dataset.zip` | Original NCBI data package downloaded from the API |
| `data/raw/retrieval.json` | Request, retrieval time and package checksum |
| `data/processed/sequences.fasta` | Task 1 normalised nucleotide sequences |
| `data/processed/metadata.csv` | Task 1 normalised NCBI metadata |
| `data/processed/manifest.json` | Task 1 record count and output checksums |
| `data/checksums/ncbi_validation.json` | Validation summary for the extracted package |
| `data/checksums/ncbi_extracted.sha256` | Checksums for files extracted from the raw ZIP |
| `data/checksums/dataset.sha256` | Checksums for the raw ZIP and Task 1 processed files |

The expanded `data/raw/extracted/` directory is created by the download stage
when needed, but is not duplicated in this submission archive: the original
NCBI ZIP is included and the client re-extracts it before validation and
normalisation. This avoids packaging a second copy of hundreds of megabytes of
source data.

Task 1 processed file checksums:

| File | Bytes | SHA-256 |
|---|---:|---|
| `sequences.fasta` | 146,731,669 | `b02f7a2ebf188467b92ef2be03e40556d458bed05b655c3db8aeb0c1c9f3b3d1` |
| `metadata.csv` | 53,165,248 | `e90499ba688d6ca37f00617ec68151425b1b71242a9f25686da9194b6aa5e202` |

Task 2 output checksums are recorded in `data/processed/manifest.json`.

## 4. Metadata fields and limitations

The Task 1 metadata includes accession, sequence length, completeness,
annotation flags, gene/protein counts, source database, release/update dates,
organism name, taxonomic identifier, isolate, lineage and submitter where
available. Missing NCBI fields are left blank. There is no consistently
populated, standardised genotype field in the supplied metadata, so the
pipeline does not produce a genotype comparison.

The dataset is a time-specific NCBI snapshot. Future downloads may differ in
record count, sequence content and metadata. The family-level request includes
related hepadnaviruses; downstream HBV analysis must continue to apply Task 2's
HBV taxonomic filter.

## 5. Reproduction

From the repository root:

```bash
# Retrieve the default family-level dataset (taxon 10404)
./download.sh

# Retrieve HBV-only records directly (taxon 10407)
./download.sh --taxon 10407

# Use a local FASTA instead of the NCBI API
./download.sh --file path/to/hbv.fasta
```

The default pipeline uses the family-level request and filters to HBV in Task 2.
To use the already bundled Task 1 snapshot without making a new network request,
run the later stages with `./run.sh --skip-download`.
