# HBV-Genome-Pipeline

A reproducible bioinformatics pipeline for analysing **Hepatitis B Virus (HBV)**
genomic data using Shell and Python.

## Status

| Phase | Scope | State |
|---|---|---|
| 1 | Data collection (NCBI API + local file) | ✅ implemented |
| 2 | Preprocessing | ✅ implemented |
| 3 | Sequence comparison / mutation detection | not started |
| 4 | Statistics and figures | not started |
| 5 | Report | not started |

## Requirements

- Python ≥ 3.8 (tested on 3.12.3) — **standard library only**
- bash ≥ 4

No third-party Python packages are required for Phases 1 or 2.

## Quick start

```bash
# Download the HBV/Hepadnaviridae dataset from NCBI (default source)
./download.sh

# Hepatitis B virus only
./download.sh --taxon 10407

# Use a local FASTA / .gz / .zip instead
./download.sh --file path/to/hbv.fasta
```

Outputs land in `data-collection/data/`, and the full provenance is documented in
[`data-collection/dataset_documentation.md`](data-collection/dataset_documentation.md).

## Task 2 — Preprocessing

Run the preprocessing stage from the repository root:

```bash
python3 preprocess.py
```

Task 2 reads the three files in `data-collection/data/processed/`, verifies the
Task 1 manifest hashes, and writes the cleaned HBV dataset to `data/processed/`
and audit reports to `results/preprocessing/`. It preserves partial sequences
and IUPAC ambiguity symbols, audits accession and sequence duplicates, and does
not perform alignment. See [`PREPROCESSING.md`](PREPROCESSING.md) for the rules,
outputs, and command-line options.

## Tests

```bash
bash tests/run_tests.sh                 # Task 1 + Task 2 suites, offline
bash data-collection/tests/run_tests.sh # Task 1 suite only
python3 -m unittest discover -s tests/preprocessing -v # Task 2 only
```

No `pytest` or third-party Python package is required. Tests use temporary
fixtures and do not contact NCBI or modify the actual Task 1 dataset.

## Repository layout

```
download.sh                  entry point (delegates to data-collection/)
data-collection/            Phase 1 — data acquisition
  download.sh               orchestrator (NCBI default, --file alternative)
  dataset_documentation.md  provenance, scope, checksums
  README.md                 detailed setup and troubleshooting
  scripts/validate_dataset.sh
  src/ncbi_client.py        NCBI Datasets REST client
  src/file_loader.py        local + NCBI normalisation
  src/dataset_validator.py  integrity / format / size checks
  data/{raw,processed,checksums}/
preprocess.py                Phase 2 — validated HBV preprocessing
analyze_mutations.py         Phase 3 (empty)
statistics.py                Phase 4 (empty)
report/report.md             Phase 5 (not started)
tests/preprocessing/         Phase 2 offline unittest suite
PREPROCESSING.md             Phase 2 policies and output documentation
```

See [`AGENTS.md`](AGENTS.md) for conventions used by contributors and AI agents.

## Dataset

The current Task 1 files bundled with this project were retrieved for taxon
`10407` (Hepatitis B virus) and contain **140,855 sequences · 129,365,332 nt**.
The machine-readable provenance is in
`data-collection/data/raw/retrieval.json` and
`data-collection/data/processed/manifest.json`. Some historical Task 1
documentation files still describe an earlier family-wide taxon `10404` snapshot;
see [`PREPROCESSING.md`](PREPROCESSING.md) for this discrepancy.
