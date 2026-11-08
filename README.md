# HBV-Genome-Pipeline

A reproducible bioinformatics pipeline for analysing **Hepatitis B Virus (HBV)**
genomic data using Shell and Python.

## Status

| Phase | Scope | State |
|---|---|---|
| 1 | Data collection (NCBI API + local file) | ✅ implemented |
| 2 | Preprocessing | ✅ implemented |
| 3 | Sequence comparison / mutation detection | ✅ implemented |
| 4 | Statistics and figures | ✅ implemented |
| 5 | Report | ✅ implemented |

## Requirements

- Python ≥ 3.10 (tested on 3.12.3)
- bash ≥ 4
- MAFFT ≥ 7.487 for Task 3 (system executable)
- pandas and matplotlib for Task 4

Task 1–3 use the Python standard library; Task 4 additionally requires pandas and
matplotlib.

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

## Task 3 — Sequence comparison and mutation detection

After Task 2, install MAFFT (`sudo apt update && sudo apt install mafft`) and run:

```bash
mafft --version
python3 analyze_mutations.py
```

Task 3 aligns eligible samples to the included `NC_003977.2` RefSeq using
MAFFT `--6merpair --addfragments` (500 sequences per batch by default, with
automatic batch splitting after an alignment failure), then writes
reference-coordinate substitutions,
insertions, deletions and uncertain ambiguous-base calls to
`results/mutations.csv`. Skipped sequences and run provenance are recorded in
`results/mutation_analysis_skipped.csv` and `results/mutation_analysis.json`.
See [`TASK3.md`](TASK3.md) for the alignment policy, thresholds, output schema,
known limitations and options.

## Tests

```bash
bash tests/run_tests.sh                 # Task 1 + Task 2 + Task 3 + Task 4 suites, offline
bash data-collection/tests/run_tests.sh # Task 1 suite only
python3 -m unittest discover -s tests/preprocessing -v # Task 2 only
python3 -m unittest discover -s tests/task3 -v         # Task 3 only
python3 -m unittest discover -s tests/task4 -v         # Task 4 only
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
analyze_mutations.py         Phase 3 — reference-aligned mutation detection
statistics.py                Phase 4 — statistics and figures
report/report.md             Phase 5 (not started)
tests/preprocessing/         Phase 2 offline unittest suite
tests/task3/                  Phase 3 offline unittest suite
TASK3.md                     Phase 3 method, usage and limitations
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

## Running the pipeline

Run all stages in order from the repository root:

```bash
./run.sh
```

For an offline integration run when input data has already been collected:

```bash
./run.sh --clean --skip-download
```

The `--clean` option removes generated CSV files and figures before running.
The orchestrator uses `set -euo pipefail`, checks required stage files, and
stops immediately when a stage fails. Task 3 requires MAFFT, and Task 4 requires
pandas and matplotlib.

## Integration test

The offline integration test runs the pipeline from a temporary clean copy and
checks that invalid command-line options fail:

```bash
bash tests/run_integration_test.sh
```
