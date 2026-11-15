# HBV Genome Pipeline

A reproducible, multi-stage bioinformatics workflow for collecting HBV nucleotide
records, validating and preprocessing sequences, comparing them with a reference
genome, and summarising observed nucleotide variation.

## Project status

| Stage | Scope | Status |
|---|---|---|
| 1 | NCBI data collection and normalisation | Implemented |
| 2 | HBV taxonomy filtering, sequence validation and audit reports | Implemented |
| 3 | Reference-anchored alignment and mutation detection | Implemented |
| 4 | Mutation statistics and figures | Implemented |
| 5 | Results report and self-contained HTML dashboard | Implemented |

## Requirements

- Python **3.11 or newer** for the complete project.
- Bash 4 or newer.
- MAFFT **7.487 or newer** for Task 3.
- The Python packages pinned in `requirements.txt` (`pandas` and `matplotlib`) for
  Task 4 and the complete test suite.

Tasks 1–3 use the Python standard library. MAFFT is a separate system executable,
not a Python package.

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r requirements.txt
```

Install MAFFT on Ubuntu/Debian and confirm that it is available on `PATH`:

```bash
sudo apt update
sudo apt install mafft
mafft --version
```

## Run the pipeline

From the repository root:

```bash
./run.sh
```

By default, Task 1 retrieves the NCBI *Hepadnaviridae* family dataset (taxon
10404); Task 2 then retains records classified as Hepatitis B virus. The
repository includes a fixed input snapshot, so a rerun using those bundled Task
1 files can skip the network download:

```bash
./run.sh --skip-download
```

To retrieve HBV records directly from NCBI instead of the family-level dataset:

```bash
./download.sh --taxon 10407
```

Task 1 also accepts a local FASTA, gzip-compressed FASTA or ZIP containing FASTA:

```bash
./download.sh --file path/to/sequences.fasta
```

`run.sh --clean` removes the main mutation/summary CSVs and figure files before
running; it does not remove the bundled input data. The orchestrator stops when
a stage exits unsuccessfully.

## Pipeline stages

### Task 1 — Data collection

`download.sh` delegates to `data-collection/download.sh`. The default NCBI
pathway uses the NCBI Datasets Virus REST API; a local-file pathway is also
available. The stage validates the source package and writes normalised
`sequences.fasta`, `metadata.csv` and `manifest.json` files under
`data-collection/data/processed/`. Provenance and SHA-256 checksums are recorded
in `data-collection/data/raw/retrieval.json` and `data-collection/data/checksums/`.
See [`data-collection/README.md`](data-collection/README.md) and
[`data-collection/dataset_documentation.md`](data-collection/dataset_documentation.md).

### Task 2 — Preprocessing

```bash
python3 preprocess.py
```

Task 2 verifies the Task 1 manifest, joins FASTA and metadata by versioned
accession, filters to HBV taxonomy, and produces the cleaned data under
`data/processed/` plus audit reports under `results/preprocessing/`. It retains
valid partial sequences and IUPAC ambiguity symbols, distinguishes duplicate
accessions from identical sequences, and does not perform alignment. See
[`PREPROCESSING.md`](PREPROCESSING.md).

### Task 3 — Sequence comparison and mutation detection

```bash
python3 analyze_mutations.py
```

Task 3 uses MAFFT `--6merpair --addfragments` to align eligible sample fragments
to the included `NC_003977.2` reference in batches. It calls substitutions,
insertions, deletions and uncertain ambiguous-base events in 1-based reference
coordinates. Partial terminal coverage is not called as a deletion. Run
thresholds, output columns and known limitations are described in
[`TASK3.md`](TASK3.md).

### Task 4 — Statistics and figures

```bash
python3 statistics.py
```

The stage streams the large mutation CSV in chunks and writes
`results/summary.csv`, `results/mutation_types.csv`,
`results/statistics_summary.json` and figures under `results/figures/`.
Uncertain IUPAC calls remain in the mutation table but are excluded from
confirmed-mutation frequency statistics. Genotype comparison is not produced
because the supplied metadata does not include a usable genotype column.

### Task 5 — Report and dashboard

After Tasks 1–4 have completed, generate the self-contained HTML dashboard:

```bash
python3 report/make_dashboard.py
```

The dashboard is written to `report/dashboard.html`. The written project
summary is also available in [`report/report.md`](report/report.md), and the
submission report is provided as `report/final_project_report.docx`.

## Tests

Install the requirements first, then run the complete offline suite:

```bash
bash tests/run_tests.sh
```

The runner covers Task 1, Task 2, Task 3, Task 4 and pipeline integration. It
uses temporary fixtures and does not contact NCBI or overwrite the bundled
production dataset. The real-MAFFT integration test is skipped when MAFFT is
not installed; Task 3's other tests use a deterministic MAFFT test double.

Individual suites:

```bash
bash data-collection/tests/run_tests.sh
python3 -m unittest discover -s tests/preprocessing -p 'test_*.py' -v
python3 -m unittest discover -s tests/task3 -p 'test_*.py' -v
python3 -m unittest discover -s tests/task4 -p 'test_*.py' -v
python3 -m unittest discover -s tests/integration -p 'test_*.py' -v
```

## Bundled dataset and completed run

The included NCBI snapshot was retrieved on **14 November 2026** using taxon
10404 (*Hepadnaviridae*): 142,288 input records and 132,005,421 nucleotide
positions. Task 2 excluded 1,399 non-HBV records and retained **140,889 HBV
records** (129,414,793 nucleotide positions). Of these, 127,420 are marked
`PARTIAL` and 13,469 `COMPLETE`; valid IUPAC ambiguity symbols are preserved.

The included Task 3 run used MAFFT v7.525 and analysed **128,406 sequences**.
It produced **10,567,649 mutation events**, of which 10,098,062 are definitive
and 469,587 are uncertain ambiguous-base calls. Task 4 reports 16,235 distinct
confirmed mutations. These values describe the bundled snapshot and its
recorded run; they are not fixed expectations for a future NCBI retrieval.

## Repository layout

```text
README.md
PREPROCESSING.md                 Task 2 rules and outputs
TASK3.md                          Task 3 method and limitations
SOURCES.md                        external data/tool provenance
download.sh                       Task 1 entry point
run.sh                            stage orchestrator
preprocess.py                     Task 2 implementation
analyze_mutations.py              Task 3 implementation
statistics.py                     Task 4 implementation
data-collection/                  Task 1 source, tests and data
  data/processed/                 bundled Task 1 input snapshot
data/processed/                   cleaned Task 2 dataset
results/                           audit reports, mutation table and figures
report/                             Markdown report and HTML dashboard generator
  final_project_report.docx       four-page submission report
  report.md                       detailed results summary
  dashboard.html                  self-contained dashboard
  make_dashboard.py               dashboard generator
tests/                              Task 2–4 and integration tests
```
