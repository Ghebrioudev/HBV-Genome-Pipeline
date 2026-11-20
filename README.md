# HBV Genome Pipeline — Projet « Enquête sur une maladie »

A reproducible pipeline that downloads every public *Hepatitis B virus* (HBV)
nucleotide sequence from NCBI Virus, filters and validates them, aligns them to
the reference genome **NC_003977.2**, detects nucleotide mutations and
summarises how often each mutation occurs. One command, `./run.sh`, goes from
the raw data to figures, tables, an HTML dashboard and a dated archive.

| | |
|---|---|
| **Disease** | Hepatitis B (chronic HBV infection) |
| **Research question** | Which nucleotide positions of the HBV genome are most frequently mutated, relative to the reference genome NC_003977.2, among the HBV sequences publicly available in NCBI Virus? |
| **Dataset** | NCBI Datasets *virus genome* package, taxon 10404 (*Hepadnaviridae*), ≈ 406 MB uncompressed — see [SOURCES.md](SOURCES.md) |
| **Pattern (motif)** | A *mutation* = a definitive difference from NC_003977.2 at a 1-based reference position: `position:REF>ALT` (substitution, insertion or deletion). A *recurrent mutation* is one carried by many analysed sequences; its frequency is `n_sequences / analysed sequences`. |
| **Report** | [rapport.pdf](rapport.pdf) (4 pages) — source in [report/rapport.md](report/rapport.md) |

---

## 1. Quick start (clean Ubuntu/Debian machine)

```bash
# 1. system packages: Python >= 3.11 with venv, and MAFFT >= 7.487
sudo apt update
sudo apt install -y python3 python3-venv mafft

# 2. isolated Python environment (.venv) from requirements.txt
bash scripts/setup_env.sh

# 3. run the whole pipeline (downloads ~33 MB, writes ~1.5 GB, ~45 min)
./run.sh
```

That is all. `run.sh` automatically uses `.venv/bin/python`; there is no need
to activate the environment. Re-running `./run.sh` at any time (even after an
interruption with Ctrl-C) works without any manual cleanup.

Tested on a clean `ubuntu:24.04` container (Python 3.12.3, MAFFT 7.505).

### Requirements

| Component | Version | Installed by |
|---|---|---|
| Bash, coreutils (`sha256sum`, `md5sum`, `sort`, `cut`, `join`…), `grep`, `sed`, `awk`, `tar`, `gzip` | Bash ≥ 4.4 | present on every Linux distribution |
| Python | ≥ 3.11 | `apt install python3 python3-venv` |
| pandas, matplotlib | pinned in `requirements.txt` | `bash scripts/setup_env.sh` |
| MAFFT | ≥ 7.487 | `apt install mafft` |
| Network access to `api.ncbi.nlm.nih.gov` | — | only for the download stage |

Disk space: about 2 GB free (raw package 406 MB, cleaned data 213 MB,
mutation table ≈ 800 MB). Tasks 1–3 use only the Python standard library.

## 2. Usage

```text
./run.sh [-h] [-s] [-c] [-f] [-n] [-t TAXON] [-j THREADS] [-k TOP]
```

| Option | Meaning |
|---|---|
| `-h` | show the help and exit |
| `-s` | skip the download and reuse the data in `data-collection/data/` (checksums are still verified) |
| `-c` | delete previously generated results before running |
| `-f` | force a fresh NCBI download even if a valid cached package exists |
| `-n` | do not create the dated tar archive |
| `-t TAXON` | NCBI taxonomy id to download (default `10404`, *Hepadnaviridae*; HBV records are selected afterwards) |
| `-j THREADS` | MAFFT threads per alignment batch (default `1`, as in the reference run) |
| `-k TOP` | number of top mutations plotted and cross-checked in shell (default `20`) |

Options are parsed with `getopts`; the long aliases `--help`,
`--skip-download`, `--clean` and `--force` are also accepted. Environment
variables: `PYTHON` (interpreter), `NCBI_API_KEY` (optional NCBI key).

**Exit codes:** `0` success · `1` a stage failed · `2` bad option · `3` missing
file · `4` missing dependency · `5` checksum mismatch · `130`/`143` interrupted
(Ctrl-C / SIGTERM). Every error is printed with a timestamp and an explicit
message, e.g. `[ERROR] required file does not exist: …/dataset.sha256 — run without -s to download the dataset`.

Examples:

```bash
./run.sh                 # full run, download if needed
./run.sh -s -j 4         # reuse the downloaded data, 4 MAFFT threads
./run.sh -c -f           # wipe results, fresh NCBI snapshot
```

## 3. What `run.sh` does

| # | Stage | Implementation | Main outputs |
|---|---|---|---|
| 1 | Download the NCBI package (cached; `-f` to refresh) and normalise it | `data-collection/download.sh` → `src/ncbi_client.py`, `dataset_validator.py`, `file_loader.py` | `data-collection/data/raw/`, `data-collection/data/processed/`, SHA-256 lists in `data-collection/data/checksums/` |
| 2 | **Integrity**: NCBI `md5sum.txt` checked with `md5sum -c`, project lists with `sha256sum -c`; raw data set **read-only** | `verify_md5_file`, `verify_sha256_file`, `set_readonly` in `scripts/lib.sh` | log lines `OK: n files match` |
| 3 | **Shell profiling** of the sequences with pipes (`grep`, `cut`, `sed`, `sort`, `uniq`, `awk`, `tr`) | `scripts/shell_analysis.sh -m sequences` | `results/shell/sequence_profile.tsv`, `length_distribution.tsv`, `header_annotations.tsv`, `base_composition.tsv` |
| 4 | HBV taxonomy filter, sequence validation, duplicate audit | `preprocess.py` (argparse) | `data/processed/hbv_sequences.fasta`, `hbv_metadata.csv`, `manifest.json`; `results/preprocessing/*.csv`, `cleaning_report.json` |
| 5 | MAFFT `--6merpair --addfragments` alignment to NC_003977.2 and mutation calling | `analyze_mutations.py` (argparse) | `results/mutations.csv`, `mutation_analysis_skipped.csv`, `mutation_analysis.json` |
| 6 | Statistics and figures | `statistics.py` (argparse, pandas, matplotlib) | `results/summary.csv`, `mutation_types.csv`, `statistics_summary.json`, `results/figures/*.png` |
| 7 | **Shell cross-check**: the top mutations are recounted with a pipe and compared with Python | `scripts/shell_analysis.sh -m mutations` | `results/shell/top_mutations_shell.tsv`, `crosscheck_top_mutations.tsv`, `events_by_type.tsv` |
| 8 | Self-contained HTML dashboard | `report/make_dashboard.py` | `report/dashboard.html` |
| 9 | SHA-256 of every result, then verified | `write_sha256_manifest` + `verify_sha256_file` | `results/checksums.sha256` |
| 10 | Dated archive of the results, listed back with `tar -tzf` | `archive_results` | `archives/hbv_results_YYYYMMDD_HHMMSS.tar.gz` |

Each run writes a **timestamped log** `logs/run_YYYYMMDD_HHMMSS.log` (every
line prefixed by `[YYYY-MM-DD HH:MM:SS]`, with stage durations).

### Robustness

* `set -euo pipefail` everywhere; an `ERR` trap prints the failing line and command.
* `trap` on `EXIT`, `INT` and `TERM`: the running stage (and MAFFT) is stopped and
  all temporary files are deleted — the private `mktemp -d` directory (all
  children write there through `TMPDIR`), half-written staging directories and
  partial archives.
* Stages run in the background and `run.sh` waits for them, so Ctrl-C is
  handled immediately even during a long MAFFT run.
* Restart without manual cleanup: every stage overwrites its outputs, data are
  promoted from staging directories only when complete, and a stale PID file
  from a killed run is ignored. Two simultaneous runs (e.g. cron + manual) are refused.

### Example of a pipe chain on the data (stage 7)

```bash
tail -n +2 results/mutations.csv | cut -d, -f3-5,7 | grep ',definitive$' \
  | awk -F, '{ print $1 ":" $2 ">" $3 }' | sort | uniq -c | sort -k1,1nr | awk 'NR <= 20'
```

## 4. Files produced

| File | Content |
|---|---|
| `results/summary.csv` | one row per distinct confirmed mutation: `mutation, position, ref, alt, type, n_sequences, frequency` (sorted by frequency) |
| `results/mutation_types.csv` | number of events and sequences per mutation type |
| `results/statistics_summary.json` | denominator, number of confirmed / uncertain calls, number of distinct mutations |
| `results/mutation_analysis.json` | Task 3 run record: MAFFT version, parameters, counts, skip reasons |
| `results/mutations.csv` | every mutation event (≈ 10 M rows, ≈ 800 MB; not archived, its SHA-256 is in `checksums.sha256`) |
| `results/mutation_analysis_skipped.csv` | sequences not analysed and why |
| `results/preprocessing/` | `cleaning_report.json` and audit CSVs (excluded non-HBV, rejected, duplicates, identical sequences…) |
| `results/figures/top_mutations.png` | the 20 most frequent mutations |
| `results/figures/mutations_along_genome.png` | mutation density along the 3,182-nt genome (50-nt bins) |
| `results/figures/mutation_types.png` | events per mutation type |
| `results/shell/*.tsv` | shell profiles and the shell-vs-Python cross-check |
| `results/checksums.sha256` | SHA-256 of all result files |
| `report/dashboard.html` | interactive, self-contained dashboard (open in a browser) |
| `logs/run_*.log` | timestamped log of each run |
| `archives/hbv_results_*.tar.gz` | dated archive of `results/` + dashboard |
| `data/processed/` | cleaned HBV FASTA + metadata + manifest |

## 5. Directory structure

```text
HBV-Genome-Pipeline/
├── run.sh                    single entry point (getopts, trap, logging, archive)
├── download.sh               shortcut to data-collection/download.sh
├── preprocess.py             stage 4 — HBV filter and validation
├── analyze_mutations.py      stage 5 — MAFFT alignment, mutation calling
├── statistics.py             stage 6 — statistics and figures
├── requirements.txt          pinned Python dependencies
├── scripts/                  Bash helpers
│   ├── lib.sh                reusable functions (log, die, require_*, run_stage, checksums, archive…)
│   ├── shell_analysis.sh     grep/sed/awk/cut/sort/uniq pipelines (stages 3 and 7)
│   ├── setup_env.sh          creates .venv and checks system tools
│   ├── install_cron.sh       installs/removes the weekly cron job
│   ├── cron/hbv-pipeline.cron  crontab template
│   └── make_submission.sh    builds Projet_Maladie_Groupe_XX.zip with LIEN_GIT.txt
├── data-collection/          stage 1 (download, validation, normalisation) + its tests
│   ├── download.sh  src/  scripts/  tests/
│   └── data/                 created by run.sh: raw/ (read-only), processed/, checksums/
├── data/processed/           cleaned HBV dataset (stage 4)
├── results/                  all numerical results, figures, shell outputs, checksums
├── report/                   make_dashboard.py, dashboard.html, rapport.md, report.md
├── logs/                     timestamped run logs
├── archives/                 dated tar.gz archives of the results
├── tests/                    unit + integration tests
├── README.md  SOURCES.md  PREPROCESSING.md  TASK3.md  rapport.pdf
```

`data/`, `results/`, `logs/` and `archives/` are (re)created by `run.sh`.

## 6. Scheduled run (cron)

The template `scripts/cron/hbv-pipeline.cron` runs the pipeline **every Sunday
at 03:00** with a fresh NCBI snapshot (`./run.sh -f`). Install it with:

```bash
bash scripts/install_cron.sh -p    # preview the crontab line
bash scripts/install_cron.sh       # install (idempotent)
crontab -l                         # check
bash scripts/install_cron.sh -r    # remove
```

The installed line is
`0 3 * * 0  cd "/path/to/HBV-Genome-Pipeline" && ./run.sh -f >> "/path/to/HBV-Genome-Pipeline/logs/cron.log" 2>&1`.
Each scheduled run writes its own `logs/run_*.log` and a new dated archive;
`run.sh` refuses to start if a previous run is still in progress.

## 7. Permissions

* All scripts are executable (`chmod +x`, mode `100755` in git).
* After verification, `data-collection/data/raw/` (the NCBI package) is made
  **read-only** (`chmod -R a-w`). `run.sh` restores write permission only for
  the download stage and locks the data again right after.

## 8. Tests

```bash
bash tests/run_tests.sh               # all suites (≈ 110 tests, offline, ~30 s)
bash tests/run_integration_test.sh    # only the end-to-end run.sh tests
```

* `data-collection/tests/` — NCBI client (mocked HTTP), validator, loader, shell scripts.
* `tests/preprocessing/` — taxonomy filter, IUPAC handling, duplicates, manifest checks.
* `tests/task3/` — known-answer alignments, ambiguous bases, batch-failure recovery
  (MAFFT test double), plus one real-MAFFT test when MAFFT is installed.
* `tests/task4/` — frequencies, denominator, streaming aggregation, figures.
* `tests/integration/` — runs the real `run.sh` on a tiny fixture with a MAFFT
  stub: all stages, getopts errors, checksum tampering (exit 5), read-only raw
  data, log and archive creation, SIGTERM cleanup (exit 143) and restart.

The tests never contact NCBI and never touch the real dataset.

## 9. Results of the reference run

The reference snapshot (NCBI, taxon 10404, retrieved 14 November 2026) gave:
142,288 records → **140,889 HBV records** after filtering → **128,406
sequences analysed** → **10,567,649 mutation events** (10,098,062 definitive)
→ **16,235 distinct confirmed mutations**. Most frequent: `530:T>C` (61.5 %),
`535:A>C` (57.2 %), `520:C>A` (56.0 %), `576:C>A` (45.6 %), `522:G>A` (45.1 %).
NCBI Virus grows continuously, so a new download gives slightly different
counts; `SOURCES.md` and `results/mutation_analysis.json` record what was used.

## 10. Limitations

* Most records (≈ 90 %) are **partial** genomes; frequencies are computed over
  the sequences that were analysed, whatever region they cover.
* HBV has 8+ genotypes (A–H); the reference NC_003977.2 is one strain (ayw).
  Many "frequent mutations" are therefore **genotype differences**, not
  disease-acquired mutations. NCBI metadata has no usable genotype column, so
  no per-genotype comparison is made.
* Alignment uses heuristic thresholds (≥ 100 nt, ≥ 80 aligned bases, ≥ 70 %
  identity), keeps input orientation (no reverse-complement search) and a
  linear coordinate system although the genome is circular.
* Identical sequences submitted under different accessions are all counted.
* Positions are NC_003977.2 coordinates; they are not translated into the
  classical HBV numbering, genes or amino acids, and no clinical
  interpretation is made.

Details: [PREPROCESSING.md](PREPROCESSING.md), [TASK3.md](TASK3.md),
[data-collection/README.md](data-collection/README.md).

## 11. Work distribution

> **TODO (group):** replace with the real names and the group number before submission.

| Member | Main responsibility |
|---|---|
| Student 1 — *Name* | Stage 1: NCBI data collection, validation, checksums (`data-collection/`) |
| Student 2 — *Name* | Stage 2: preprocessing and HBV filtering (`preprocess.py`) |
| Student 3 — *Name* | Stage 3: MAFFT alignment and mutation detection (`analyze_mutations.py`) |
| Student 4 — *Name* | Stage 4: statistics, figures, dashboard (`statistics.py`, `report/`) |
| Student 5 — *Name* | Orchestration and system administration (`run.sh`, `scripts/`, cron, tests, report) |

## 12. Submission

```bash
bash scripts/make_submission.sh -g XX     # -> dist/Projet_Maladie_Groupe_XX.zip
```

The archive contains `LIEN_GIT.txt` (repository URL + `HEAD` commit hash), the
code, `README.md`, `SOURCES.md` and `rapport.pdf`, but **no data**: `run.sh`
downloads and verifies the dataset itself. Commit and push before building it.
