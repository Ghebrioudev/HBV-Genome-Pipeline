# AGENTS.md — HBV-Genome-Pipeline

Instructions for AI agents working in this repository. Read this first.

## 1. What this project is

A reproducible **bioinformatics pipeline** for analysing **Hepatitis B Virus (HBV)**
genomic data, written in **Shell + Python 3**.

The pipeline is meant to go from raw public HBV sequences to a written analysis report:

```
download  →  preprocess  →  mutation analysis  →  statistics/figures  →  report
download.sh   preprocess.py   analyze_mutations.py    statistics.py        report/report.md
                             run.sh orchestrates all of the above
```

## 2. Current state (important)

This repository is an **initialized scaffold, not a working pipeline.**

- All entry-point scripts exist but are **empty (0 bytes)**: `download.sh`,
  `preprocess.py`, `analyze_mutations.py`, `statistics.py`, `run.sh`.
- `SOURCES.md` and `report/report.md` are also empty.
- Output files exist but are **empty placeholders**: `results/mutations.csv`,
  `results/summary.csv`, `results/figures/` (only `.gitkeep`).
- `data/` contains only `.gitkeep` — no sequence data is committed.
- No tests, no `requirements.txt`, no CI, no lint/format config.
- Only two commits exist (`chore: initialize HBV project structure`).

**Consequence:** there is no established code style, no API surface, and no test suite
to imitate. When you implement a script, you are defining the convention. Keep it simple,
documented, and consistent with the layout below.

## 3. Repository layout

| Path | Purpose | Committed? |
|---|---|---|
| `download.sh` | Fetch reference + sample HBV sequences from public sources into `data/` | yes |
| `preprocess.py` | Clean/filter/align sequences, produce analysis-ready inputs | yes |
| `analyze_mutations.py` | Detect mutations/variants relative to reference | yes |
| `statistics.py` | Aggregate stats + generate figures into `results/figures/` | yes |
| `run.sh` | Orchestrator: runs the stages in order, stops on error | yes |
| `data/` | Raw + intermediate sequence data | **no** (only `.gitkeep`) |
| `results/` | Machine-readable outputs (`mutations.csv`, `summary.csv`) | headers only |
| `results/figures/` | Generated plots | no (only `.gitkeep`) |
| `report/report.md` | Human-readable final report | yes |
| `SOURCES.md` | Provenance: URLs, access dates, database versions | yes |

## 4. Ground rules

1. **Never commit data.** `data/` is git-ignored in spirit and must stay out of git.
   Do not `git add data/*` even if `.gitignore` is missing — add a `.gitignore` entry instead
   (`data/*`, `!data/.gitkeep`, plus `__pycache__/`, `*.pyc`, `.venv/`).
2. **`results/` is generated, never hand-edited.** Regenerate by running the pipeline.
   Keep the files tracked (headers only) so the expected schema is visible.
3. **Provenance matters.** Any downloaded sequence must be recorded in `SOURCES.md`
   (source URL, accession IDs, retrieval date, database version).
4. **Reproducibility.** No hard-coded absolute paths; resolve paths relative to the repo
   root. Pin tool/package versions when you add dependencies.
5. **Do not delete empty placeholder files** (`.gitkeep`, empty CSVs) unless you are
   replacing them with real generated content.
6. **Scripts must fail loudly.** Use `set -euo pipefail` in every shell script;
   raise/exit non-zero on missing inputs in Python.
7. **Science first.** Never fabricate sequences, mutation counts, or statistics.
   If data is unavailable, report the failure — do not invent results for the report.

## 5. Coding conventions

- **Shell:** `#!/usr/bin/env bash`, `set -euo pipefail`, quote variables, print a clear
  usage/help block, and echo stage progress to stderr.
- **Python:** Python 3.12 target; standard library first, `pandas`/`matplotlib` only when
  needed. Use `argparse` for CLI args, `pathlib.Path` for paths, `main()` + `if __name__ ==
  "__main__"` entry point, and module-level docstrings describing inputs/outputs.
- Keep each stage independently runnable (`python preprocess.py --input ...`) *and*
  callable from `run.sh`.
- Prefer small, single-responsibility functions with type hints on public ones.

## 6. Expected output schemas (define these if you touch them)

`results/mutations.csv` — one row per detected variant:

```
sample_id,position,reference_base,alternate_base,amino_acid_change,region,gene,effect
```

`results/summary.csv` — one row per sample:

```
sample_id,genotype,length,num_mutations,mutations_per_kb,dominant_region
```

If you need to change a schema, update this file **and** `report/report.md` in the same change.

## 7. Branching model

`main` is the integration branch. Feature branches already exist on `origin`:

- `feature/data-collection` → `download.sh`, `SOURCES.md`
- `feature/data-preprocessing` → `preprocess.py`
- `feature/sequence-comparison` → `analyze_mutations.py`
- `feature/statistics-visualization` → `statistics.py`, `results/figures/`
- `feature/pipeline-integration` → `run.sh`, `report/report.md`

Work on the branch matching your stage. Do not push to `main` directly, and do not
rewrite or force-push shared branches.

## 8. Environment

- Python 3.12.3, GNU bash 5.2.21 available.
- **Not installed:** `blastn`, `makeblastdb`, `samtools`, `mafft`, `minimap2`,
  `pandas`, `matplotlib`, `biopython`.
- If a stage needs one of these, add it to `requirements.txt` and/or document the system
  dependency in `README.md` with install instructions — do not silently assume it exists.
- No virtualenv is present; create `.venv/` locally if needed (must be git-ignored).

## 9. How to verify a change

There is no test suite yet. Until one exists:

```bash
bash -n download.sh run.sh          # shell syntax check
python3 -m py_compile *.py          # python syntax check
bash run.sh                         # end-to-end smoke run
```

Then inspect `results/mutations.csv`, `results/summary.csv`, and `results/figures/`.
When you add tests, put them in `tests/` and prefer `pytest`.

## 10. Open decisions (not yet settled)

- **Alignment strategy** for `analyze_mutations.py` — pure Python position-by-position
  comparison (no dependencies, works today) vs. installing MAFFT/BLAST (stronger
  multiple alignment, needs network + `apt`). Do not assume one; confirm before building
  that stage.
- **Data source / accessions** for `download.sh` — no specific NCBI accessions or
  reference genome have been chosen yet.
- **CSV schemas** in section 6 are proposals, not yet implemented or validated.

## 11. Definition of done for a stage

- Script is non-empty, runnable, and documented (usage + inputs + outputs).
- Reads only from `data/` (or documented inputs) and writes only to `results/`,
  `results/figures/`, or `report/`.
- Failure modes exit non-zero with an actionable message.
- `SOURCES.md` / `README.md` updated if data sources or dependencies changed.
- No generated data or large files staged for commit.
