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

**Phase 1 (data collection) is implemented; the analysis stages are still empty.**

- ✅ `data-collection/` — NCBI retrieval, local-file pathway, validation,
  checksums, documentation. `download.sh` at the root delegates to it.
- ⬜ Empty (0 bytes): `preprocess.py`, `analyze_mutations.py`, `statistics.py`,
  `run.sh`, `report/report.md`.
- `results/mutations.csv`, `results/summary.csv` remain empty placeholders.
- `data-collection/data/raw/` and `processed/` hold a real 405 MB dataset and are
  git-ignored; `data-collection/data/checksums/` is committed.
- No tests directory, no CI, no lint config.

**Consequence:** there is no established test suite to imitate. When you implement
an analysis stage, follow the conventions in `data-collection/src/` (stdlib-only,
`argparse`, `main()`, explicit non-zero exits, atomic output promotion).

## 3. Repository layout

| Path | Purpose | Committed? |
|---|---|---|
| `download.sh` | Root entry point; delegates to `data-collection/download.sh` | yes |
| `data-collection/download.sh` | Stage orchestrator: NCBI default, `--file` alternative | yes |
| `data-collection/src/ncbi_client.py` | NCBI Datasets REST client (retry, atomic download) | yes |
| `data-collection/src/file_loader.py` | Local + NCBI normalisation to `sequences.fasta` / `metadata.csv` | yes |
| `data-collection/src/dataset_validator.py` | ZIP/FASTA/size checks + SHA-256 | yes |
| `data-collection/scripts/validate_dataset.sh` | Standalone validation entry point | yes |
| `data-collection/dataset_documentation.md` | Dataset provenance, scope, checksums | yes |
| `data-collection/data/raw/`, `processed/` | Downloaded + normalised data | **no** (git-ignored) |
| `data-collection/data/checksums/` | SHA-256 lists and validation reports | yes |
| `preprocess.py` | Clean/filter/align sequences (Phase 2) | yes |
| `analyze_mutations.py` | Detect mutations/variants vs. reference (Phase 3) | yes |
| `statistics.py` | Aggregate stats + figures into `results/figures/` (Phase 4) | yes |
| `run.sh` | Orchestrator for the analysis stages | yes |
| `data/` | Legacy root data dir (only `.gitkeep`) | **no** |
| `results/` | Machine-readable outputs (`mutations.csv`, `summary.csv`) | headers only |
| `results/figures/` | Generated plots | no (only `.gitkeep`) |
| `report/report.md` | Human-readable final report | yes |
| `SOURCES.md` | Provenance: URLs, accessions, retrieval dates | yes |

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
- `datasets` CLI v18.39.0 is installed at `~/.local/bin/datasets` (optional; the
  Phase 1 code uses only the standard library and never calls it).
- If a stage needs one of these, add it to `requirements.txt` and/or document the system
  dependency in `README.md` with install instructions — do not silently assume it exists.
- No virtualenv is present; create `.venv/` locally if needed (must be git-ignored).

## 9. How to verify a change

There is no test suite yet. Until one exists:

```bash
bash -n download.sh run.sh data-collection/download.sh data-collection/scripts/validate_dataset.sh
python3 -m py_compile *.py data-collection/src/*.py
bash run.sh                         # analysis pipeline (empty until Phase 2)

# Phase 1 smoke test (small taxon, no size gate)
./download.sh --taxon 10406 --min-mb 0
./download.sh --file data-collection/tests/data/sample.fasta
```

Then inspect `results/mutations.csv`, `results/summary.csv`, and `results/figures/`.
When you add tests, put them in `tests/` and prefer `pytest`.

## 10. Open decisions (not yet settled)

- **Alignment strategy** for `analyze_mutations.py` — pure Python position-by-position
  comparison (no dependencies, works today) vs. installing MAFFT/BLAST (stronger
  multiple alignment, needs network + `apt`). Do not assume one; confirm before building
  that stage.
- **Taxonomic scope for downstream stages** — the Phase 1 dataset is the
  *Hepadnaviridae* family (taxid 10404) and contains ~1 % non-HBV records. Decide
  whether analysis filters to HBV (`lineage` contains `Hepatitis B virus`) or keeps
  the family.
- **CSV schemas** in section 6 are proposals, not yet implemented or validated.
  Phase 1 outputs (`sequences.fasta`, `metadata.csv`) are defined in
  `data-collection/dataset_documentation.md` §7.

## 11. Definition of done for a stage

- Script is non-empty, runnable, and documented (usage + inputs + outputs).
- Reads only from `data/` (or documented inputs) and writes only to `results/`,
  `results/figures/`, or `report/`.
- Failure modes exit non-zero with an actionable message.
- `SOURCES.md` / `README.md` updated if data sources or dependencies changed.
- No generated data or large files staged for commit.
