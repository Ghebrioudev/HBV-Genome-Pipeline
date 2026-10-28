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
- Phase 1 has a dependency-free `unittest` suite: `data-collection/tests/`
  (57 tests, offline, ~0.2 s). There is still no CI and no lint config.

**Consequence:** Phase 1 established the conventions. When you implement an
analysis stage, follow `data-collection/src/` (stdlib-only, `argparse`, `main()`,
explicit non-zero exits, atomic output promotion) and add `test_*.py` modules
alongside the existing offline suite.

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
| `data-collection/tests/data/` | Real 2-record fixtures for the local pathway | yes |
| `data-collection/tests/test_*.py` | Phase 1 unittest suite (offline, stdlib only) | yes |
| `data-collection/tests/run_tests.sh` | Test runner | yes |
| `data/` | Legacy root data dir (only `.gitkeep`) | **no** |
| `results/` | Machine-readable outputs (`mutations.csv`, `summary.csv`) | headers only |
| `results/figures/` | Generated plots | no (only `.gitkeep`) |
| `report/report.md` | Human-readable final report | yes |
| `SOURCES.md` | Provenance: URLs, accessions, retrieval dates | yes |

### Phase 1 output contract (consume this in later stages)

Everything downstream must read from **`data-collection/data/processed/`** and
must not re-implement retrieval:

| File | Format | Notes |
|---|---|---|
| `sequences.fasta` | FASTA | 142,254 records, 131,955,960 nt, 140 MB |
| `metadata.csv` | CSV | 142,254 rows, same order as the FASTA |
| `manifest.json` | JSON | record count + SHA-256 of the outputs |

`metadata.csv` columns, in order: `accession, length, completeness, is_annotated,
gene_count, protein_count, source_database, release_date, update_date,
organism_name, tax_id, isolate, lineage, submitter`.

Key facts agents must not re-derive by hand:

- The dataset is **family-wide** (*Hepadnaviridae*, taxid 10404): 140,866 records
  (99.0 %) are HBV, 1,388 (1.0 %) are other hepadnaviruses. Filter with
  `'Hepatitis B virus' in row['lineage']` if a stage needs HBV only.
- 128,154 records are `PARTIAL`, 14,100 are `COMPLETE`.
- The HBV reference genome for mutation calling is `NC_003977.2`.
- Full provenance and checksums: `data-collection/dataset_documentation.md`.

## 4. Ground rules

1. **Never commit data.** Raw and normalised data must stay out of git:
   `data-collection/data/raw/`, `data-collection/data/processed/`, and `data/*`.
   `.gitignore` already covers these; verify with `git status --short --ignored`
   before committing. Only `.gitkeep` and `data-collection/data/checksums/` are tracked.
   Do not force-add ignored data with `git add -f`.
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
8. **`download.sh` replaces the processed dataset.** Both the NCBI and `--file`
   pathways promote their result into `data-collection/data/processed/`, overwriting
   the previous dataset (atomically, only on success). Never run
   `./download.sh --file ...` on a whim: it will swap the 142k-record dataset for
   the supplied file. Use `data-collection/tests/data/` fixtures only when you
   intend to do that, and re-run `./download.sh` afterwards to restore the full set.
9. **Do not add network calls to tests or analysis stages.** Only
   `src/ncbi_client.py` talks to the network. Everything else must work offline on
   the local files.

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

### Remote write access (currently blocked)

- The git remote is `https://github.com/Ghebrioudev/HBV-Genome-Pipeline.git`, owned
  by **Ghebrioudev**. The configured account (`berkal749`) has `push: false`, so
  every push — including creating new branches — returns HTTP 403.
- Work locally and commit; pushing is deferred until access changes (switch `gh`
  account, add the account as a collaborator, or push to a fork).
- `origin/feature/data-collection` is an **unrelated history** containing only
  `README.md` and is 1 commit "behind" while the local branch is ahead. Once write
  access exists, syncing requires:
  `git push --force-with-lease origin feature/data-collection`.
- `origin/main` is 2 commits behind the local scaffold and never contains
  `AGENTS.md` (by design).

## 8. Environment

- Python 3.12.3, GNU bash 5.2.21 available.
- **Not installed:** `blastn`, `makeblastdb`, `samtools`, `mafft`, `minimap2`,
  `pandas`, `matplotlib`, `biopython`.
- `datasets` CLI v18.39.0 is installed at `~/.local/bin/datasets` (optional; the
  Phase 1 code uses only the standard library and never calls it).
- If a stage needs one of these, add it to `requirements.txt` and/or document the system
  dependency in `README.md` with install instructions — do not silently assume it exists.
- No virtualenv is present; create `.venv/` locally if needed (must be git-ignored).
- Data is already on disk: `data-collection/data/` is ~630 MB total (436 MB raw,
  191 MB processed). No network access is needed to work on analysis stages.

## 9. How to verify a change

There is no test suite yet. Until one exists:

```bash
bash -n download.sh run.sh data-collection/download.sh data-collection/scripts/validate_dataset.sh
python3 -m py_compile *.py data-collection/src/*.py

# Phase 1 test suite (stdlib unittest, offline)
data-collection/tests/run_tests.sh

# dataset is intact
(cd data-collection && sha256sum -c data/checksums/dataset.sha256)
```

The suite covers the validator, the local and NCBI loaders, the HTTP client
(mocked — **no network**), and the shell CLI contract. It must stay green; add a
test with every behaviour change.

> `./download.sh --file <x>` **overwrites** `data-collection/data/processed/`.
> Use `dataset_validator.py` (or the fixture tests) for read-only checks; only run
> `download.sh` when you actually mean to replace the dataset.

Then inspect `results/mutations.csv`, `results/summary.csv`, and `results/figures/`.
Add new tests as `test_*.py` modules in `data-collection/tests/` (Phase 1) or a
root `tests/` directory (later stages). Use `unittest`: `pytest` is **not
installed** and must not be assumed.

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
- Reads only from `data-collection/data/processed/` (or documented inputs) and
  writes only to `results/`, `results/figures/`, or `report/`.
- Failure modes exit non-zero with an actionable message.
- Outputs are written atomically or staged, so a failed run does not corrupt a
  previous good result.
- `SOURCES.md` / `README.md` / `dataset_documentation.md` updated if data sources,
  schemas or dependencies changed.
- No generated data or large files staged for commit; `git status --short --ignored`
  shows data as ignored.
