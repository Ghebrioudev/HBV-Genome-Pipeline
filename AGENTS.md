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

**Phases 1–4 (data collection, preprocessing, mutation detection, and statistics) are implemented; Phase 5 remains unimplemented.**

- ✅ `data-collection/` — NCBI retrieval, local-file pathway, validation,
  checksums, documentation. `download.sh` at the root delegates to it.
- ✅ `preprocess.py` — Phase 2 preprocessing (stdlib-only); reads Task 1 processed inputs and writes cleaned HBV data plus audit reports.
- ✅ `analyze_mutations.py` — Task 3 reference-anchored alignment and mutation calling.
- ✅ `statistics.py` — Task 4 mutation statistics and figures; `run.sh` orchestrates the pipeline.
- `tests/preprocessing/` — Phase 2 offline `unittest` coverage; `tests/task3/` — Phase 3 offline tests; `tests/task4/` — Phase 4 offline tests; `tests/integration/` — clean-start fixture that runs all four stages through `run.sh` with a deterministic MAFFT stub; `tests/run_tests.sh` runs all five suites (102 tests).
- `results/mutations.csv` is generated locally by Task 3 (~813 MB, 10.5 M rows)
  and is git-ignored; regenerate it with `python3 analyze_mutations.py`. Task 4
  writes summary statistics and figures to `results/`.
- The bundled Task 1 snapshot occupies approximately 431 MB under
  `data-collection/data/raw/` and 188 MB under `data-collection/data/processed/`;
  both input directories are git-ignored. Sizes vary with the retrieved snapshot.
  `data-collection/data/checksums/` is committed.
- Phase 1 has a dependency-free `unittest` suite: `data-collection/tests/` (57 tests).
  Phase 2 has 20 offline tests in `tests/preprocessing/`; Task 3 and Task 4 have
  focused offline suites. There is still no CI or lint config.

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
| `preprocess.py` | Taxonomy filtering, sequence validation, duplicate detection and reporting (Phase 2; no alignment) | yes |
| `PREPROCESSING.md` | Phase 2 scientific rules, outputs and execution | yes |
| `tests/preprocessing/` | Phase 2 offline unit tests | yes |
| `tests/run_tests.sh` | Runs the Phase 1–4 test suites | yes |
| `analyze_mutations.py` | MAFFT reference alignment and mutation calling (Phase 3) | yes |
| `TASK3.md` | Phase 3 method, inputs, outputs, execution and limitations | yes |
| `tests/task3/` | Phase 3 offline unit tests and optional MAFFT integration check | yes |
| `statistics.py` | Aggregate stats + figures into `results/figures/` (Phase 4) | yes |
| `run.sh` | Orchestrator for the analysis stages | yes |
| `data-collection/tests/data/` | Real 2-record fixtures for the local pathway | yes |
| `data-collection/tests/test_*.py` | Phase 1 unittest suite (offline, stdlib only) | yes |
| `data-collection/tests/run_tests.sh` | Test runner | yes |
| `data/processed/` | Generated Task 2 HBV FASTA, metadata and manifest (git-ignored) | **no** |
| `results/` | Task 2 audit reports; Task 3 mutation table; Task 4 summaries and figures | mixed |
| `results/figures/` | Generated plots | no (only `.gitkeep`) |
| `report/report.md` | Human-readable final report draft; Phase 5 is not started | yes |
| `SOURCES.md` | Provenance: URLs, accessions, retrieval dates | yes |

**Archive note:** the supplied Task 2 ZIP is a source/data snapshot and does not
include a `.git/` directory. Branch, remote, and commit state must therefore be
verified in the actual Git clone; this archive cannot establish or change it.

### Phase 1 output contract (consume this in later stages)

Everything downstream must read from **`data-collection/data/processed/`** and
must not re-implement retrieval:

| File | Format | Notes |
|---|---|---|
| `sequences.fasta` | FASTA | Current bundled snapshot: 140,855 records, 129,365,332 nt |
| `metadata.csv` | CSV | 140,855 rows; accession order is not guaranteed to match FASTA order |
| `manifest.json` | JSON | record count + SHA-256 of the outputs |

`metadata.csv` columns, in order: `accession, length, completeness, is_annotated,
gene_count, protein_count, source_database, release_date, update_date,
organism_name, tax_id, isolate, lineage, submitter`.

Key facts agents must not re-derive by hand:

- The earlier documented retrieval was family-wide (*Hepadnaviridae*, taxid 10404).
  The **current files bundled in the Task 2 archive are different**: `retrieval.json`
  specifies taxon 10407 and the Task 1 manifest reports 140,855 records. For this
  current snapshot, 127,395 records are `PARTIAL` and 13,460 are `COMPLETE`.
  The current machine-readable retrieval and processed manifest are authoritative.
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

`results/mutations.csv` — one row per detected variant/event:

```
sample_id,reference_accession,position,reference_base,observed_base,variation_type,call_status,aligned_bases,percent_identity,completeness
```

Coordinates are 1-based reference positions; insertions use a between-base anchor such as `125^126`. See `TASK3.md` for the precise event and uncertainty rules.

`results/summary.csv` — one row per distinct mutation, written by Task 4:

```
mutation,position,ref,alt,type,n_sequences,frequency
```

`position` keeps the insertion anchor form (`125^126`). `n_sequences` is the
number of distinct sequences carrying the mutation and `frequency` is
`n_sequences / denominator` (the denominator is metadata sequences minus
skipped sequences unless `--total-sequences` is given).

`results/mutation_types.csv` — one row per variant type
(`type,n_mutations,n_sequences`); `results/statistics_summary.json` — denominator,
call counts and the genotype-comparison flag. Figures are written to
`results/figures/` and are not committed.

If a downstream schema changes, update this file. Reflect the corresponding interpretation in `report/report.md` when Task 5 is implemented; Task 3 does not write the final report.

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

```bash
bash -n download.sh run.sh data-collection/download.sh data-collection/scripts/validate_dataset.sh
python3 -m py_compile *.py data-collection/src/*.py

# Phase 1 test suite (stdlib unittest, offline)
bash tests/run_tests.sh

# dataset is intact
(cd data-collection && sha256sum -c data/checksums/dataset.sha256)
```

The Task 1 suite covers the validator, local and NCBI loaders, the HTTP client
(mocked — **no network**), and the shell CLI contract. Task 2 has its own offline
unit tests in `tests/preprocessing/`. Both suites must stay green; add tests for
every behaviour change.

> `./download.sh --file <x>` **overwrites** `data-collection/data/processed/`.
> Use `dataset_validator.py` (or the fixture tests) for read-only checks; only run
> `download.sh` when you actually mean to replace the dataset.

Then inspect `results/mutations.csv`, `results/summary.csv`, and `results/figures/`.
Add new tests as `test_*.py` modules in `data-collection/tests/` (Phase 1) or a
root `tests/` directory (later stages). Use `unittest`: `pytest` is **not
installed** and must not be assumed.

## 10. Open decisions (not yet settled)

- **Task 3 scope implemented:** whole-reference comparison to `NC_003977.2` with MAFFT `--addfragments`; partial records are eligible by configurable thresholds, terminal gaps are treated as uncovered, and IUPAC ambiguity calls are marked uncertain. See `TASK3.md`.

The Task 2 taxonomic scope is resolved: retain taxid 10407 or records whose NCBI
lineage contains the exact `Hepatitis B virus` taxon; report non-HBV family members
and unresolved taxonomy separately.
- **CSV schemas** in section 6 are proposals, not yet implemented or validated.
  Phase 1 outputs (`sequences.fasta`, `metadata.csv`) are defined in
  `data-collection/dataset_documentation.md` §7.

## 11. Definition of done for a stage

- Script is non-empty, runnable, and documented (usage + inputs + outputs).
- Task 2 reads only from `data-collection/data/processed/` and writes its cleaned
  dataset to `data/processed/` plus audit reports to `results/preprocessing/`.
- Failure modes exit non-zero with an actionable message.
- Outputs are written atomically or staged, so a failed run does not corrupt a
  previous good result.
- `SOURCES.md` / `README.md` / `dataset_documentation.md` updated if data sources,
  schemas or dependencies changed.
- No generated data or large files staged for commit; `git status --short --ignored`
  shows data as ignored.
