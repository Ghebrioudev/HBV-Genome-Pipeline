# Plan compliance audit — HBV-Genome-Pipeline

Audit date: 2026-11-08 · Branch: `feature/pipeline-integration-final`
Purpose: check the delivered code and artifacts against the five-task project
plan (Rayan, Khiro, Achraf, A9ba, New Member) and list what is still missing.

Status legend: **MET** = requirement satisfied and evidenced · **PARTIAL** =
implemented but not evidenced, not reproducible, or conditionally skipped ·
**MISSING** = required deliverable absent.

---

## 1. Source plan (verbatim scope)

| # | Owner | Task | Required bullets | Deliverables |
|---|---|---|---|---|
| 1 | Rayan | Data Collection | public HBV dataset on NCBI; size within 100 MB–1 GB; download + integrity script; document source/metadata | `download.sh`, dataset documentation |
| 2 | Khiro | Data Preprocessing | inspect format; check invalid/incomplete/duplicate; filter/prepare; record counts before/after | `preprocess.py`, cleaning report |
| 3 | Achraf | Sequence Comparison | align to a reference; identify substitutions and other variations; record positions + sequence IDs; export results | `analyze_mutations.py`, mutation table |
| 4 | A9ba | Statistics & Visualization | mutation frequencies; genotype comparison *if reliable metadata*; graphs + summary tables; interpret main patterns without clinical claims | `statistics.py`, CSV summaries, graphs |
| 5 | New Member | Integration & Documentation | `run.sh` in correct order; scripts work together + handle errors; methods/reproducibility/limitations sections; assemble structure and test from a clean start | `run.sh`, integration tests, report sections |

## 2. Audit method

Each requirement is checked against an artifact, not a claim. The whole audit is
offline and reproducible from the repository root:

```bash
# 1. deliverables exist and are non-empty
ls -l download.sh preprocess.py analyze_mutations.py statistics.py run.sh \
      results/mutations.csv results/summary.csv
ls -l results/figures/ results/preprocessing/

# 2. tests (contract-level behaviour)
PATH="$PWD/.venv/bin:$PATH" bash tests/run_tests.sh

# 3. size window (Task 1)
du -sm data-collection/data/raw data-collection/data/processed
python3 data-collection/src/dataset_validator.py --help   # MIN_MB=100, MAX_MB=1024

# 4. integrity
(cd data-collection && sha256sum -c data/checksums/dataset.sha256)

# 5. cleaning before/after
python3 -c "import json;c=json.load(open('results/preprocessing/cleaning_report.json'))['counts'];print(c['input_fasta_records'],c['retained_records'])"

# 6. orchestration order and failure behaviour
./run.sh --help; ./run.sh --bogus; echo $?   # 0 then 2

# 7. end-to-end (requires MAFFT)
command -v mafft && ./run.sh --skip-download
```

Pass criteria: non-empty deliverable, test suite green, size inside the window,
checksums verify, before/after counts recorded, unknown option exits non-zero,
clean run produces mutations + summaries + figures without manual steps.

## 3. Result matrix

### Task 1 — Data Collection · **MET (one provenance defect)**

| Requirement | Evidence | Status |
|---|---|---|
| Public HBV dataset from NCBI | `SOURCES.md` (Hepadnaviridae taxid 10404 / HBV 10407), `data/raw/retrieval.json` | MET |
| 100 MB–1 GB | raw 436 MB, processed 191 MB; `dataset_validator.py` gates 100–1024 MB (hard) | MET |
| Download + integrity script | `download.sh`, SHA-256 lists, `scripts/validate_dataset.sh` | MET |
| Source + metadata documented | `SOURCES.md`, `dataset_documentation.md` §1–§12 | MET |
| Committed checksums match the data on disk | `sha256sum -c` → **4 FAILED** | **MISSING** |

The committed checksums describe a different snapshot (taxid 10407, 140,855
records) than the files on this machine (taxid 10404, 142,254 records).

### Task 2 — Data Preprocessing · **MET (with scientific caveats)**

| Requirement | Evidence | Status |
|---|---|---|
| Inspect downloaded format | `dataset_documentation.md` §6–§7 | MET |
| Check invalid / incomplete / duplicate | `results/preprocessing/{rejected_records,review_records,duplicate_accessions,identical_sequences}.csv` | MET |
| Filter and prepare | taxonomy filter + IUPAC validation in `preprocess.py` | MET |
| Record counts before/after | `cleaning_report.json`: `input_fasta_records` 140,855 → `retained_records` 140,855 | PARTIAL |
| Cleaning report present | `results/preprocessing/cleaning_report.json` | MET |

Caveats:
- The committed run removed **0** records. Because the snapshot was already
  HBV-only, the filter/dedup logic is exercised only by unit tests, so the
  "before vs after" evidence is degenerate.
- **31,288 extra identical sequences** in 11,912 groups are retained (different
  accessions, same sequence). They remain in the Task 3/4 denominators and can
  bias mutation frequencies.
- 127,395 of 140,855 records are `PARTIAL`. Retaining them is a documented
  choice; eligibility is deferred to Task 3 thresholds.

### Task 3 — Sequence Comparison · **MISSING deliverable**

| Requirement | Evidence | Status |
|---|---|---|
| Align to a reference | `NC_003977.2` via MAFFT `--6merpair --addfragments` (`TASK3.md`) | MET (code) |
| Identify substitutions + other variations | `substitution`, `deletion`, `insertion`, `ambiguous` in `analyze_mutations.py` | MET (code) |
| Record positions + sequence IDs | schema `sample_id,position,reference_base,observed_base,...` | MET (code) |
| Export comparison results | `results/mutations.csv` = **0 bytes** | **MISSING** |
| Reproducible on this machine | `mafft` not installed | **MISSING** |

Also:
- `TASK3.md` states "the dataset was not aligned here and no biological mutation
  counts are claimed", yet the committed `results/mutation_analysis.json` claims
  10,563,024 rows analysed from `/home/ubuntu/SysProject/HBV-Preprocessing-Push`.
  The committed JSON and absent CSV are mutually inconsistent (the script
  publishes both together).
- `results/mutation_analysis_skipped.csv` (611 KB) is committed without its CSV
  counterpart.

### Task 4 — Statistics & Visualization · **PARTIAL**

| Requirement | Evidence | Status |
|---|---|---|
| Mutation frequencies | `compute_frequencies()` + denominator logic in `statistics.py` | MET (code) |
| Genotype comparison *if reliable metadata* | code supports it; `metadata.csv` has **no genotype column** → skipped with a warning | PARTIAL |
| Graphs + summary tables | `plot_*` functions exist; `results/figures/` holds only `.gitkeep`, `results/summary.csv` is 0 bytes | **MISSING** |
| Interpret main patterns, no clinical claims | no interpretation anywhere; `report/report.md` explicitly defers it | **MISSING** |

Schema defect: `statistics.py` writes `summary.csv` as **one row per mutation**
(`mutation,position,ref,alt,type,n_sequences,frequency`), but `AGENTS.md` §6
still defines it as **one row per sample**. The schema changed without updating
the contract.

### Task 5 — Integration & Documentation · **PARTIAL**

| Requirement | Evidence | Status |
|---|---|---|
| `run.sh` in correct order | download → preprocess → mutations → statistics | MET |
| Scripts work together + handle errors | `set -euo pipefail`, `require_file`, unknown option exits 2; 99 tests green | MET |
| Methods / reproducibility / limitations | `report/report.md` has all three | PARTIAL |
| Final structure assembled | 31 source/test/doc files | MET |
| Test from a clean start | `tests/run_integration_test.sh` runs `PYTHON=true run.sh`, so **no stage actually executes** | PARTIAL |

Fixes already applied during this audit: `.gitignore` restored (the merge had
reduced it to `__pycache__/`, exposing the 630 MB dataset), executable bits
restored on `run.sh`, `tests/run_integration_test.sh`, `preprocess.py`,
`analyze_mutations.py`, and `requirements.txt` added.

## 4. Missing things (prioritised)

### P0 — blocks a defensible delivery
1. **Mutation table** (`results/mutations.csv`) — required Task 3 deliverable,
   currently 0 bytes.
2. **MAFFT installed** — Phase 3 cannot run or be reproduced without it.
3. **One canonical dataset snapshot** — code/docs/checksums must all describe
   the same fetch; today 142,254 (disk, taxid 10404) vs 140,855 (docs, taxid
   10407).
4. **One real end-to-end run** from a clean start, producing mutations, summary
   CSVs and figures.

### P1 — required by the plan, absent or incomplete
5. **Task 4 outputs** — `summary.csv`, `mutation_types.csv`,
   `statistics_summary.json`, `results/figures/*.png`.
6. **Results + interpretation section** in `report/report.md` (patterns only, no
   clinical claims).
7. **Genotype statement** — either derive genotype from `organism_name`/`lineage`
   or state in the report that the metadata cannot support it.
8. **Reconcile `summary.csv` schema** between `statistics.py` and `AGENTS.md` §6.

### P2 — quality and provenance
9. **Strengthen the integration test** to execute the real stdlib stages on a
   small fixture (`PYTHON=true` proves nothing about the pipeline).
10. **Resolve the Task 3 provenance contradiction** — remove the foreign-machine
    JSON, or regenerate it here and delete the stale "no counts claimed" note.
11. **Document the duplicate-sequence policy** (31,288 redundant records) and,
    ideally, report a sensitivity check with and without them.
12. **Stop committing large generated outputs** (`identical_sequences.csv`,
    6.7 MB / 43k rows) per the "results are generated" rule.
13. Optional, outside the plan: CI and a lint configuration.

## 5. Remediation order

```bash
# 1. choose one snapshot, then re-validate
./download.sh --taxon 10407        # overwrites data-collection/data/processed/
(cd data-collection && sha256sum -c data/checksums/dataset.sha256)

# 2. install the Phase 3 dependency
sudo apt install mafft && mafft --version

# 3. full run
./run.sh --skip-download

# 4. confirm deliverables
wc -l results/mutations.csv results/summary.csv
ls -l results/figures/
```

Then write §Results in `report/report.md` from `results/summary.csv` and the
figures, and add the genotype/duplicate-sequence limitations.

## 6. Re-audit checklist

- [ ] `sha256sum -c` passes for the canonical snapshot
- [ ] `mafft --version` prints ≥ 7.487
- [ ] `./run.sh --skip-download` exits 0 from a clean `results/`
- [ ] `results/mutations.csv` has a header + rows
- [ ] `results/summary.csv` matches the documented schema
- [ ] `results/figures/` contains the three PNGs
- [ ] `report/report.md` has Methods, Results, Reproducibility, Limitations
- [ ] `bash tests/run_tests.sh` stays green (99 tests)
