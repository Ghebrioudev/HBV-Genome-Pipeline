# Task 3 — Sequence Comparison: completion plan

Branch: `anti-sprint` (from `feature/pipeline-integration-final` @ `7b00d07`)
Owner: Achraf · Deliverables: `analyze_mutations.py` + mutation table

## 1. Current state

The **code is complete and tested** (18 Task 3 tests, 1 real-MAFFT test skipped).
What is missing is the **execution evidence and the required output file**:

| Item | State |
|---|---|
| `analyze_mutations.py` (843 lines) | present, tested |
| Reference `NC_003977.2` in input FASTA | present (line 333150 of `sequences.fasta`) |
| `results/mutations.csv` | **0 bytes** |
| `results/mutation_analysis.json` | present, but from `/home/ubuntu/...`, not this repo |
| `results/mutation_analysis_skipped.csv` | present (611 KB), same foreign run |
| `data/processed/` (Task 2 output) | **absent** — preprocessing never ran here |
| `mafft` executable | **not installed** |

Two gates block a real run:

- **D1 — dataset snapshot.** Disk holds 142,254 records (taxid 10404, family-wide);
  docs/checksums describe 140,855 (taxid 10407, HBV-only). Pick one.
- **D2 — MAFFT.** Needs `sudo apt install mafft` (≥ 7.487).

## 2. Definition of done

Task 3 is done when, on the canonical dataset:

1. `python3 preprocess.py` has produced `data/processed/hbv_sequences.fasta` and
   `hbv_metadata.csv`.
2. `python3 analyze_mutations.py` exits 0 and atomically writes all three outputs.
3. `results/mutations.csv` has the 11-column schema and plausible row counts.
4. `results/mutation_analysis.json` records this machine, this MAFFT version, and
   counts that reconcile with `mutation_analysis_skipped.csv`.
5. A documented validation pass (gates in §6) is recorded.
6. `TASK3.md`, `SOURCES.md`, `AGENTS.md` describe this actual run, not a foreign one.
7. No generated multi-GB artifact is committed; regeneration is documented.

## 3. Decisions required before execution

### D1 — canonical snapshot

| Option | Cost | Consequence |
|---|---|---|
| **A. Re-download taxid 10407 (HBV only)** — *recommended* | ~31 MB, ~1 min, then re-run preprocess | Matches committed docs, checksums and the prior 140,855-record results; no taxonomy filtering needed |
| B. Keep local taxid 10404 (142,254, family-wide) | free | Must regenerate all checksums and rewrite the dataset docs; ~1% non-HBV filtered by Task 2 |

Either way the reference `NC_003977.2` is present **(verified)**, so Task 3 is
unaffected by the choice; only the counts change.

### D2 — MAFFT

`sudo apt install mafft` → candidate 7.505 (≥ 7.487, supported).
Verify with `mafft --version`; the version is captured in the JSON report.

### D3 — commit policy for the mutation table

The full table is ~10.5 M rows / **~1.2 GB**, over GitHub's 100 MB file limit.
**Recommended:** commit the header-only `mutations.csv`, the JSON report and the
skipped CSV; keep the full table generated locally and document regeneration in
`TASK3.md`. Do not use `git add -f`.

## 4. Phase plan

### Phase 0 — preflight (no output changes)

```bash
cd /home/ant-caffien/Desktop/HBV-Genome-Pipeline
mafft --version                      # must be >= 7.487
du -sh data-collection/data/processed
df -h .                              # need ~2 GB free for the CSV
python3 -m py_compile analyze_mutations.py
PATH="$PWD/.venv/bin:$PATH" bash tests/run_tests.sh
```

Gate: all green. Fix environment before touching results.

### Phase 1 — resolve the dataset and run Task 2

```bash
# Option A
./download.sh --taxon 10407
(cd data-collection && sha256sum -c data/checksums/dataset.sha256)

# Task 2 (required Task 3 input)
PATH="$PWD/.venv/bin:$PATH" python3 preprocess.py
```

Gate: `data/processed/hbv_sequences.fasta`, `hbv_metadata.csv` and `manifest.json`
exist; `cleaning_report.json` records `input_fasta_records` / `retained_records`;
`NC_003977.2` survives preprocessing (`grep -c '^>NC_003977.2' data/processed/hbv_sequences.fasta` → 1).

### Phase 2 — provenance cleanup

- Quarantine the foreign run: move `results/mutation_analysis.json` and
  `results/mutation_analysis_skipped.csv` aside (do not keep them as evidence).
- Correct the stale claims:
  - `TASK3.md` "no biological mutation counts are claimed" → replace with the
    real run's validation note.
  - `SOURCES.md` "Batch size 1,000" → **500** (code default `--batch-size 500`).
  - `AGENTS.md` §2 "results/mutations.csv contains the Task 3 output schema" →
    state the actual state after the run.
- Keep the committed 0-byte `mutations.csv` header semantics.

### Phase 3 — smoke run on a small subset

Build a 200-record slice plus matching metadata in a temp dir, and run with a
separate output dir. This catches schema/parsing problems in minutes instead of
after a ~40-minute full run.

```bash
TMP=$(mktemp -d)
# keep the reference + first 200 non-reference records, and their metadata rows
python3 - "$TMP" <<'PY'
import csv, sys
from pathlib import Path
tmp = Path(sys.argv[1])
keep, write = [], False
for line in open("data/processed/hbv_sequences.fasta"):
    if line.startswith(">"):
        write = line.startswith(">NC_003977.2") or len(keep) < 201
        if write: keep.append(line)
    elif write:
        keep.append(line)
(tmp / "smoke.fasta").write_text("".join(keep))
PY
python3 analyze_mutations.py --input "$TMP/smoke.fasta" \
    --metadata data/processed/hbv_metadata.csv \
    --output "$TMP/mutations.csv" --skipped "$TMP/skipped.csv" \
    --report "$TMP/report.json" --batch-size 50
head -1 "$TMP/mutations.csv"
```

Gate: header matches the 11 columns, rows parse, JSON counts reconcile, exit 0.

### Phase 4 — full run

```bash
./run.sh --skip-download          # preprocess + analyze + statistics
# or, Task 3 only:
python3 analyze_mutations.py --batch-size 500 --threads 1
```

Run with `--threads 1` for a deterministic, reproducible artifact; use more
threads only after the first successful run. Expected ~35–45 min on this
hardware, 262 batches, memory bounded by one 500-sequence batch (safe under
~2.4 GB free).

### Phase 5 — validation gates

See §6. Record results in `TASK3.md`.

### Phase 6 — publish

- Update `report/report.md` methods/limitations with the real numbers.
- `SOURCES.md`: MAFFT version, batch size, run date.
- Commit code/docs/tests; leave the 1.2 GB CSV untracked per D3.
- Regenerate `data-collection/data/checksums/` if the snapshot changed (D1-A).

## 5. Runtime and resources

| Resource | Estimate | Basis |
|---|---|---|
| Wall clock | 35–45 min | prior run: 01:25:31 → 02:01:34 |
| Batches | ~262 | 130,837 eligible ÷ 500 |
| Peak memory | < 1 GB | streaming CSV writes, one batch in memory |
| Disk | ~1.2 GB CSV + ~0.6 GB skipped | 10.5 M rows × ~110 B |
| MAFFT threads | 1 (default) | determinism |

## 6. Validation gates (must all pass)

1. **Schema** — `results/mutations.csv` header equals
   `sample_id,reference_accession,position,reference_base,observed_base,variation_type,call_status,aligned_bases,percent_identity,completeness`.
2. **Reference integrity** — no row has `sample_id == NC_003977.2`; every
   `reference_accession == NC_003977.2`.
3. **Coordinates** — `1 <= position <= 3182`; insertions use `N^N+1`.
4. **Enumerations** — `variation_type ∈ {substitution,deletion,insertion,ambiguous}`,
   `call_status ∈ {definitive,uncertain}`.
5. **Arithmetic** — `analyzed + identity_skips + aligned_base_skips == eligible`;
   `input == eligible + below_length + other_skips + reference`.
6. **Reconciliation** — JSON `mutation_rows` equals `wc -l - 1` of the CSV;
   skipped CSV row count equals JSON `skipped_sequences`.
7. **Spot check** — manually re-align 3 sequences (`mafft --6merpair --addfragments`)
   and confirm their calls at 5 random positions.
8. **Determinism** — rerun the Phase 3 subset with the same settings; files are
   byte-identical.
9. **Failure behaviour** — point `--mafft` at a missing binary: exit non-zero,
   actionable message, previous outputs untouched.

## 7. Risks

| Risk | Mitigation |
|---|---|
| MAFFT segfault on a batch | built-in halving retry; 5 consecutive singleton failures abort loudly |
| Run dies mid-way | outputs are staged `.tmp` + `os.replace`, so the previous good result survives |
| Partial/short records dominate | documented thresholds (100 nt / 80 aligned / 70 % identity); record skip reasons |
| Circular genome | documented limitation; sequences rotated across the origin may not align |
| Reverse-complement records | documented limitation; may fail identity and be skipped |
| 1.2 GB artifact in git | D3: header-only commit + documented regeneration |

## 8. Acceptance checklist

- [ ] `mafft --version` ≥ 7.487
- [ ] Task 2 outputs present for the canonical snapshot
- [ ] smoke run passes all gates
- [ ] full run exits 0
- [ ] gates 1–9 pass and are recorded
- [ ] `TASK3.md` / `SOURCES.md` / `AGENTS.md` match the real run
- [ ] `bash tests/run_tests.sh` still green (99 tests)
- [ ] no >100 MB file staged for commit
