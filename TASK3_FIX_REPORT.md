# Task 3 — Debugging and Fix Report

## What the screenshot shows

There were two separate issues:

1. The first command was typed with the wrong script name (`analyse_mutations.py`;
   the assignment text also refers to `analyze_mutation.py`). The actual file is
   `analyze_mutations.py` (American spelling, plural `mutations`).
2. After the correct script was executed, MAFFT 7.525 exited with a segmentation
   fault in its native `pairlocalalign` routine while processing batch 2. This is
   an external MAFFT process crash, not a Python syntax error.

The MAFFT log visible in the screenshot says that `--auto` selected the
`multipair` strategy (`use multipair`) and that MAFFT crashed at `pairlocalalign`.
The official MAFFT documentation describes `--multipair` as high-cost in time
and memory and `--6merpair` as a low-cost option for `--addfragments`.

The screenshot is not enough to prove whether the immediate trigger was memory
pressure, a particular sequence, or a native MAFFT defect. The correction avoids
the expensive path indicated in the log and adds recovery for batch-specific
failures rather than claiming a more specific cause than the evidence supports.

## Why the CSV files only had headers

The project ZIP contained header-only CSV templates and no
`results/mutation_analysis.json` completion report. The script writes to
temporary files and replaces the final outputs only after all batches finish.
When MAFFT failed during batch 2, the temporary files were deleted and the
pre-existing header-only templates remained. This means the files were not the
result of a completed zero-mutation analysis.

The updated script removes header-only Task 3 templates when there is no
completion JSON report. A truly completed run that finds zero mutation rows will
still preserve its CSV header and will also have a completion report.

## Corrections

- Replaced `--auto --addfragments` with `--6merpair --addfragments`.
- Set the default batch size to 500 (from 1,000) to reduce resource pressure.
- Added adaptive retries: if a batch fails, it is split into two smaller batches
  and retried. A single record that still fails is listed with reason
  `mafft_alignment_failed`; five consecutive single-record failures stop the
  run as a likely system-wide issue.
- Added cleanup of unverified, header-only output templates when no successful
  JSON report exists.
- Preserved the reference, 1-based reference coordinates, sequence identifiers,
  partial-coverage logic, IUPAC uncertainty calls, and CSV schema.

## Dataset inspection

The supplied Task 2 FASTA was parsed without changing it:

- FASTA records: 140,855 (including the `NC_003977.2` reference).
- Reference records: 1; reference length: 3,182 nt.
- Non-reference records at least 100 nt long: 130,837 eligible by the initial
  length rule.
- Non-reference records shorter than 100 nt: 10,017; these are skipped by the
  existing documented minimum-length rule.

These are input/eligibility counts, not completed alignment counts. The user’s
screenshot confirms that batch 1 finished and reported 980 analyzed sequences
before batch 2 crashed. The remaining dataset was not completed.

## Validation actually performed in this environment

The complete test runner was executed after the changes:

- Task 1: 57 passed.
- Task 2: 20 passed.
- Task 3: 17 passed, 1 skipped.
- Total: 94 passed, 1 skipped.

The skipped test is the real-MAFFT integration test. MAFFT is not installed in
this build environment, and installing it here was not successful. The new
Task 3 tests use a fake MAFFT command to verify output generation, the required
`--6merpair` argument, and automatic batch splitting when a simulated MAFFT
segmentation fault occurs. These tests validate the Python recovery logic, but
they do not substitute for running MAFFT against the full biological dataset.

Therefore no full-dataset mutation count is claimed in this fix report. The
actual counts will be written to `results/mutation_analysis.json` only after a
successful run with real MAFFT.

## Run on Ubuntu

From the project root:

```bash
mafft --version
python3 analyze_mutations.py
```

The corrected default command invokes `mafft --6merpair --addfragments` in
batches of 500. You can inspect progress in the terminal. At the end, check that
`results/mutation_analysis.json` exists and reports a completed run before
interpreting either CSV. The outputs are:

- `results/mutations.csv`
- `results/mutation_analysis_skipped.csv`
- `results/mutation_analysis.json`

Task 1 and Task 2 implementations and data were not edited.
