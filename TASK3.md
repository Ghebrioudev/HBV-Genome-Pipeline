# Task 3 — Sequence comparison and mutation detection

Task 3 consumes the cleaned HBV FASTA and metadata produced by Task 2. It does
not download data, repeat preprocessing or generate Task 4 statistics/figures.

## Reference and scope

The reference is **NC_003977.2**, the NCBI RefSeq genome for *Hepatitis B virus*
(strain ayw), length **3,182 nt**. The reference record is included in
`data/processed/hbv_sequences.fasta`; the script verifies that it exists and
contains only A/C/G/T before starting. Mutation positions use 1-based
coordinates on the ungapped reference sequence.

Reference record: <https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2>

## Alignment method

Task 3 uses **MAFFT 7** in reference-anchored `--addfragments` mode. The
reference is supplied as a one-record alignment and eligible sample sequences
are added as fragments. The default batch size is 500 sequences and the default
thread count is one. This bounds the resources required for the bundled dataset.

The script selects `--6merpair` rather than `--auto`. In the initial supplied
MAFFT 7.525 run, `--auto` selected the higher-cost `multipair`/`pairlocalalign`
path and the process exited with a segmentation fault during batch 2. The
immediate trigger was not established from the available evidence. The
corrected configuration uses `--6merpair`; if an alignment batch still fails,
Task 3 recursively splits it into smaller batches. A singleton that continues
to fail is recorded as skipped, and five consecutive singleton failures stop
the run as a likely system-wide alignment problem.

The command is equivalent to:

```bash
mafft --6merpair --addfragments fragments.fasta --thread 1 reference.fasta
```

No `--keeplength` option is used, so insertion columns can be retained. After
each alignment, removing gaps from the reference row must reproduce the
original reference exactly. Queries are compared to the reference only after
alignment, never by raw FASTA character index. The run report records the
MAFFT version, method, batch size and thread count.

### Install MAFFT

On Ubuntu/Debian:

```bash
sudo apt update
sudo apt install mafft
mafft --version
```

Use MAFFT **7.487 or newer**. If `mafft` is not on `PATH`, pass its executable
path using `--mafft /path/to/mafft`.

## Variant definitions and coordinates

- `substitution`: a concrete A/C/G/T base differs from the reference.
- `deletion`: a query alignment gap occurs between reference positions covered
  by observed sample bases.
- `insertion`: one or more sample bases occur in reference-gap columns; each
  consecutive run is represented as one event.
- `ambiguous`: an IUPAC ambiguity symbol prevents a definitive base call. Its
  `call_status` is `uncertain`.

Substitution and deletion positions are 1-based reference positions. An
insertion between reference positions 125 and 126 is written `125^126`.
For insertions, `reference_base` is `-`; for deletions, `observed_base` is `-`.
Each deleted reference base has its own row; a grouped insertion can contain
multiple letters in `observed_base`.

Partial sequence ends are treated as **uncovered sequence**, not deletions.
Only internal alignment gaps within the span supported by observed sample bases
are deletion calls. Terminal insertions are not called because both reference
flanks are not observed. Ambiguous bases may contribute to the covered span but
are excluded from the A/C/G/T identity calculation.

## Eligibility defaults

The default run skips sequences shorter than **100 nt**. After alignment, a
sequence must contain at least **80 concrete aligned bases** and have at least
**70% pairwise identity** to the reference. These are configurable analysis
thresholds, not biological constants. Partial records that meet the thresholds
remain eligible; full-genome coverage is not required. All skip decisions are
recorded in `results/mutation_analysis_skipped.csv`.

Run from the repository root after Task 2:

```bash
python3 analyze_mutations.py
```

Optional parameters:

```bash
python3 analyze_mutations.py \
  --min-sequence-length 100 \
  --min-aligned-bases 80 \
  --min-identity 70 \
  --batch-size 500 \
  --threads 1
```

Use `python3 analyze_mutations.py --help` for all supported options and custom
input/output paths. Task 3 reads only `data/processed/hbv_sequences.fasta` and
`data/processed/hbv_metadata.csv`; it does not modify Task 1 or Task 2 inputs.
Outputs are staged and atomically published after a successful run.

## Outputs

| Output | Purpose |
|---|---|
| `results/mutations.csv` | One row per variant event/base, including sample/reference accessions, position, bases, variation type, call certainty, aligned-base count, identity and completeness label |
| `results/mutation_analysis_skipped.csv` | Accession, sequence length, completeness and skip reason for each skipped record |
| `results/mutation_analysis.json` | Run provenance, reference, MAFFT version/parameters, thresholds and counts by variant and skip reason |

Task 4 reads the mutation table and creates `results/summary.csv`,
`results/mutation_types.csv`, `results/statistics_summary.json` and figures.
Those Task 4 files are not read or written by Task 3.

## Completed run recorded in this archive

The bundled run completed on **14 November 2026** using MAFFT **v7.525**,
500-sequence batches and one thread (262 batches). Its Task 2 input contained
140,889 records, including the reference. The run reported:

- **130,871 eligible sample sequences** after the 100-nt minimum-length rule.
- **128,406 analysed sequences** and **12,482 skipped records**. The skipped
  reasons were 10,017 below 100 nt, 2,464 below 70% identity and 1 with fewer
  than 80 concrete aligned bases after alignment.
- **10,567,649 mutation rows:** 10,098,062 definitive and 469,587 uncertain.
- Variant rows: 9,097,676 substitutions, 976,187 deletions, 24,199 insertions
  and 469,587 ambiguous calls.

The result counts reconcile with `results/mutation_analysis.json`,
`results/mutation_analysis_skipped.csv` and `results/mutations.csv` in the
archive. Task 4 uses the 128,406 analysed sequences as its frequency denominator
and excludes uncertain calls from confirmed-mutation frequencies.

## Tests

```bash
python3 -m unittest discover -s tests/task3 -p 'test_*.py' -v
```

The offline tests include known-answer alignment fixtures and a deterministic
MAFFT test double for output generation and batch-retry behavior. A real-MAFFT
integration test runs when MAFFT is installed and is skipped otherwise.

## Known limitations

- The reference is treated as a linear 1–3,182 coordinate path. HBV is circular;
  origin-spanning sequences may not align correctly.
- Input orientation is preserved; no reverse-complement search is performed.
  Opposite-orientation records may fail the identity threshold and be skipped.
- The 100-nt, 80-base and 70%-identity thresholds are heuristics and should be
  reviewed before biological interpretation.
- Alignment-dependent indel placement can remain uncertain in repetitive or
  low-complexity regions.
- The analysis reports nucleotide events only. It does not infer genotype,
  phylogeny, amino-acid effects, resistance, clinical significance or treatment
  implications.
- Partial sequences predominate in the input; mutation frequencies count
  accessions rather than unique sequence groups, so identical submissions can
  influence frequencies.

No clinical interpretation is made.
