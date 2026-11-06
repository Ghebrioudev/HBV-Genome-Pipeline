# Task 3 — Sequence Comparison and Mutation Detection

Task 3 consumes the cleaned HBV FASTA and metadata produced by Task 2. It does
not download data, repeat preprocessing, generate statistics/figures, or write
the final project report.

## Reference and scope

The reference is **NC_003977.2**, the NCBI RefSeq genome for *Hepatitis B virus*
(strain ayw), length **3,182 nt**. This exact accession/version is already
included in `data/processed/hbv_sequences.fasta`; the script verifies that it
exists and contains only A/C/G/T before starting. The reference coordinates in
the mutation table are 1-based positions on the ungapped reference sequence.

Reference record: <https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2>

## Alignment method

The script uses **MAFFT 7** and its reference-anchored `--addfragments` mode.
The reference is supplied as the existing one-record alignment, and eligible
sample sequences are added as fragments. MAFFT is run in batches (500
sequences by default) to bound the resources required for a 140k-record dataset.
The script explicitly selects `--6merpair` instead of `--auto`: the supplied
MAFFT 7.525 run showed `--auto` selecting the expensive `multipair`/`pairlocalalign`
path and then crashing with a segmentation fault in batch 2. MAFFT documents
`--multipair` as high-cost in time and memory and `--6merpair` as its low-cost
alternative for fragment addition. If an alignment batch still fails, Task 3
recursively splits it into smaller batches; a single record that continues to
fail is reported in the skipped-record CSV. Five consecutive singleton failures
stop the run as a likely system-wide alignment problem.
The reference row is checked after every alignment: removing alignment gaps
must reproduce the original reference sequence exactly. Query and reference
rows are then compared only after alignment, never by raw FASTA character
index.

The invoked command is equivalent to:

```bash
mafft --6merpair --addfragments fragments.fasta --thread 1 reference.fasta
```

`--6merpair` uses a lower-cost k-mer distance strategy for fragment addition.
No `--keeplength` option is used, so MAFFT can retain insertion columns
rather than silently dropping sample insertions. The program records the actual
MAFFT version, CLI method, batch size and thread count in
`results/mutation_analysis.json`.

### Install MAFFT on Ubuntu

```bash
sudo apt update
sudo apt install mafft
mafft --version
```

Use MAFFT **7.487 or newer**; the upstream MAFFT documentation warns about a
`--auto` issue in releases 7.463–7.486. The upstream software page currently
lists 7.526. If using a package manager with an older package, install a newer
MAFFT build using the official Linux instructions:
<https://mafft.cbrc.jp/alignment/software/linux.html>.

No Python packages beyond the standard library are required. If `mafft` is
not on `PATH`, pass its executable path using `--mafft /path/to/mafft`.

## How variants are identified

The script maps aligned columns back to the ungapped reference coordinate
system and emits:

- `substitution`: a concrete A/C/G/T base differs from the reference.
- `deletion`: a query alignment gap occurs between reference positions that
  are actually covered by sample bases.
- `insertion`: one or more sample bases occur in reference-gap columns. A run
  of inserted bases is represented as one event.
- `ambiguous`: an IUPAC ambiguity symbol prevents a definitive base call. Its
  `call_status` is `uncertain`, never `definitive`.

For substitutions and deletions, `position` is the 1-based reference position.
An insertion between reference positions 125 and 126 is written as `125^126`.
The `reference_base` is `-` for an insertion, and `observed_base` is `-` for a
deletion. Each deleted reference base has its own row. A grouped insertion can
have multiple letters in `observed_base`.

Partial sequence ends are treated as **uncovered sequence**, not deletions.
Only internal alignment gaps within the span supported by observed sample bases
are deletion calls. Terminal insertions are not called because both reference
flanks are not observed. Ambiguous bases may contribute to the covered span,
but they are excluded from the A/C/G/T identity calculation.

## Eligibility defaults

To avoid making unreliable placements from extremely short fragments, the
default run skips sequences shorter than 100 nt. After alignment, a sequence
must have at least 80 concrete aligned base pairs and at least 70% pairwise
identity against the reference to pass the default callability filter. These
are explicit, configurable analysis thresholds, not biological constants.
All skip decisions are retained in `results/mutation_analysis_skipped.csv`.
Partial records of sufficient length are still eligible; they are not required
to cover the whole genome.

## Run Task 3

Run from the project root after Task 2 has completed:

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

Use `python3 analyze_mutations.py --help` to see all options and custom input/output
paths. Task 3 reads only the Task 2 files `data/processed/hbv_sequences.fasta`
and `data/processed/hbv_metadata.csv`; it does not modify them or any Task 1
files. Each output file is atomically replaced only after the complete run
succeeds.

## Outputs

| Output | Purpose |
|---|---|
| `results/mutations.csv` | One row per variant event/base, including sample and reference accessions, position, reference/observed bases, variation type, call certainty, aligned-base count, identity and Task 2 completeness label. |
| `results/mutation_analysis_skipped.csv` | One row per record skipped, with accession, sequence length, completeness and reason. |
| `results/mutation_analysis.json` | Run provenance, reference, MAFFT version/parameters, thresholds and counts by variation type and skip reason. |

The Task 4 placeholder `results/summary.csv` is not read or written by Task 3.
No genotype, amino-acid effect, clinical significance, statistics, charts or
interpretation is produced here.

## Tests

Run the offline Task 3 unit tests:

```bash
python3 -m unittest discover -s tests/task3 -p 'test_*.py' -v
```

The tests use a tiny MAFFT executable test double for the end-to-end output
contract, plus direct, known-answer aligned fixtures for substitutions,
insertions, deletions, partial coverage, ambiguity, thresholds and invalid
inputs. A real-MAFFT integration check runs automatically when MAFFT is
installed and is skipped otherwise.

**Validation note for this delivered archive:** MAFFT was not installed in the
build environment. The offline tests passed, but the 140,855-record project
dataset was not aligned here and no biological mutation counts are claimed.
Install MAFFT and execute the command above to create the real CSV and JSON
analysis outputs.

## Known limitations

- The reference is treated as a linear 1..3182 coordinate path. HBV is circular;
  a sequence rotated across the reference origin may not align correctly.
- The script preserves input orientation; it does not search reverse complements.
  Opposite-orientation records may fail the identity threshold and be skipped.
- The 100-nt, 80-aligned-base and 70%-identity defaults are transparent
  heuristics. They should be reviewed with the project supervisor/biologist
  before interpreting biological results; CLI options can change them.
- MAFFT alignment is not a phylogenetic, genotype, resistance, protein-effect or
  clinical analysis. Indel placement in repetitive/low-complexity regions can
  remain alignment-dependent.
- The supplied snapshot includes many short and partial records. A full run may
  require substantial CPU time; batch size and thread count are configurable.

## Failed and incomplete runs

The script publishes the CSV and JSON outputs only after a complete run. If a
run fails, non-empty prior results are left untouched. Header-only CSV templates
without a matching `mutation_analysis.json` are removed at the start of a new
valid run, because they are not evidence of a completed analysis. A successful
run that genuinely finds zero rows will retain its header-only CSV together with
a completed JSON report.
