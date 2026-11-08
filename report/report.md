# HBV Genome Pipeline Report

## Methods

The pipeline is orchestrated by `run.sh`, which executes the stages in this
order:

1. Data download (`download.sh`)
2. Sequence preprocessing (`preprocess.py`)
3. Mutation analysis (`analyze_mutations.py`)
4. Statistics and visualization (`statistics.py`)

The orchestrator uses repository-relative paths, creates the expected output
directories, and stops on the first failed command. `--skip-download` supports
offline reruns using data already present in `data/`, while `--clean` removes
generated CSV and figure files before execution.

## Results

The canonical run used the 142,254-record family-wide NCBI snapshot, of which
Task 2 retained **140,855 Hepatitis B virus records** (1,399 non-HBV
hepadnaviruses excluded). Task 3 aligned 130,837 eligible records (≥100 nt) to
`NC_003977.2` with MAFFT 7.525; **128,372 were analysed** and 12,482 were
skipped (10,017 below the 100-nt minimum, 2,464 below 70 % identity, 1 with
insufficient aligned bases).

The comparison produced **10,562,861 variant events** (10,093,413 definitive,
469,448 uncertain). Uncertain IUPAC calls are excluded from the frequency
statistics, leaving 16,215 distinct mutations across 128,372 analysed
sequences; 128,032 sequences (99.7 %) carry at least one definitive variant.

| Variant type | Calls | Share |
|---|---|---|
| Substitution | 9,093,190 | 86.1 % |
| Deletion | 976,187 | 9.2 % |
| Insertion | 24,036 | 0.2 % |
| Ambiguous (uncertain) | 469,448 | 4.4 % |

The most frequent substitutions are at reference positions 530 (61.5 % of
sequences), 535 (57.2 %), 520 (56.0 %), 576 (45.6 %) and 522 (45.1 %). The top
five cluster in the 520–576 window, and further frequent calls occur at 192,
292, 632 and 659. The mutation table carries no gene or codon annotation, so no
amino-acid change or protein effect is assigned here.

Figures in `results/figures/`: `top_mutations.png` (most frequent mutations),
`mutations_along_genome.png` (calls per 50-nt bin) and `mutation_types.png`.

**Genotype comparison was not performed.** The Task 1 metadata has no genotype
column (it is not normalised by NCBI), so the conditional Task 4 comparison is
skipped and `genotype_comparison.csv` is not produced.
## Reproducibility

Run the pipeline from the repository root with:

```bash
./run.sh
```

To test the integration layer without downloading data:

```bash
bash tests/run_integration_test.sh
```

The integration test copies the project to a temporary directory, runs a clean
offline pipeline execution, and verifies that invalid options fail.

## Current limitations

- **Partial sequences dominate.** 127,395 of 140,855 records are flagged
  `PARTIAL`; terminal unobserved sequence is not called as deletion.
- **Duplicate sequences are retained.** 31,288 records share an identical
  sequence with another accession (11,912 groups). Frequencies count sequences,
  not unique genomes, so closely related or redundant deposits can inflate them.
- **Thresholds are heuristics.** The 100-nt minimum, 80 aligned bases and 70 %
  identity defaults are analysis choices, not biological constants.
- **Linear reference.** HBV is circular; a sequence rotated across the origin
  may not align correctly. Input orientation is preserved, so
  reverse-complemented records may fail the identity threshold.
- **No functional annotation.** The table records nucleotide events only; it
  does not assign genes, codons or amino-acid changes.
- **No genotype comparison.** The NCBI metadata does not carry a normalised
  genotype column.
- **Uncertain calls are excluded** from frequencies (469,448 IUPAC-ambiguous
  events) but remain in the mutation table.
- **The integration test is orchestration-only**: it substitutes `PYTHON=true`,
  so it verifies stage ordering and option handling, not the numeric stages.

No clinical or resistance interpretation is made.