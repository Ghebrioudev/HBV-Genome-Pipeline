# HBV Genome Pipeline — Project Report

## 1. Objective

The project implements a reproducible workflow for acquiring HBV-related
nucleotide records, preparing a validated HBV dataset, identifying
reference-relative nucleotide variation, and summarising the resulting calls.
The analysis is based on the NCBI RefSeq genome `NC_003977.2` and is intended as
a bioinformatics data-processing and descriptive analysis workflow. It does not
provide clinical or treatment recommendations.

## 2. Methods and implementation

### Data collection and preprocessing

Task 1 retrieves a nucleotide genome package from the NCBI Datasets Virus v2
API. The bundled snapshot was retrieved on 14 November 2026 for taxon 10404
(*Hepadnaviridae*) and contained 142,288 records (132,005,421 nucleotide
positions). Task 2 verifies the Task 1 SHA-256 manifest, joins FASTA and metadata
by versioned accession, and filters records using HBV taxonomy and lineage
information. It preserves valid partial sequences and IUPAC ambiguity symbols,
records exclusions and validation decisions, and distinguishes repeated
accessions from different accessions with identical sequence strings.

### Alignment and mutation calling

Task 3 aligns eligible sample fragments against the 3,182-nt `NC_003977.2`
reference using MAFFT v7.525 with `--6merpair --addfragments`, 500-record batches
and one thread. It identifies substitutions, insertions, deletions and
uncertain ambiguous-base calls in 1-based reference coordinates. Partial
terminal coverage is treated as unobserved sequence rather than a deletion.
The default filters require a minimum sequence length of 100 nt, at least 80
concrete aligned bases and at least 70% pairwise identity.

### Statistics and reporting

Task 4 streams the mutation table in chunks, calculates per-mutation sequence
frequencies and mutation-type counts, and generates three figures. Ambiguous
calls are retained in the event table but excluded from confirmed-mutation
frequency statistics. Task 5 generates a self-contained HTML dashboard from the
machine-readable run summaries, CSV tables and figures.

## 3. Results

Task 2 excluded 1,399 non-HBV hepadnavirus records and retained **140,889 HBV
records**, representing 129,414,793 nucleotide positions. The retained set
contains 127,420 partial and 13,469 complete records. It includes 508,170 valid
IUPAC ambiguity positions across 18,285 records (0.392668% of retained
nucleotide positions). Exact sequence comparison identified 11,913 groups
containing 43,203 records; the 31,290 additional members were retained because
they have different accessions.

In the completed Task 3 run, 130,871 sample records met the initial 100-nt
length threshold. **128,406 sequences were analysed**; 12,482 were skipped
(10,017 below the length threshold, 2,464 below the 70% identity threshold and
one with insufficient aligned bases). The run produced **10,567,649 mutation
events**, comprising 10,098,062 definitive calls and 469,587 uncertain calls.
The event types were 9,097,676 substitutions, 976,187 deletions, 24,199
insertions and 469,587 ambiguous calls.

Task 4 used 128,406 analysed sequences as its frequency denominator and
identified **16,235 distinct confirmed mutations**. The five most frequent
substitutions were 530:T>C (61.5%), 535:A>C (57.2%), 520:C>A (56.0%),
576:C>A (45.6%) and 522:G>A (45.1%). Frequencies count accessions rather than
unique sequence strings.

## 4. Testing and reproducibility

The complete offline test runner executed 103 tests across five suites: **102 passed and one was
skipped**. The skipped test is the real-MAFFT integration check because a MAFFT
executable was not installed in the independent review environment. The other
Task 3 tests use a deterministic test double and include known-answer alignment
fixtures and simulated batch-failure recovery. Task 1, Task 2, Task 4 and
pipeline integration tests all passed.

The archive contains a completed full-dataset result record identifying MAFFT
v7.525 and the run parameters. The full biological analysis was not rerun during
this final review because the review environment did not have MAFFT installed.
The supplied results were checked for row-count and summary reconciliation,
and the bundled input/output manifests were checked against their file hashes.

The pipeline can be reproduced from the repository root after installing
MAFFT and the Python environment (`bash scripts/setup_env.sh`):

```bash
./run.sh -s     # reuse the downloaded snapshot (checksums verified)
./run.sh        # or: download from NCBI first
```

`run.sh` verifies the data checksums, reruns preprocessing, alignment and
statistics, regenerates the HTML dashboard, checksums the results and writes
a dated archive in `archives/`.

## 5. Limitations and conclusion

The dataset is dominated by partial sequences. Task 3 uses explicit heuristic
thresholds, keeps input orientation, and maps calls to a linear reference
coordinate system even though HBV is circular; origin-spanning or reverse-
complement records may therefore be misaligned or skipped. Repeated identical
sequences under different accessions are retained, so accession-based
frequencies may be influenced by repeated submissions. The project does not
perform genotype, phylogenetic, codon, amino-acid, resistance or clinical
interpretation because the required validated annotations are not available
in the supplied inputs.

The completed workflow provides auditable dataset provenance, validated
preprocessing, reference-relative mutation calls, descriptive statistics and
reproducible report artifacts. Its outputs should be interpreted within the
stated data-quality and alignment limitations.
