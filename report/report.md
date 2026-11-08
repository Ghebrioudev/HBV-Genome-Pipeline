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

The stage scripts are currently placeholders while the other project members
complete their work. Therefore, this integration layer verifies command
ordering, path handling, error propagation, and clean-start behavior; it does
not claim that mutation tables, summary statistics, figures, or clinical
interpretations have been produced. Those outputs must be validated after the
corresponding analysis stages are implemented.