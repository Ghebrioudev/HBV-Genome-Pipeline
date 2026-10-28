# Test fixtures

`sample.fasta`, `sample.fasta.gz` and `sample.zip` contain the **first two records
of the real Phase 1 dataset** (`data/processed/sequences.fasta`), so they are real
HBV sequences rather than synthetic ones.

They exist only to exercise the local-file pathway and the validator without
touching the full 405 MB dataset.

```bash
# validation only (does not modify the processed dataset)
python3 src/dataset_validator.py --dataset tests/data/sample.fasta --min-mb 0 --size-policy warn
python3 src/dataset_validator.py --dataset tests/data/sample.fasta.gz --min-mb 0 --size-policy warn
python3 src/dataset_validator.py --dataset tests/data/sample.zip --min-mb 0 --size-policy warn
```

Regenerate from the real dataset:

```bash
awk 'BEGIN{n=0} /^>/{n++} n<=2{print}' data/processed/sequences.fasta > tests/data/sample.fasta
gzip -c tests/data/sample.fasta > tests/data/sample.fasta.gz
(cd tests/data && zip -q sample.zip sample.fasta)
```
