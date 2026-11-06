# Tests

Phase 1 test suite — standard library `unittest` only, fully offline.

```bash
./run_tests.sh                    # all tests
./run_tests.sh -v
./run_tests.sh test_ncbi_client   # one module
```

| Module | What it protects |
|---|---|
| `test_dataset_validator.py` | FASTA / `.gz` / `.zip` parsing, IUPAC validation, empty records, the 100 MB – 1 GB size policy (`hard` vs `warn`), SHA-256 and portable checksum paths |
| `test_file_loader.py` | local `.fasta`/`.gz`/`.zip` → `sequences.fasta` + `metadata.csv`; NCBI report mapping (camelCase and snake_case); column order; missing-package errors |
| `test_ncbi_client.py` | request payload shape, ZIP validation, safe extraction, path-traversal rejection, retry/backoff on 429, fail-fast on 404, atomic `.part` handling. `urllib` is mocked — **no network** |
| `test_shell_scripts.py` | `bash -n`, executable bits, `--help`, unknown arguments, missing input, and that a failed `download.sh` never modifies `data/processed/` |

## Fixtures

`data/sample.fasta`, `data/sample.fasta.gz` and `data/sample.zip` contain the
**first two records of the real Phase 1 dataset**, so they are genuine HBV
sequences, not synthetic ones.

```bash
# validate without replacing the processed dataset
python3 ../src/dataset_validator.py --dataset data/sample.fasta --min-mb 0 --size-policy warn
```

Regenerate from the full dataset:

```bash
awk 'BEGIN{n=0} /^>/{n++} n<=2{print}' ../data/processed/sequences.fasta > data/sample.fasta
gzip -c data/sample.fasta > data/sample.fasta.gz
(cd data && zip -q sample.zip sample.fasta)
```

## Adding tests

- Put `test_*.py` next to these files; discovery picks them up automatically.
- Import the modules under test through `helpers.py`, which puts `../src` on
  `sys.path`.
- Never perform network I/O; mock `ncbi_client.urllib.request.urlopen`.
- Never write into `data/processed/`; use `tempfile.TemporaryDirectory()`.
