#!/usr/bin/env bash
# Clean-start integration test: runs preprocess -> analyze -> statistics on a
# tiny fixture through run.sh. MAFFT is replaced by a deterministic stub, so
# this stays offline and fast.
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

exec python3 -m unittest discover -s tests/integration -p 'test_*.py' -v
