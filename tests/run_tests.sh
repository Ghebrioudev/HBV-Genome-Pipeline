#!/usr/bin/env bash
# Run the original Task 1 suite and the Task 2 preprocessing suite.
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

command -v python3 >/dev/null 2>&1 || {
  echo "[tests] ERROR: python3 is required" >&2
  exit 1
}

python3 -m unittest discover -s data-collection/tests -p 'test_*.py' -v
python3 -m unittest discover -s tests/preprocessing -p 'test_*.py' -v
