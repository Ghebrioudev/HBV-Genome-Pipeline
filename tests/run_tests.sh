#!/usr/bin/env bash
# Run the offline test suites for Tasks 1–4 and pipeline integration.
# Uses .venv/bin/python when present (created by scripts/setup_env.sh),
# otherwise $PYTHON or python3.
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

if [[ -z "${PYTHON:-}" && -x "$BASE_DIR/.venv/bin/python" ]]; then
  PYTHON="$BASE_DIR/.venv/bin/python"
fi
PYTHON="${PYTHON:-python3}"
export PYTHON

command -v "$PYTHON" >/dev/null 2>&1 || {
  echo "[tests] ERROR: Python interpreter not found: $PYTHON" >&2
  exit 1
}

"$PYTHON" -m unittest discover -s data-collection/tests -p 'test_*.py' -v
"$PYTHON" -m unittest discover -s tests/preprocessing -p 'test_*.py' -v
"$PYTHON" -m unittest discover -s tests/task3 -p 'test_*.py' -v
"$PYTHON" -m unittest discover -s tests/task4 -p 'test_*.py' -v
"$PYTHON" -m unittest discover -s tests/integration -p 'test_*.py' -v
