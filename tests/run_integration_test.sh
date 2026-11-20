#!/usr/bin/env bash
# Clean-start integration test: runs the whole run.sh pipeline (with -s) on a
# tiny fixture. MAFFT is replaced by a deterministic stub, so this stays
# offline and fast. Also covers getopts, checksums, logs, archive, trap and
# restart behaviour.
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

if [[ -z "${PYTHON:-}" && -x "$BASE_DIR/.venv/bin/python" ]]; then
  PYTHON="$BASE_DIR/.venv/bin/python"
fi
PYTHON="${PYTHON:-python3}"

exec "$PYTHON" -m unittest discover -s tests/integration -p 'test_*.py' -v
