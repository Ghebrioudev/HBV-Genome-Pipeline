#!/usr/bin/env bash
#
# Run the Phase 1 test suite.
#
# No network access and no third-party packages are required.
#
# Usage:
#   ./tests/run_tests.sh            # all tests
#   ./tests/run_tests.sh -v         # verbose
#   ./tests/run_tests.sh test_file_loader   # a single module
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$BASE_DIR"

command -v python3 >/dev/null 2>&1 || { echo "[tests] ERROR: python3 is required" >&2; exit 1; }

if [[ $# -gt 0 && "$1" != -* ]]; then
  module="$1"; shift
  exec python3 -m unittest discover -s tests -p "${module}.py" -v "$@"
fi

exec python3 -m unittest discover -s tests -p 'test_*.py' -v "$@"
