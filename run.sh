#!/usr/bin/env bash
#
# Run the HBV pipeline stages in their documented order.
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Prefer a project-local virtual environment when present, so the Phase 4
# dependencies are found without the caller having to activate it.
if [[ -z "${PYTHON:-}" && -x "$BASE_DIR/.venv/bin/python" ]]; then
    PYTHON="$BASE_DIR/.venv/bin/python"
else
    PYTHON="${PYTHON:-python3}"
fi
SKIP_DOWNLOAD=0
CLEAN_RESULTS=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [--skip-download] [--clean]

Options:
  --skip-download  Start with the data already present in data/.
  --clean          Remove generated result files before running.
EOF
}

for argument in "$@"; do
    case "$argument" in
        --skip-download)
            SKIP_DOWNLOAD=1
            ;;
        --clean)
            CLEAN_RESULTS=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            printf 'Error: unknown option: %s\n\n' "$argument" >&2
            usage >&2
            exit 2
            ;;
    esac
done

require_file() {
    if [[ ! -f "$BASE_DIR/$1" ]]; then
        printf 'Error: required pipeline file is missing: %s\n' "$1" >&2
        exit 1
    fi
}

run_stage() {
    local name="$1"
    shift
    printf '\n==> %s\n' "$name"
    "$@"
}

require_file "download.sh"
require_file "preprocess.py"
require_file "analyze_mutations.py"
require_file "statistics.py"

if ! command -v "$PYTHON" >/dev/null 2>&1; then
    printf 'Error: Python interpreter not found: %s\n' "$PYTHON" >&2
    exit 1
fi

if (( CLEAN_RESULTS )); then
    rm -f "$BASE_DIR/results/mutations.csv" "$BASE_DIR/results/summary.csv"
    find "$BASE_DIR/results/figures" -mindepth 1 -maxdepth 1 -type f -delete
fi

mkdir -p "$BASE_DIR/data" "$BASE_DIR/results/figures" "$BASE_DIR/report"

if (( ! SKIP_DOWNLOAD )); then
    run_stage "Download data" bash "$BASE_DIR/download.sh"
else
    printf '\n==> Download data (skipped)\n'
fi

run_stage "Preprocess sequences" "$PYTHON" "$BASE_DIR/preprocess.py"
run_stage "Analyze mutations" "$PYTHON" "$BASE_DIR/analyze_mutations.py"
run_stage "Generate statistics" "$PYTHON" "$BASE_DIR/statistics.py"

printf '\nPipeline completed successfully.\n'