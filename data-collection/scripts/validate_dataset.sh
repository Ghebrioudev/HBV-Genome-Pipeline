#!/usr/bin/env bash
#
# Standalone dataset validation (Task B / integrity checks).
#
# Usage:
#   ./scripts/validate_dataset.sh --dataset path/to/hbv.fasta
#   ./scripts/validate_dataset.sh --dataset data/raw/extracted --strict
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATASET=""
POLICY="warn"
MIN_MB="100"
MAX_MB="1024"
REPORT=""

usage() {
  cat <<'EOF'
Options:
  --dataset PATH   File or directory to validate (required)
  --strict         Treat out-of-range size as a failure (hard policy)
  --min-mb N       Minimum size in MB (default: 100)
  --max-mb N       Maximum size in MB (default: 1024)
  --report PATH    Write the JSON validation report to this path
  -h, --help       Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset) DATASET="${2:?--dataset needs a value}"; shift 2 ;;
    --strict)  POLICY="hard"; shift ;;
    --min-mb)  MIN_MB="${2:?--min-mb needs a value}"; shift 2 ;;
    --max-mb)  MAX_MB="${2:?--max-mb needs a value}"; shift 2 ;;
    --report)  REPORT="${2:?--report needs a value}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[validate] ERROR: unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -n "$DATASET" ]] || { usage >&2; exit 2; }
command -v python3 >/dev/null 2>&1 || { echo "[validate] ERROR: python3 is required" >&2; exit 1; }

ARGS=(--dataset "$DATASET" --min-mb "$MIN_MB" --max-mb "$MAX_MB" --size-policy "$POLICY")
[[ -n "$REPORT" ]] && ARGS+=(--report-json "$REPORT")

exec python3 "$BASE_DIR/src/dataset_validator.py" "${ARGS[@]}"
