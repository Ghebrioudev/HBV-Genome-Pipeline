#!/usr/bin/env bash
#
# HBV data collection — default source is the NCBI Datasets virus API,
# with an optional local-file pathway.
#
# Usage:
#   ./download.sh                                  # NCBI, taxon 10404 (Hepadnaviridae)
#   ./download.sh --taxon 10407                    # NCBI, Hepatitis B virus only
#   ./download.sh --file path/to/hbv.fasta         # local dataset
#   ./download.sh --force                          # re-download even if valid
#   ./download.sh --help
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$BASE_DIR/src"
RAW_DIR="$BASE_DIR/data/raw"
PROC_DIR="$BASE_DIR/data/processed"
CHK_DIR="$BASE_DIR/data/checksums"

TAXON="${HBV_TAXON:-10404}"
INCLUDE="genome"
INPUT_FILE=""
FORCE=""
MIN_MB="100"
MAX_MB="1024"
API_KEY="${NCBI_API_KEY:-}"

info()  { printf '[download] %s\n' "$*" >&2; }
fail()  { printf '[download] ERROR: %s\n' "$*" >&2; exit 1; }

usage() {
  sed -n '3,14p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  cat <<'EOF'

Options:
  --taxon ID         NCBI taxonomy id/name (default: 10404 Hepadnaviridae)
  --include LIST     Comma separated sequence types (default: genome)
  --file PATH        Use a local FASTA/.gz/.zip instead of the NCBI API
  --out-dir DIR      Override the raw download directory
  --min-mb N         Minimum dataset size in MB (default: 100)
  --max-mb N         Maximum dataset size in MB (default: 1024)
  --api-key KEY      Optional NCBI API key (or set NCBI_API_KEY)
  --force            Re-download even when a valid package already exists
  -h, --help         Show this help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxon)      TAXON="${2:?--taxon needs a value}"; shift 2 ;;
    --include)    INCLUDE="${2:?--include needs a value}"; shift 2 ;;
    --file)       INPUT_FILE="${2:?--file needs a value}"; shift 2 ;;
    --out-dir)    RAW_DIR="${2:?--out-dir needs a value}"; shift 2 ;;
    --min-mb)     MIN_MB="${2:?--min-mb needs a value}"; shift 2 ;;
    --max-mb)     MAX_MB="${2:?--max-mb needs a value}"; shift 2 ;;
    --api-key)    API_KEY="${2:?--api-key needs a value}"; shift 2 ;;
    --force)      FORCE="--force"; shift ;;
    -h|--help)    usage; exit 0 ;;
    *)            fail "unknown argument: $1 (try --help)" ;;
  esac
done

# ---------------------------------------------------------------- dependencies
command -v python3 >/dev/null 2>&1 || fail "python3 is required but not installed"
python3 -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' \
  || fail "python3 >= 3.8 is required (found $(python3 -V 2>&1))"
[[ -f "$SRC_DIR/ncbi_client.py" ]]       || fail "missing $SRC_DIR/ncbi_client.py"
[[ -f "$SRC_DIR/dataset_validator.py" ]] || fail "missing $SRC_DIR/dataset_validator.py"
[[ -f "$SRC_DIR/file_loader.py" ]]       || fail "missing $SRC_DIR/file_loader.py"

mkdir -p "$RAW_DIR" "$PROC_DIR" "$CHK_DIR"

# ---------------------------------------------------------------------- source
if [[ -n "$INPUT_FILE" ]]; then
  [[ -f "$INPUT_FILE" ]] || fail "input file not found: $INPUT_FILE"
  info "source: local file ($INPUT_FILE)"

  info "step 1/3 — validating local dataset (size policy: warn)"
  python3 "$SRC_DIR/dataset_validator.py" \
    --dataset "$INPUT_FILE" \
    --min-mb "$MIN_MB" --max-mb "$MAX_MB" \
    --size-policy warn \
    --report-json "$CHK_DIR/local_validation.json" \
    || fail "local dataset failed validation"

  info "step 2/3 — normalising local dataset"
  python3 "$SRC_DIR/file_loader.py" \
    --input "$INPUT_FILE" \
    --out-dir "$PROC_DIR.staging" \
    --workdir "$RAW_DIR/local" \
    || fail "failed to normalise the local dataset"

  SOURCE_DESCRIPTION="local_file:$(basename "$INPUT_FILE")"
  ORIGINAL="$INPUT_FILE"
else
  info "source: NCBI Datasets virus API (taxon $TAXON)"

  info "step 1/4 — retrieving data package from NCBI"
  python3 "$SRC_DIR/ncbi_client.py" \
    --taxon "$TAXON" \
    --include "$INCLUDE" \
    --out-dir "$RAW_DIR" \
    ${API_KEY:+--api-key "$API_KEY"} \
    ${FORCE:+$FORCE} \
    || fail "NCBI download failed"

  info "step 2/4 — validating downloaded package (size policy: hard)"
  python3 "$SRC_DIR/dataset_validator.py" \
    --dataset "$RAW_DIR/extracted" \
    --min-mb "$MIN_MB" --max-mb "$MAX_MB" \
    --size-policy hard \
    --report-json "$CHK_DIR/ncbi_validation.json" \
    --checksums "$CHK_DIR/ncbi_extracted.sha256" \
    --checksums-base "$BASE_DIR" \
    || fail "downloaded dataset failed validation"

  info "step 3/4 — normalising NCBI package"
  python3 "$SRC_DIR/file_loader.py" \
    --package "$RAW_DIR/extracted" \
    --out-dir "$PROC_DIR.staging" \
    || fail "failed to normalise the NCBI package"

  SOURCE_DESCRIPTION="ncbi_datasets_virus:taxon=$TAXON"
  ORIGINAL="$RAW_DIR/ncbi_dataset.zip"
fi

# Promote the new processed output only after everything succeeded, so a failed
# run never destroys a previously good dataset.
info "step — promoting normalised dataset"
rm -rf "$PROC_DIR.old"
if [[ -d "$PROC_DIR" ]]; then mv "$PROC_DIR" "$PROC_DIR.old"; fi
mv "$PROC_DIR.staging" "$PROC_DIR"
rm -rf "$PROC_DIR.old"

# ------------------------------------------------------------------ checksums
# Paths are stored relative to the data-collection directory so the committed
# checksum files stay portable across machines.
info "step — writing checksums"
python3 - "$BASE_DIR" "$ORIGINAL" "$PROC_DIR" "$CHK_DIR/dataset.sha256" "$SOURCE_DESCRIPTION" <<'PY'
import hashlib
import os
import sys
from datetime import datetime, timezone

base, original, proc_dir, destination, source = sys.argv[1:6]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


files = [original]
files += sorted(
    os.path.join(proc_dir, name)
    for name in os.listdir(proc_dir)
    if os.path.isfile(os.path.join(proc_dir, name))
)

with open(destination, "w", encoding="utf-8") as handle:
    handle.write(f"# generated on {datetime.now(timezone.utc).isoformat(timespec='seconds')}\n")
    handle.write(f"# source: {source}\n")
    for path in files:
        handle.write(f"{sha256(path)}  {os.path.relpath(path, base)}\n")
PY

# -------------------------------------------------------------------- summary
read -r RECORDS PROC_MB SRC_MB <<<"$(python3 - "$PROC_DIR" "$ORIGINAL" <<'PY'
import json, os, sys

proc_dir, original = sys.argv[1], sys.argv[2]

def size_of(path):
    if os.path.isfile(path):
        return os.path.getsize(path)
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            total += os.path.getsize(os.path.join(root, name))
    return total

manifest = json.load(open(os.path.join(proc_dir, "manifest.json"), encoding="utf-8"))
print(manifest["records"], f"{size_of(proc_dir) / 1048576:.1f}", f"{size_of(original) / 1048576:.1f}")
PY
)"

cat >&2 <<EOF

[download] done
  source          : $SOURCE_DESCRIPTION
  sequences       : $RECORDS records
  processed data  : $PROC_DIR ($PROC_MB MB)
  checksums       : $CHK_DIR/dataset.sha256
  validation      : $CHK_DIR/*_validation.json
EOF

exit 0
