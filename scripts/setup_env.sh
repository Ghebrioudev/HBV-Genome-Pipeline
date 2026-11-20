#!/usr/bin/env bash
#
# scripts/setup_env.sh — create the isolated Python environment (.venv) from
# requirements.txt and check the system tools needed by run.sh.
#
# Usage: bash scripts/setup_env.sh [-p PYTHON] [-h]
#   -p PYTHON   interpreter used to create .venv (default: python3, >= 3.11)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"

BOOTSTRAP_PYTHON="python3"
while getopts ":p:h" opt; do
    case "$opt" in
        p) BOOTSTRAP_PYTHON="$OPTARG" ;;
        h) sed -n '3,7p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        :) die "option -$OPTARG requires a value" "$EXIT_USAGE" ;;
        \?) die "unknown option: -$OPTARG (see -h)" "$EXIT_USAGE" ;;
    esac
done

require_cmd "$BOOTSTRAP_PYTHON" "install Python >= 3.11"
"$BOOTSTRAP_PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)' \
    || die "Python >= 3.11 is required (found $("$BOOTSTRAP_PYTHON" -V 2>&1))" "$EXIT_MISSING_DEP"
require_file "$BASE_DIR/requirements.txt"

if [[ ! -x "$BASE_DIR/.venv/bin/python" ]]; then
    log "creating virtual environment in .venv/"
    "$BOOTSTRAP_PYTHON" -m venv "$BASE_DIR/.venv" \
        || die "python3 -m venv failed (Ubuntu/Debian: sudo apt install python3-venv)" "$EXIT_MISSING_DEP"
fi
log "installing requirements.txt"
"$BASE_DIR/.venv/bin/python" -m pip install --quiet --upgrade pip
"$BASE_DIR/.venv/bin/python" -m pip install --quiet -r "$BASE_DIR/requirements.txt"

missing=0
for cmd in mafft sha256sum md5sum tar gzip awk sed grep cut sort uniq join; do
    if command -v "$cmd" >/dev/null 2>&1; then
        log "found $cmd"
    else
        warn "missing system command: $cmd"
        missing=1
    fi
done
(( missing == 0 )) || die "install the missing commands (Ubuntu/Debian: sudo apt install mafft coreutils gawk)" "$EXIT_MISSING_DEP"
log "environment ready: $("$BASE_DIR/.venv/bin/python" -V), $(mafft --version 2>&1 | head -n 1)"
log "run.sh uses .venv/bin/python automatically"
