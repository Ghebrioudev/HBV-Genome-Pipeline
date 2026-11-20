#!/usr/bin/env bash
#
# scripts/lib.sh — reusable Bash functions shared by run.sh and the helper
# scripts in scripts/. This file is meant to be *sourced*, not executed:
#
#     source "$BASE_DIR/scripts/lib.sh"
#
# Functions:
#   log / warn / die            timestamped messages, die exits non-zero
#   require_file / require_cmd  clear error + exit code 3 / 4 when missing
#   require_python_module       checks an importable Python module
#   run_stage                   runs one pipeline stage (signal-aware)
#   verify_sha256_file          verifies a "sha256sum" style checksum list
#   verify_md5_file             verifies NCBI's md5sum.txt
#   write_sha256_manifest       computes SHA-256 checksums for a file list
#   set_readonly / set_writable protects / unlocks a directory tree
#   archive_results             dated .tar.gz archive, verified after creation
#
# Exit codes used across the project:
#   0 success · 1 generic failure · 2 bad command-line usage
#   3 missing file · 4 missing dependency · 5 integrity (checksum) failure
#   130 interrupted (SIGINT) · 143 terminated (SIGTERM)

# Guard against double sourcing.
[[ -n "${_HBV_LIB_LOADED:-}" ]] && return 0
_HBV_LIB_LOADED=1

readonly EXIT_USAGE=2
readonly EXIT_MISSING_FILE=3
readonly EXIT_MISSING_DEP=4
readonly EXIT_INTEGRITY=5

# PID of the stage currently running in the background (see run_stage).
CHILD_PID=""

_timestamp() { printf '%(%Y-%m-%d %H:%M:%S)T' -1; }

log()  { printf '[%s] [INFO]  %s\n' "$(_timestamp)" "$*" >&2; }
warn() { printf '[%s] [WARN]  %s\n' "$(_timestamp)" "$*" >&2; }

# die MESSAGE [EXIT_CODE] — print a clear error and stop.
die() {
    local message="$1" code="${2:-1}"
    printf '[%s] [ERROR] %s\n' "$(_timestamp)" "$message" >&2
    exit "$code"
}

# require_file PATH [HINT]
require_file() {
    local path="$1" hint="${2:-}"
    if [[ ! -f "$path" ]]; then
        die "required file does not exist: $path${hint:+ — $hint}" "$EXIT_MISSING_FILE"
    fi
}

# require_cmd COMMAND [HINT]
require_cmd() {
    local cmd="$1" hint="${2:-}"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        die "required command not found: $cmd${hint:+ — $hint}" "$EXIT_MISSING_DEP"
    fi
}

# require_python_module PYTHON MODULE [HINT]
require_python_module() {
    local python="$1" module="$2" hint="${3:-}"
    if ! "$python" -c "import $module" >/dev/null 2>&1; then
        die "Python module '$module' is not importable by $python${hint:+ — $hint}" "$EXIT_MISSING_DEP"
    fi
}

# run_stage NAME COMMAND [ARGS...]
# The command runs in the background and the shell waits for it, so that
# SIGINT/SIGTERM traps fire immediately instead of after the stage finishes.
run_stage() {
    local name="$1"
    shift
    local started=$SECONDS status=0
    log "==> START  $name"
    "$@" &
    CHILD_PID=$!
    wait "$CHILD_PID" || status=$?
    CHILD_PID=""
    if (( status != 0 )); then
        die "stage '$name' failed with exit code $status (see the messages above)" "$status"
    fi
    log "==> DONE   $name ($(( SECONDS - started ))s)"
}

# verify_sha256_file CHECKSUM_FILE BASE_DIR
# Comment lines (#...) are removed with a pipe before sha256sum -c.
verify_sha256_file() {
    local list="$1" base="$2" failed
    require_file "$list"
    log "verifying SHA-256 checksums from ${list#"$BASE_DIR"/}"
    failed="$(cd "$base" && grep -v '^#' "$list" | sha256sum -c - 2>&1 | grep -v ': OK$' || true)"
    if [[ -n "$failed" ]]; then
        printf '%s\n' "$failed" >&2
        die "checksum verification failed for $list" "$EXIT_INTEGRITY"
    fi
    log "  OK: $(grep -vc '^#' "$list") files match"
}

# verify_md5_file MD5SUM_FILE — NCBI ships md5sum.txt inside each package.
verify_md5_file() {
    local list="$1" failed
    require_file "$list"
    log "verifying NCBI MD5 checksums from ${list#"$BASE_DIR"/}"
    failed="$(cd "$(dirname "$list")" && md5sum -c "$(basename "$list")" 2>&1 | grep -v ': OK$' || true)"
    if [[ -n "$failed" ]]; then
        printf '%s\n' "$failed" >&2
        die "NCBI MD5 verification failed for $list" "$EXIT_INTEGRITY"
    fi
    log "  OK: $(wc -l < "$list") files match"
}

# write_sha256_manifest DESTINATION BASE_DIR FILE...
# Paths are stored relative to BASE_DIR so the list is portable.
write_sha256_manifest() {
    local destination="$1" base="$2"
    shift 2
    {
        printf '# generated on %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        (cd "$base" && sha256sum "$@")
    } > "$destination"
}

# set_readonly DIR — remove write permission from every file and directory.
set_readonly() {
    [[ -d "$1" ]] || return 0
    chmod -R a-w "$1"
}

# set_writable DIR — restore owner write permission (needed before a refresh).
set_writable() {
    [[ -d "$1" ]] || return 0
    chmod -R u+w "$1"
}

# archive_results ARCHIVE_DIR BASE_DIR PATH... — create a dated tar.gz from the
# given paths (relative to BASE_DIR), then list it back to prove it is readable.
# Prints the archive path on stdout.
archive_results() {
    local archive_dir="$1" base="$2"
    shift 2
    local archive
    archive="$archive_dir/hbv_results_$(date +%Y%m%d_%H%M%S).tar.gz"
    mkdir -p "$archive_dir"
    tar -czf "$archive.part" -C "$base" --exclude='results/mutations.csv' "$@"
    tar -tzf "$archive.part" >/dev/null || die "archive is unreadable: $archive.part" "$EXIT_INTEGRITY"
    mv "$archive.part" "$archive"
    printf '%s\n' "$archive"
}
