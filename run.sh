#!/usr/bin/env bash
#
# run.sh — single entry point of the HBV Genome Pipeline.
# Runs every stage, from the raw NCBI download to the archived results:
#
#   1  download        NCBI Datasets virus package (data-collection/download.sh)
#   2  integrity       NCBI MD5 + project SHA-256 verification, raw data -> read-only
#   3  shell profile   grep/sed/awk/cut/sort/uniq pipelines on the sequences
#   4  preprocess      HBV taxonomy filter and sequence validation (preprocess.py)
#   5  mutations       MAFFT alignment to NC_003977.2 + mutation calls (analyze_mutations.py)
#   6  statistics      frequencies, CSV/JSON tables and figures (statistics.py)
#   7  shell checks    shell recount of the top mutations vs statistics.py
#   8  dashboard       self-contained HTML dashboard (report/make_dashboard.py)
#   9  checksums       SHA-256 of all results, written then verified
#  10  archive         dated tar.gz of the results in archives/
#
# Every run writes a timestamped log to logs/run_YYYYMMDD_HHMMSS.log.
#
set -euo pipefail

BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="$BASE_DIR/scripts"
if [[ ! -f "$SCRIPTS_DIR/lib.sh" ]]; then
    printf 'ERROR: required file does not exist: %s\n' "$SCRIPTS_DIR/lib.sh" >&2
    exit 3
fi
# shellcheck source=scripts/lib.sh
source "$SCRIPTS_DIR/lib.sh"

# ------------------------------------------------------------------ paths
DC_DIR="$BASE_DIR/data-collection"
RAW_DIR="$DC_DIR/data/raw"
TASK1_DIR="$DC_DIR/data/processed"
CHECKSUM_DIR="$DC_DIR/data/checksums"
CLEAN_DIR="$BASE_DIR/data/processed"
RESULTS_DIR="$BASE_DIR/results"
LOG_DIR="$BASE_DIR/logs"
ARCHIVE_DIR="$BASE_DIR/archives"
PID_FILE="$LOG_DIR/.run.pid"

# ---------------------------------------------------------------- options
SKIP_DOWNLOAD=0
CLEAN_RESULTS=0
FORCE_DOWNLOAD=0
NO_ARCHIVE=0
TAXON="${HBV_TAXON:-10404}"
THREADS=1
TOP=20

usage() {
    cat <<EOF
Usage: $(basename "$0") [-h] [-s] [-c] [-f] [-n] [-t TAXON] [-j THREADS] [-k TOP]

Runs the complete HBV Genome Pipeline (download -> integrity -> preprocessing
-> MAFFT mutation detection -> statistics/figures -> dashboard -> archive).

Options:
  -h            Show this help and exit.
  -s            Skip the NCBI download and reuse the data already in
                data-collection/data/ (its checksums are still verified).
  -c            Clean previously generated results before running.
  -f            Force a fresh NCBI download even if a valid package exists.
  -n            Do not create the dated tar archive at the end.
  -t TAXON      NCBI taxonomy id to download (default: $TAXON, Hepadnaviridae;
                Hepatitis B virus records are selected afterwards).
  -j THREADS    MAFFT threads per alignment batch (default: $THREADS).
  -k TOP        Number of top mutations plotted and cross-checked (default: $TOP).

Long aliases: --help, --skip-download, --clean, --force.

Environment:
  PYTHON        Python interpreter (default: .venv/bin/python if present, else python3)
  NCBI_API_KEY  Optional NCBI API key (higher download rate limit)

Exit codes: 0 ok, 1 stage failure, 2 bad usage, 3 missing file,
            4 missing dependency, 5 checksum failure, 130/143 interrupted.

Outputs: results/, report/dashboard.html, logs/run_*.log, archives/hbv_results_*.tar.gz
EOF
}

# Translate the documented long aliases into short options for getopts.
args=()
for argument in "$@"; do
    case "$argument" in
        --help)          args+=("-h") ;;
        --skip-download) args+=("-s") ;;
        --clean)         args+=("-c") ;;
        --force)         args+=("-f") ;;
        --*)             printf 'Error: unknown option: %s\n\n' "$argument" >&2
                         usage >&2
                         exit "$EXIT_USAGE" ;;
        *)               args+=("$argument") ;;
    esac
done
set -- "${args[@]+"${args[@]}"}"

while getopts ":hscfnt:j:k:" opt; do
    case "$opt" in
        h) usage; exit 0 ;;
        s) SKIP_DOWNLOAD=1 ;;
        c) CLEAN_RESULTS=1 ;;
        f) FORCE_DOWNLOAD=1 ;;
        n) NO_ARCHIVE=1 ;;
        t) TAXON="$OPTARG" ;;
        j) THREADS="$OPTARG" ;;
        k) TOP="$OPTARG" ;;
        :) printf 'Error: option -%s requires a value\n\n' "$OPTARG" >&2; usage >&2; exit "$EXIT_USAGE" ;;
        \?) printf 'Error: unknown option: -%s\n\n' "$OPTARG" >&2; usage >&2; exit "$EXIT_USAGE" ;;
    esac
done
shift $((OPTIND - 1))
if (( $# > 0 )); then
    printf 'Error: unexpected argument: %s\n\n' "$1" >&2
    usage >&2
    exit "$EXIT_USAGE"
fi
[[ "$TAXON" =~ ^[0-9]+$ ]]       || die "-t expects a numeric NCBI taxonomy id, got '$TAXON'" "$EXIT_USAGE"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "-j expects a positive integer, got '$THREADS'" "$EXIT_USAGE"
[[ "$TOP" =~ ^[1-9][0-9]*$ ]]     || die "-k expects a positive integer, got '$TOP'" "$EXIT_USAGE"

# -------------------------------------------- temporary files and signals
mkdir -p "$LOG_DIR"
TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/hbv_pipeline.XXXXXX")"
# Every child (Python tempfile, sort -T, ...) writes its temporary files here.
export TMPDIR="$TMP_DIR"
export PYTHONUNBUFFERED=1

cleanup() {
    local status=$?
    trap - EXIT INT TERM ERR
    if [[ -n "${CHILD_PID:-}" ]] && kill -0 "$CHILD_PID" 2>/dev/null; then
        # Stop the running stage and its own children (e.g. MAFFT).
        pkill -TERM -P "$CHILD_PID" 2>/dev/null || true
        kill -TERM "$CHILD_PID" 2>/dev/null || true
        wait "$CHILD_PID" 2>/dev/null || true
    fi
    rm -rf "$TMP_DIR"
    # Half-written staging areas left by an interrupted stage.
    rm -rf "$TASK1_DIR.staging" "$TASK1_DIR.old" "$RAW_DIR/extracted.tmp" \
           "$BASE_DIR"/.hbv-preprocess-* "$ARCHIVE_DIR"/*.part
    rm -f "$PID_FILE"
    if (( status == 0 )); then
        log "pipeline finished successfully — log: ${LOG_FILE#"$BASE_DIR"/}"
    else
        printf '[%s] [ERROR] pipeline stopped (exit code %s); temporary files removed. Log: %s\n' \
            "$(_timestamp)" "$status" "${LOG_FILE#"$BASE_DIR"/}" >&2
    fi
    exit "$status"
}
on_signal() {
    printf '[%s] [ERROR] received SIG%s — stopping and cleaning up\n' "$(_timestamp)" "$1" >&2
    [[ "$1" == INT ]] && exit 130 || exit 143
}
on_error() {
    printf '[%s] [ERROR] command failed (exit %s) at line %s: %s\n' \
        "$(_timestamp)" "$1" "$2" "$3" >&2
}
trap cleanup EXIT
trap 'on_signal INT' INT
trap 'on_signal TERM' TERM
trap 'on_error $? $LINENO "$BASH_COMMAND"' ERR

# Refuse to run twice at the same time (e.g. cron + manual run); a stale PID
# file from a killed run is ignored, so no manual cleanup is ever needed.
if [[ -f "$PID_FILE" ]] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
    trap - EXIT
    rm -rf "$TMP_DIR"
    die "another run.sh is already running (PID $(cat "$PID_FILE"))" 1
fi
echo "$$" > "$PID_FILE"

# ------------------------------------------------------- timestamped log
LOG_FILE="$LOG_DIR/run_$(date +%Y%m%d_%H%M%S).log"
# Every line written to stdout/stderr is prefixed with a timestamp and copied
# into the log file; stdout and stderr stay separate streams on the terminal.
stamp_lines() {
    local line
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" == \[20[0-9][0-9]-* ]]; then
            printf '%s\n' "$line"            # already stamped by log()/warn()/die()
        else
            printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$line"
        fi
    done
}
exec > >(stamp_lines | tee -a "$LOG_FILE") \
    2> >(stamp_lines | tee -a "$LOG_FILE" >&2)

log "HBV Genome Pipeline — run started in $BASE_DIR"
log "options: skip_download=$SKIP_DOWNLOAD clean=$CLEAN_RESULTS force=$FORCE_DOWNLOAD taxon=$TAXON threads=$THREADS top=$TOP archive=$(( ! NO_ARCHIVE ))"

# ------------------------------------------------------------ preflight
if [[ -z "${PYTHON:-}" && -x "$BASE_DIR/.venv/bin/python" ]]; then
    PYTHON="$BASE_DIR/.venv/bin/python"
fi
PYTHON="${PYTHON:-python3}"

preflight() {
    local file
    for file in preprocess.py analyze_mutations.py statistics.py \
                report/make_dashboard.py scripts/shell_analysis.sh; do
        require_file "$BASE_DIR/$file" "the project checkout is incomplete"
    done
    if (( ! SKIP_DOWNLOAD )); then
        require_file "$DC_DIR/download.sh" "the project checkout is incomplete"
    fi
    for cmd in grep sed awk cut sort uniq tr join sha256sum md5sum tar gzip; do
        require_cmd "$cmd"
    done
    require_cmd "$PYTHON" "install Python >= 3.11 or set PYTHON=/path/to/python"
    "$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)' \
        || die "Python >= 3.11 is required (found $("$PYTHON" -V 2>&1))" "$EXIT_MISSING_DEP"
    require_python_module "$PYTHON" pandas "run: bash scripts/setup_env.sh"
    require_python_module "$PYTHON" matplotlib "run: bash scripts/setup_env.sh"
    require_cmd mafft "install MAFFT >= 7.487 (Ubuntu/Debian: sudo apt install mafft)"
    log "preflight OK — $("$PYTHON" -V 2>&1), MAFFT $(mafft --version 2>&1 | head -n 1)"
}
preflight

if (( CLEAN_RESULTS )); then
    log "cleaning previously generated results (-c)"
    rm -rf "$RESULTS_DIR/figures" "$RESULTS_DIR/shell" "$RESULTS_DIR/preprocessing"
    rm -f "$RESULTS_DIR"/*.csv "$RESULTS_DIR"/*.json "$RESULTS_DIR"/*.sha256 "$BASE_DIR/report/dashboard.html"
fi
mkdir -p "$BASE_DIR/data" "$RESULTS_DIR/figures" "$RESULTS_DIR/shell" "$ARCHIVE_DIR"

# ------------------------------------------------------------- 1 download
if (( SKIP_DOWNLOAD )); then
    log "==> SKIP   download (-s): reusing data-collection/data/"
else
    # Raw data are kept read-only; unlock them only for the refresh.
    set_writable "$RAW_DIR"
    download_args=(--taxon "$TAXON")
    (( FORCE_DOWNLOAD )) && download_args+=(--force)
    run_stage "1/10 download NCBI dataset (taxon $TAXON)" \
        bash "$DC_DIR/download.sh" "${download_args[@]}"
fi

# ------------------------------------------------------------ 2 integrity
verify_integrity() {
    require_file "$TASK1_DIR/sequences.fasta" "run without -s to download the dataset"
    require_file "$TASK1_DIR/metadata.csv"    "run without -s to download the dataset"
    require_file "$CHECKSUM_DIR/dataset.sha256" "run without -s to download the dataset"
    if [[ -f "$RAW_DIR/extracted/md5sum.txt" ]]; then
        verify_md5_file "$RAW_DIR/extracted/md5sum.txt"
    else
        warn "no NCBI md5sum.txt found (local-file source?) — NCBI MD5 check skipped"
    fi
    if [[ -f "$CHECKSUM_DIR/ncbi_extracted.sha256" && -d "$RAW_DIR/extracted" ]]; then
        verify_sha256_file "$CHECKSUM_DIR/ncbi_extracted.sha256" "$DC_DIR"
    fi
    verify_sha256_file "$CHECKSUM_DIR/dataset.sha256" "$DC_DIR"
    set_readonly "$RAW_DIR"
    log "raw data protected read-only: ${RAW_DIR#"$BASE_DIR"/}"
}
run_stage "2/10 verify dataset integrity" verify_integrity

# ----------------------------------------------------- 3 shell profiling
run_stage "3/10 shell profile of the sequences" \
    bash "$SCRIPTS_DIR/shell_analysis.sh" -m sequences \
        -i "$TASK1_DIR/sequences.fasta" -o "$RESULTS_DIR/shell"

# ----------------------------------------------------------- 4..6 Python
run_stage "4/10 preprocess (HBV filter, validation)" \
    "$PYTHON" "$BASE_DIR/preprocess.py" \
        --input-dir "$TASK1_DIR" --output-dir "$CLEAN_DIR" \
        --reports-dir "$RESULTS_DIR/preprocessing"

run_stage "5/10 MAFFT alignment and mutation detection" \
    "$PYTHON" "$BASE_DIR/analyze_mutations.py" \
        --input "$CLEAN_DIR/hbv_sequences.fasta" --metadata "$CLEAN_DIR/hbv_metadata.csv" \
        --output "$RESULTS_DIR/mutations.csv" \
        --skipped "$RESULTS_DIR/mutation_analysis_skipped.csv" \
        --report "$RESULTS_DIR/mutation_analysis.json" --threads "$THREADS"

run_stage "6/10 statistics and figures" \
    "$PYTHON" "$BASE_DIR/statistics.py" \
        --mutations "$RESULTS_DIR/mutations.csv" --metadata "$CLEAN_DIR/hbv_metadata.csv" \
        --skipped "$RESULTS_DIR/mutation_analysis_skipped.csv" \
        --outdir "$RESULTS_DIR" --top "$TOP"

# --------------------------------------------------- 7 shell cross-check
run_stage "7/10 shell cross-check of the top mutations" \
    bash "$SCRIPTS_DIR/shell_analysis.sh" -m mutations \
        -i "$RESULTS_DIR/mutations.csv" -s "$RESULTS_DIR/summary.csv" \
        -o "$RESULTS_DIR/shell" -n "$TOP"

# ------------------------------------------------------------ 8 dashboard
run_stage "8/10 HTML dashboard" "$PYTHON" "$BASE_DIR/report/make_dashboard.py"

# --------------------------------------------------- 9 results checksums
results_checksums() {
    local list="$RESULTS_DIR/checksums.sha256" files=()
    # Every result file except the checksum list itself, in a stable order.
    mapfile -t files < <(cd "$BASE_DIR" && find results report/dashboard.html -type f \
        ! -name checksums.sha256 | LC_ALL=C sort)
    write_sha256_manifest "$list" "$BASE_DIR" "${files[@]}"
    verify_sha256_file "$list" "$BASE_DIR"
}
run_stage "9/10 results checksums" results_checksums

# -------------------------------------------------------------- 10 archive
if (( NO_ARCHIVE )); then
    log "==> SKIP   archive (-n)"
else
    make_archive() {
        local archive
        # The running log is still being written, so it stays in logs/ only.
        archive="$(archive_results "$ARCHIVE_DIR" "$BASE_DIR" results report/dashboard.html)"
        log "archive: ${archive#"$BASE_DIR"/} ($(du -h "$archive" | cut -f1))"
    }
    run_stage "10/10 dated results archive" make_archive
fi

log "summary: $(tr -d '\n ' < "$RESULTS_DIR/statistics_summary.json")"
