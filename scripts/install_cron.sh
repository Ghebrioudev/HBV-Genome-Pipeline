#!/usr/bin/env bash
#
# scripts/install_cron.sh — install, show or remove the weekly cron job of
# the HBV Genome Pipeline (template: scripts/cron/hbv-pipeline.cron).
#
# Usage:
#   scripts/install_cron.sh          install (or update) the job — idempotent
#   scripts/install_cron.sh -p       print the line that would be installed
#   scripts/install_cron.sh -r       remove the job
#   scripts/install_cron.sh -h       help
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"

TEMPLATE="$SCRIPT_DIR/cron/hbv-pipeline.cron"
MARKER="# hbv-genome-pipeline"
ACTION="install"

while getopts ":prh" opt; do
    case "$opt" in
        p) ACTION="print" ;;
        r) ACTION="remove" ;;
        h) sed -n '3,11p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        \?) die "unknown option: -$OPTARG (see -h)" "$EXIT_USAGE" ;;
    esac
done

require_file "$TEMPLATE"
# Pipe: template -> drop comments/blank lines -> substitute the project path
JOB="$(grep -v '^#' "$TEMPLATE" | grep -v '^[[:space:]]*$' | sed "s|@PROJECT_DIR@|$BASE_DIR|g")"

if [[ "$ACTION" == print ]]; then
    printf '%s\n' "$JOB"
    exit 0
fi

require_cmd crontab "install cron (Ubuntu/Debian: sudo apt install cron)"
CURRENT="$(crontab -l 2>/dev/null | grep -vF "$MARKER" || true)"

if [[ "$ACTION" == remove ]]; then
    printf '%s\n' "$CURRENT" | grep -v '^$' | crontab - || crontab -r
    log "cron job removed"
    exit 0
fi

mkdir -p "$BASE_DIR/logs"
{ [[ -n "$CURRENT" ]] && printf '%s\n' "$CURRENT"; printf '%s\n' "$JOB"; } | crontab -
log "cron job installed:"
crontab -l | grep -F "$MARKER" >&2
