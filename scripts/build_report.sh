#!/usr/bin/env bash
#
# scripts/build_report.sh — compile report/rapport.md into rapport.pdf
# (4 pages maximum) with pandoc + a LaTeX engine. Not part of run.sh: the
# report is a written deliverable, built once the results exist.
#
# Usage: bash scripts/build_report.sh [-e ENGINE] [-f FONT] [-h]
#   -e ENGINE   LaTeX engine (default: xelatex; needs a TeX distribution)
#   -f FONT     main font with the symbols >= and ~= (default: DejaVu Serif;
#               e.g. -f Cambria on Windows)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"

ENGINE="xelatex"
FONT="DejaVu Serif"
while getopts ":e:f:h" opt; do
    case "$opt" in
        e) ENGINE="$OPTARG" ;;
        f) FONT="$OPTARG" ;;
        h) sed -n '3,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
        :) die "option -$OPTARG requires a value" "$EXIT_USAGE" ;;
        \?) die "unknown option: -$OPTARG (see -h)" "$EXIT_USAGE" ;;
    esac
done

require_cmd pandoc "https://pandoc.org/installing.html"
require_cmd "$ENGINE" "install a TeX distribution (TeX Live / MiKTeX)"
require_file "$BASE_DIR/report/rapport.md"
for figure in top_mutations.png mutations_along_genome.png; do
    require_file "$BASE_DIR/results/figures/$figure" "run ./run.sh first"
done

cd "$BASE_DIR/report"
pandoc rapport.md --from markdown --pdf-engine="$ENGINE" \
    --resource-path=.:.. -V mainfont="$FONT" -o "$BASE_DIR/rapport.pdf"
log "rapport.pdf written ($(du -h "$BASE_DIR/rapport.pdf" | cut -f1))"
