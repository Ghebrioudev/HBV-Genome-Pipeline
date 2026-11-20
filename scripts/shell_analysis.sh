#!/usr/bin/env bash
#
# scripts/shell_analysis.sh — data profiling and cross-checks done with the
# standard text tools (grep, sed, awk, cut, sort, uniq, tr, join),
# chained with pipes. Called automatically by run.sh; can also be run alone.
#
# Usage:
#   scripts/shell_analysis.sh -m sequences -i FASTA [-o OUTDIR]
#   scripts/shell_analysis.sh -m mutations -i MUTATIONS_CSV -s SUMMARY_CSV [-o OUTDIR] [-n TOP]
#   scripts/shell_analysis.sh -h
#
# Mode "sequences" (on the Task 1 FASTA, before Python preprocessing):
#   sequence_profile.tsv          record count, duplicated accessions, total nt
#   length_distribution.tsv       number of records per 500-nt length bin
#   header_annotations.tsv        "complete genome" / "partial ..." labels in headers
#   base_composition.tsv          count of every nucleotide symbol (A,C,G,T,IUPAC...)
#
# Mode "mutations" (on results/mutations.csv, after Task 3/4):
#   events_by_type.tsv            mutation events per type and call status
#   top_mutations_shell.tsv       most frequent definitive mutations, counted in shell
#   crosscheck_top_mutations.tsv  shell counts compared with Python's summary.csv
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
# shellcheck source=lib.sh
source "$SCRIPT_DIR/lib.sh"
# Any failing command is reported with its line number before exiting.
trap 'die "shell_analysis.sh: command failed at line $LINENO: $BASH_COMMAND"' ERR

MODE="" INPUT="" SUMMARY="" OUTDIR="$BASE_DIR/results/shell" TOP=20

usage() { sed -n '3,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while getopts ":m:i:s:o:n:h" opt; do
    case "$opt" in
        m) MODE="$OPTARG" ;;
        i) INPUT="$OPTARG" ;;
        s) SUMMARY="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        n) TOP="$OPTARG" ;;
        h) usage; exit 0 ;;
        :) die "option -$OPTARG requires a value (see -h)" "$EXIT_USAGE" ;;
        \?) die "unknown option: -$OPTARG (see -h)" "$EXIT_USAGE" ;;
    esac
done

[[ "$TOP" =~ ^[0-9]+$ ]] || die "-n expects a positive integer, got '$TOP'" "$EXIT_USAGE"
[[ -n "$INPUT" ]] || die "missing -i INPUT (see -h)" "$EXIT_USAGE"
require_file "$INPUT"
mkdir -p "$OUTDIR"

# Sorting is locale-independent and its temporary files go to $TMPDIR, which
# run.sh points at its own temporary directory (removed by its trap).
export LC_ALL=C
SORT_TMP="${TMPDIR:-/tmp}"

profile_sequences() {
    local fasta="$1" records duplicates total_nt

    # Pipe 1: headers -> accession -> sorted -> duplicated ones -> count
    records="$(grep -c '^>' "$fasta" || true)"
    duplicates="$(grep '^>' "$fasta" | cut -d' ' -f1 | sed 's/^>//' \
        | sort -T "$SORT_TMP" | uniq -d | wc -l | tr -d ' ')"
    total_nt="$(grep -v '^>' "$fasta" | tr -d '\n\r' | wc -c | tr -d ' ')"
    {
        printf 'metric\tvalue\n'
        printf 'fasta_records\t%s\n' "$records"
        printf 'duplicated_accessions\t%s\n' "$duplicates"
        printf 'total_nucleotides\t%s\n' "$total_nt"
    } > "$OUTDIR/sequence_profile.tsv"

    # Pipe 2: per-record length (awk) -> 500 nt bin (awk) -> sort -> uniq -c -> TSV
    {
        printf 'length_bin_start\trecords\n'
        awk '/^>/ { if (seen) print len; len = 0; seen = 1; next }
             { len += length($0) }
             END { if (seen) print len }' "$fasta" \
            | awk '{ print int($1 / 500) * 500 }' \
            | sort -T "$SORT_TMP" -n | uniq -c \
            | awk '{ printf "%s\t%s\n", $2, $1 }'
    } > "$OUTDIR/length_distribution.tsv"

    # Pipe 3: headers -> description -> completeness label -> lower case -> count
    {
        printf 'records\theader_label\n'
        grep '^>' "$fasta" | cut -d' ' -f2- \
            | { grep -oiE '(complete|partial) (genome|cds|sequence)' || true; } \
            | tr '[:upper:]' '[:lower:]' \
            | sort -T "$SORT_TMP" | uniq -c | sort -k1,1nr \
            | sed -E 's/^ *([0-9]+) /\1\t/'
    } > "$OUTDIR/header_annotations.tsv"

    # Pipe 4: sequence lines -> upper case -> per-symbol counts (awk gsub) -> sorted
    {
        printf 'symbol\tcount\n'
        grep -v '^>' "$fasta" | tr -d '\r' | tr '[:lower:]' '[:upper:]' \
            | awk 'BEGIN { n = split("A C G T R Y S W K M B D H V N U -", sym, " ") }
                   { line = $0
                     for (i = 1; i <= n; i++) count[sym[i]] += gsub(sym[i], "", line)
                     other += length(line) }
                   END { for (i = 1; i <= n; i++) if (count[sym[i]] > 0) print sym[i] "\t" count[sym[i]]
                         if (other > 0) print "other\t" other }' \
            | sort -t "$(printf '\t')" -k2,2nr
    } > "$OUTDIR/base_composition.tsv"

    log "  sequences: $records records, $duplicates duplicated accessions, $total_nt nt"
}

profile_mutations() {
    local mutations="$1" summary="$2" header
    require_file "$summary"

    header="$(head -n 1 "$mutations")"
    [[ "$header" == sample_id,reference_accession,position,reference_base,observed_base,variation_type,call_status* ]] \
        || die "unexpected header in $mutations: $header" 1

    # Pipe 5: events per (type, call_status)
    {
        printf 'variation_type\tcall_status\tevents\n'
        tail -n +2 "$mutations" | cut -d, -f6,7 \
            | sort -T "$SORT_TMP" | uniq -c | sort -k1,1nr \
            | awk '{ split($2, f, ","); printf "%s\t%s\t%s\n", f[1], f[2], $1 }'
    } > "$OUTDIR/events_by_type.tsv"

    # Pipe 6: definitive calls -> "pos:REF>ALT" key -> count -> top N
    {
        printf 'mutation\tevents\n'
        tail -n +2 "$mutations" | cut -d, -f3-5,7 \
            | { grep ',definitive$' || true; } \
            | awk -F, '{ print $1 ":" $2 ">" $3 }' \
            | sort -T "$SORT_TMP" | uniq -c | sort -k1,1nr -k2,2 \
            | awk -v top="$TOP" 'NR <= top { printf "%s\t%s\n", $2, $1 }'
    # (no "head" here: it would close the pipe early, and with pipefail the
    #  resulting SIGPIPE of the upstream sort would abort the script.)
    } > "$OUTDIR/top_mutations_shell.tsv"

    # Pipe 7: compare with Python (summary.csv columns: mutation,...,n_sequences,...)
    {
        printf 'mutation\tshell_events\tpython_n_sequences\tmatch\n'
        join -t "$(printf '\t')" \
            <(tail -n +2 "$OUTDIR/top_mutations_shell.tsv" | sort -T "$SORT_TMP" -t "$(printf '\t')" -k1,1) \
            <(tail -n +2 "$summary" | cut -d, -f1,6 | tr ',' '\t' | sort -T "$SORT_TMP" -t "$(printf '\t')" -k1,1) \
            | awk -F'\t' '{ printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($2 == $3 ? "yes" : "NO") }' \
            | sort -t "$(printf '\t')" -k2,2nr
    } > "$OUTDIR/crosscheck_top_mutations.tsv"

    local mismatches
    mismatches="$(tail -n +2 "$OUTDIR/crosscheck_top_mutations.tsv" | cut -f4 | { grep -c '^NO$' || true; })"
    if (( mismatches > 0 )); then
        warn "  $mismatches of the top $TOP mutations differ between shell and Python counts"
    else
        log "  top $TOP mutations: shell counts match statistics.py"
    fi
}

case "$MODE" in
    sequences) log "shell profile of $INPUT"; profile_sequences "$INPUT" ;;
    mutations) [[ -n "$SUMMARY" ]] || die "mode 'mutations' needs -s SUMMARY_CSV" "$EXIT_USAGE"
               log "shell analysis of $INPUT"; profile_mutations "$INPUT" "$SUMMARY" ;;
    "")        die "missing -m MODE (sequences|mutations), see -h" "$EXIT_USAGE" ;;
    *)         die "unknown mode '$MODE' (expected sequences|mutations)" "$EXIT_USAGE" ;;
esac
log "shell outputs written to ${OUTDIR#"$BASE_DIR"/}/"
