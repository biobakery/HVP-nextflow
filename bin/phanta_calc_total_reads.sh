#!/usr/bin/env bash
#
# Per-sample read accounting for total_reads.tsv.
#
# Adapted from Phanta's pipeline_scripts/calc_total_reads.sh at commit ecbfcf9.
# The only change is that the two input files are named explicitly instead of
# being globbed out of a shared classification directory, which is what makes
# the task isolatable. Deliberately still shell: the arithmetic is integer
# arithmetic in bash and bc, and reimplementing it in Python risks a different
# string formatting against a byte-identity bar.
#
# Emits one tab-separated row:
#   sample  Tot_Samp_Reads  Unassigned_Step_One  Assigned_Step_Three  Unassigned_Step_Three
#
# Exits 64 with "No reads in sample" when the Kraken report has neither a root
# nor an unclassified line. That exit code is part of the contract.
set -uo pipefail

usage() {
    echo "usage: $0 --sample NAME --krak-report FILE --bracken-table FILE --out FILE" >&2
    exit 2
}

SAMPLE=""; KRAK=""; BRACKEN=""; OUT=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)        SAMPLE="$2"; shift 2 ;;
        --krak-report)   KRAK="$2";   shift 2 ;;
        --bracken-table) BRACKEN="$2"; shift 2 ;;
        --out)           OUT="$2";    shift 2 ;;
        -h|--help)       usage ;;
        *) echo "unknown argument: $1" >&2; usage ;;
    esac
done
[[ -n "$SAMPLE" && -n "$KRAK" && -n "$BRACKEN" && -n "$OUT" ]] || usage

if ! command -v bc >/dev/null 2>&1; then
    echo "Error: 'bc' calculator is not installed." >&2
    exit 1
fi

# Total reads in the sample: everything Kraken2 put under root, plus everything
# it could not classify at all.
root=$(grep 'root' "$KRAK" | head -n 1 | cut -f 2 | xargs)
unclassified_krak=$(grep 'unclassified' "$KRAK" | head -n 1 | cut -f 2 | xargs)

if [[ -z $root ]] && [[ -z $unclassified_krak ]]; then
    echo "No reads in sample"
    exit 64
elif [[ -z $root ]]; then
    tot_reads=$unclassified_krak
elif [[ -z $unclassified_krak ]]; then
    tot_reads=$root
else
    tot_reads=$((root + unclassified_krak))
fi

# Bracken-classified reads: column 6 of the Bracken table, summed. An ineligible
# sample has an empty table here, so this comes back empty and the arithmetic
# below treats it as zero -- which is the behaviour prep_to_merge_normed relies
# on when it maps an empty count to 0.
tot_reads_bracken=$(cut -f 6 "$BRACKEN" | tail -n +2 | paste -sd+ | bc)

unclassified_brack="$((tot_reads - tot_reads_bracken))"

echo -e "$SAMPLE\t${tot_reads}\t${unclassified_krak}\t${tot_reads_bracken}\t${unclassified_brack}" > "$OUT"
