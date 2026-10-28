#!/usr/bin/env bash
#
# Compare two Phanta output directories, and exit nonzero on any difference.
#
# A thin wrapper around `diff -r`, deliberately. An earlier plan called for a
# bespoke comparison tool with cell-level diff artifacts and bounded previews;
# for four TSVs and a handful of text files, `diff -r` reports more usefully and
# is a tool everyone already knows.
#
# Used for three different jobs (see the migration spec, "Comparison model"):
#   1. the routine fixture regression, against a committed baseline
#   2. validating a new dataset, upstream vs this pipeline
#   3. checking a baseline refresh reproduces itself
#
# Two things are always excluded:
#   *.krak            the per-read assignment file. It is the only artifact
#                     whose content depends on Kraken2's thread scheduling, so
#                     it differs run to run even on identical input. Not
#                     published either.
#   pipeline_info/    provenance and Nextflow's own reports. The database
#                     manifest records a path and a timestamp, and the
#                     timeline/report/trace/dag files carry the launch time in
#                     their names, so none of it is byte-stable by design.
#
# There are no numeric tolerances. Byte-identity is the bar; investigate every
# mismatch rather than widening a threshold.
set -uo pipefail

usage() {
    echo "usage: $0 [--include-krak] <expected_dir> <actual_dir>" >&2
    echo >&2
    echo "  Compares the two trees recursively, reports every missing," >&2
    echo "  unexpected and differing file, and exits 1 if there are any." >&2
    exit 2
}

INCLUDE_KRAK=0
ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --include-krak) INCLUDE_KRAK=1; shift ;;
        -h|--help) usage ;;
        -*) echo "unknown option: $1" >&2; usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[[ ${#ARGS[@]} -eq 2 ]] || usage

EXPECTED="${ARGS[0]}"
ACTUAL="${ARGS[1]}"

for d in "$EXPECTED" "$ACTUAL"; do
    [[ -d "$d" ]] || { echo "ERROR: not a directory: $d" >&2; exit 2; }
done

EXCLUDES=(--exclude='pipeline_info')
[[ "$INCLUDE_KRAK" == 0 ]] && EXCLUDES+=(--exclude='*.krak')

echo "comparing:"
echo "  expected  $EXPECTED"
echo "  actual    $ACTUAL"
[[ "$INCLUDE_KRAK" == 0 ]] && echo "  excluding *.krak and pipeline_info/"
echo

# -r walks the whole tree rather than stopping at the first mismatch, which is
# the point: one run should report every difference, not the first one.
output=$(diff -r "${EXCLUDES[@]}" "$EXPECTED" "$ACTUAL" 2>&1)
status=$?

if [[ $status -eq 0 ]]; then
    count=$(find "$EXPECTED" -type f ! -name '*.krak' ! -path '*/pipeline_info/*' | wc -l | tr -d ' ')
    echo "IDENTICAL: $count files match."
    exit 0
fi

echo "$output"
echo
echo "DIFFERENCES FOUND. Byte-identity is the acceptance bar, so each of these"
echo "needs explaining rather than tolerating. Common causes:"
echo "  * a different database   -- compare pipeline_info/phanta_db_manifest.txt"
echo "  * a different container  -- pandas float formatting is version-sensitive"
echo "  * a reordered samplesheet -- merged-table columns follow it, and"
echo "                               total_reads.tsv rows are sorted by name"
echo
echo "Comparing against an UPSTREAM Phanta output directory? These differences"
echo "are structural and expected -- they are publishing-contract choices, not"
echo "content. Verified byte-identical where both sides have the file:"
echo "  * <samp>.krak.report.filtered.bracken.temp        upstream only"
echo "  * <samp>.krak.report_bracken_species.filtered.temp  upstream only"
echo "  * <samp>.krak.report.filtered.bracken.scaled.temp   upstream only"
echo "      -- upstream publishes its three zero-byte canary files; this"
echo "         pipeline replaces them with an explicit status record."
echo "  * classification/intermediate/                     upstream only"
echo "  * <samp>.krak.report.filtering_decisions.txt        ours only"
echo "      -- two halves of one difference. Upstream MOVES the filtering"
echo "         decisions and the raw reports into intermediate/ after the"
echo "         fact; this pipeline publishes the decisions at classification/"
echo "         for every sample and keeps the rest in the work directory."
echo "  * pipeline_completed.txt                           upstream only"
echo "      -- open question 5 in the migration spec; not emitted."
echo
echo "To compare only what both sides publish in the same place:"
echo "  diff -r -x '*.krak' -x pipeline_info -x '*.temp' \\"
echo "       -x pipeline_completed.txt -x intermediate <expected> <actual>"
exit 1
