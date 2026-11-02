#!/bin/bash
# HVP integration suite. The version under test is nextflow.config's
# manifest.version, printed in the summary below.
#
# Run from anywhere: bash tests/run_tests.sh
# Under SLURM:       sbatch tests/submit_tests.sh
# Requires: the hutlab module environment loaded, and Java 21 on PATH.
#
# What CI cannot check -- real tools, real databases, a cluster. The graph
# itself (every --components / --plugin selection wires up) is checked in CI by
# `nextflow run -preview`, and the registry by bin/make_diagrams.py --check.
#
#   test              components                         input
#   ----------------- ---------------------------------- -----------------------
#   baqlava_se        kneaddata,baqlava                  bundled single-end
#   baqlava_noqc_pe   baqlava (no QC: pairs concatenated) bundled paired-end
#   placeholders      kneaddata,host_scrubber,phanta     bundled single-end
#   assembly_pe       kneaddata,assembly                 simulated MAG fixture
#   plugin_genomad    --plugin genomad (adds assembly)   simulated MAG fixture
#   bad_component     --components nope                  (fails at launch)
#
# The runs are independent, so their nextflow drivers run in parallel; each
# gets its own work directory. Checks run after all of them have finished.

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_BASE="${REPO_ROOT}/tests/results"
WORK_BASE="${REPO_ROOT}/tests/work"
NF="${NF:-/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin/nextflow}"
PROFILE="-profile harvard_rc"
SE="${REPO_ROOT}/tests/data/single_end_rawfastq"
PE="${REPO_ROOT}/tests/data/rawfastq"
# The assembly fixture is simulated from two real genomes, because the bundled
# reads cannot produce a MAG; regenerate it with
# `python3 ~/biobakery_assembly_test/make_fixture.py --output input_pe`.
ASM_FIXTURE="${ASSEMBLY_FIXTURE:-${HOME}/biobakery_assembly_test}/input_pe"
VERSION="$(sed -n "s/^ *version *= *'\(.*\)'.*/\1/p" "${REPO_ROOT}/nextflow.config")"

mkdir -p "${RESULTS_BASE}" "${WORK_BASE}"

pass=0; fail=0

# Run one nextflow invocation. Named log, named work directory, and a
# .status file so the parallel runs can be judged after the fact.
nf_run() {
    local name="$1"; shift
    local out="${RESULTS_BASE}/${name}"
    local launch="${WORK_BASE}/${name}_launch"
    # the .status file goes too: a leftover from a previous run reads as this
    # run's result while this one is still going
    rm -rf "${out}" "${WORK_BASE}/${name}" "${launch}"
    rm -f  "${RESULTS_BASE}/${name}.status"
    mkdir -p "${launch}"
    # Each driver launches from its own directory: concurrent runs sharing one
    # launch directory share .nextflow/ -- the history file and the cache -- and
    # tread on each other.
    (
        cd "${launch}" || exit 1
        "${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
            --outdir "${out}" \
            -work-dir "${WORK_BASE}/${name}" \
            "$@"
    ) > "${RESULTS_BASE}/${name}.log" 2>&1
    echo "$?" > "${RESULTS_BASE}/${name}.status"
}

# Did the run exit 0?
check_run() {
    local name="$1"
    local status
    status="$(cat "${RESULTS_BASE}/${name}.status" 2>/dev/null || echo missing)"
    if [ "${status}" = "0" ]; then
        echo "[PASS] ${name}: nextflow completed"
        ((pass++))
    else
        echo "[FAIL] ${name}: nextflow exited ${status}"
        grep -E "^ERROR|Caused by:|Command error:" -A3 "${RESULTS_BASE}/${name}.log" | head -20 >&2
        ((fail++))
    fi
}

# Are the named output patterns all present under the run's output folder?
check_files() {
    local name="$1"; shift
    local dir="${RESULTS_BASE}/${name}"
    local ok=true
    for p in "$@"; do
        if ! find "${dir}" -name "${p}" 2>/dev/null | grep -q .; then
            echo "  MISSING in ${name}: ${p}" >&2
            ok=false
        fi
    done
    if $ok; then
        echo "[PASS] ${name}: expected outputs present"
        ((pass++))
    else
        echo "[FAIL] ${name}: expected outputs missing"
        ((fail++))
    fi
}

note() { echo ""; echo "=== $* ==="; }

# Does the run's log contain this text?
check_log() {
    local name="$1" text="$2" what="$3"
    if grep -qF -- "${text}" "${RESULTS_BASE}/${name}.log"; then
        echo "[PASS] ${name}: ${what}"
        ((pass++))
    else
        echo "[FAIL] ${name}: ${what} (no '${text}' in the log)"
        ((fail++))
    fi
}

# ── Launch ───────────────────────────────────────────────────────────────────
note "Launching"

# The tutorial sample has no species MetaPhlAn can detect, so BAQLaVa's
# bacterial depletion has nothing to work from: bypass it.
nf_run baqlava_se --readsdir "${SE}" --single_end true \
    --components kneaddata,baqlava --baqlava_bypass_depletion true &

nf_run baqlava_noqc_pe --readsdir "${PE}" \
    --components baqlava --baqlava_humann false --baqlava_bypass_depletion true &

nf_run placeholders --readsdir "${SE}" --single_end true \
    --components kneaddata,host_scrubber,phanta &

if [ -d "${ASM_FIXTURE}" ]; then
    nf_run assembly_pe    --readsdir "${ASM_FIXTURE}" --components kneaddata,assembly &
    nf_run plugin_genomad --readsdir "${ASM_FIXTURE}" --plugin genomad &
fi

nf_run bad_component --readsdir "${SE}" --components kneaddata,nope --log_versions false &

wait

# ── Checks ───────────────────────────────────────────────────────────────────
note "BAQLaVa, single-end, with KneadData"
check_run   baqlava_se
check_files baqlava_se "*_BAQLaVa_profile.txt" "merged_metaphlan_profiles.tsv" \
                       "*_genefamilies*.tsv" "kneaddata_read_count_table.tsv" "pipeline_log.txt"
check_log   baqlava_se "[hvp] Components: kneaddata -> baqlava" "components resolved"

note "BAQLaVa, paired-end, no QC"
check_run   baqlava_noqc_pe
check_files baqlava_noqc_pe "*_merged.fastq*" "*_BAQLaVa_profile.txt"

note "Placeholders run, warn, and produce nothing"
check_run   placeholders
check_log   placeholders "host_scrubber: skeleton component" "host_scrubber warns"
check_log   placeholders "phanta: skeleton component" "phanta warns"

if [ -d "${ASM_FIXTURE}" ]; then
    note "Assembly component"
    check_run   assembly_pe
    check_files assembly_pe "*.final.contigs.fa" "*.contig_depths.txt" "SGB_info.tsv" "final_profile.tsv"

    note "Plugin mode resolves dependencies"
    check_run   plugin_genomad
    check_log   plugin_genomad "adding 'assembly'" "assembly added for genomad"
    check_files plugin_genomad "*.final.contigs.fa"
else
    note "Assembly and plugin tests SKIPPED: no fixture at ${ASM_FIXTURE}"
fi

note "Unknown component fails at launch"
if [ "$(cat "${RESULTS_BASE}/bad_component.status")" != "0" ]; then
    echo "[PASS] bad_component: run refused"
    ((pass++))
else
    echo "[FAIL] bad_component: run was accepted"
    ((fail++))
fi
check_log bad_component "unknown component(s): nope" "error names the component"

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  biobakery/hvp ${VERSION}"
echo "  Results: ${pass} passed, ${fail} failed"
echo "  Logs:    ${RESULTS_BASE}/"
echo "========================================"
[ "${fail}" -eq 0 ]
