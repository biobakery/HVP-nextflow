#!/bin/bash
# biobakery-workflows-nextflow v0.0.2 — integration test suite
# Run from the repo root: bash test/run_tests.sh
# Requires: hutlab module loaded, Java 21 in PATH

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_BASE="${REPO_ROOT}/test/results"
NF="/n/huttenhower_lab/tools/nextflow/24.10.4/bin/nextflow"
PROFILE="-profile harvard_rc -params-file ${REPO_ROOT}/conf/harvard_rc.yaml"
COMMON="--run_viral_profiling false --run_strain_profiling false"

mkdir -p "${RESULTS_BASE}"

pass=0; fail=0

check() {
    local label="$1"; local result_dir="$2"; shift 2
    local patterns=("$@")
    local ok=true
    for p in "${patterns[@]}"; do
        if ! find "${result_dir}" -name "${p}" | grep -q .; then
            echo "  MISSING: ${p}" >&2
            ok=false
        fi
    done
    if $ok; then
        echo "[PASS] ${label}"
        ((pass++)) || true
    else
        echo "[FAIL] ${label}"
        ((fail++)) || true
    fi
}

# ── Test 1: version_log only ──────────────────────────────────────────────────
echo ""
echo "=== Test 1: version_log (local executor, no SLURM jobs) ==="
OUT="${RESULTS_BASE}/test1_versionlog"
rm -rf "${OUT}" "${REPO_ROOT}/work_test1" "${REPO_ROOT}/.nextflow_test1"
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx \
    --readsdir "${REPO_ROOT}/test/single_end_rawfastq" \
    --filepattern "*.fastq.gz" \
    --paired_end false \
    --outdir "${OUT}" \
    --run_qc false \
    --run_taxonomic_profiling false \
    --run_functional_profiling false \
    ${COMMON} \
    --log_versions true \
    -work-dir "${REPO_ROOT}/work_test1" \
    2>&1 | tee "${RESULTS_BASE}/test1.log" | grep -E "Launching|Completed|ERROR|WARN" || true

check "version_log produces pipeline_log.txt" "${OUT}/pipeline_info" "pipeline_log.txt"

if [ -f "${OUT}/pipeline_info/pipeline_log.txt" ]; then
    echo "  --- pipeline_log.txt content ---"
    grep -E "Tool Versions|kneaddata:|metaphlan:|humann:|baqlava:|workflow:|humann_version:|metaphlan_version:" \
        "${OUT}/pipeline_info/pipeline_log.txt" | head -20
    echo "  ---------------------------------"
fi

# ── Test 2: MGX single-end, MetaPhlAn only (no QC, no HUMAnN) ────────────────
echo ""
echo "=== Test 2: MGX single-end — MetaPhlAn only (SLURM) ==="
OUT="${RESULTS_BASE}/test2_se_metaphlan"
rm -rf "${OUT}" "${REPO_ROOT}/work_test2"
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx \
    --readsdir "${REPO_ROOT}/test/single_end_rawfastq" \
    --filepattern "*.fastq.gz" \
    --paired_end false \
    --outdir "${OUT}" \
    --run_qc false \
    --run_functional_profiling false \
    ${COMMON} \
    --log_versions true \
    -work-dir "${REPO_ROOT}/work_test2" \
    2>&1 | tee "${RESULTS_BASE}/test2.log" | grep -E "Launching|Completed|ERROR|WARN|process" || true

check "single-end MetaPhlAn profiles" "${OUT}/metaphlan" "*_profile_*.tsv"
check "single-end merged MetaPhlAn table" "${OUT}/metaphlan" "merged_metaphlan_profiles.tsv"
check "single-end pipeline_log" "${OUT}/pipeline_info" "pipeline_log.txt"

# ── Test 3: MGX paired-end, QC + MetaPhlAn (no HUMAnN) ───────────────────────
echo ""
echo "=== Test 3: MGX paired-end — KneadData + MetaPhlAn (SLURM) ==="
OUT="${RESULTS_BASE}/test3_pe_qc_metaphlan"
rm -rf "${OUT}" "${REPO_ROOT}/work_test3"
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx \
    --readsdir "${REPO_ROOT}/test/rawfastq" \
    --paired_end true \
    --outdir "${OUT}" \
    --run_qc true \
    --run_functional_profiling false \
    ${COMMON} \
    --log_versions true \
    -work-dir "${REPO_ROOT}/work_test3" \
    2>&1 | tee "${RESULTS_BASE}/test3.log" | grep -E "Launching|Completed|ERROR|WARN|process" || true

check "paired-end kneaddata output" "${OUT}/kneaddata" "*_kneaddata.log"
check "paired-end MetaPhlAn profiles" "${OUT}/metaphlan" "*_profile_*.tsv"
check "paired-end merged MetaPhlAn table" "${OUT}/metaphlan" "merged_metaphlan_profiles.tsv"

# ── Test 4: toggle guard — viral without taxonomic must error ─────────────────
echo ""
echo "=== Test 4: toggle guard — run_baqlava=true without taxonomic_profiling ==="
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx \
    --readsdir "${REPO_ROOT}/test/single_end_rawfastq" \
    --filepattern "*.fastq.gz" \
    --paired_end false \
    --outdir /tmp/nf_guard_test \
    --run_qc false \
    --run_taxonomic_profiling false \
    --run_functional_profiling false \
    --run_baqlava true \
    --log_versions false \
    2>&1 | grep "ERROR" | head -2 | tee "${RESULTS_BASE}/test4.log" || true

if grep -q "requires run_taxonomic_profiling" "${RESULTS_BASE}/test4.log"; then
    echo "[PASS] viral-without-taxonomic guard fires correctly"
    ((pass++)) || true
else
    echo "[FAIL] viral-without-taxonomic guard did NOT fire"
    ((fail++)) || true
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  Results: ${pass} passed, ${fail} failed"
echo "  Logs:    ${RESULTS_BASE}/"
echo "========================================"
[ "${fail}" -eq 0 ]
