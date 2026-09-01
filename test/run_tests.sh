#!/bin/bash
# biobakery-workflows-nextflow v0.0.2 — integration test suite
# Run from the repo root: bash test/run_tests.sh
# Requires: hutlab module loaded, Java 21 in PATH

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RESULTS_BASE="${REPO_ROOT}/test/results"
NF="/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin/nextflow"
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
echo "=== Test 4: toggle guard — run_viral_profiling=true without taxonomic_profiling ==="
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx \
    --readsdir "${REPO_ROOT}/test/single_end_rawfastq" \
    --filepattern "*.fastq.gz" \
    --paired_end false \
    --outdir /tmp/nf_guard_test \
    --run_qc false \
    --run_taxonomic_profiling false \
    --run_functional_profiling false \
    --run_viral_profiling true \
    --log_versions false \
    2>&1 | grep "ERROR" | head -2 | tee "${RESULTS_BASE}/test4.log" || true

if grep -qE "run_taxonomic_profiling = (false, but these require it|true)" "${RESULTS_BASE}/test4.log"; then
    echo "[PASS] viral-without-taxonomic guard fires correctly"
    ((pass++)) || true
else
    echo "[FAIL] viral-without-taxonomic guard did NOT fire"
    ((fail++)) || true
fi


# ── Test 5: MTX — KneadData (3 reference DBs) + MetaPhlAn ─────────────────────
echo ""
echo "=== Test 5: MTX paired-end — KneadData x3 DBs + MetaPhlAn (SLURM) ==="
OUT="${RESULTS_BASE}/test5_mtx"
rm -rf "${OUT}" "${REPO_ROOT}/work_test5"
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mtx \
    --readsdir "${REPO_ROOT}/test/rawfastq" \
    --outdir "${OUT}" \
    --run_functional_profiling false \
    ${COMMON} \
    --log_versions false \
    -work-dir "${REPO_ROOT}/work_test5" \
    2>&1 | tee "${RESULTS_BASE}/test5.log" | grep -E "Launching|Completed|ERROR|WARN" || true

check "MTX produces cleaned reads and MetaPhlAn profiles" "${OUT}" \
    "*_kneaddata.log" "*_profile_*.tsv" "merged_metaphlan_profiles.tsv"

# The metatranscriptome database set is the point of this workflow: KneadData
# should have written contaminant files for the host genome, the host mRNA and
# the rRNA index, not just the host genome.
mtx_dbs=$(find "${REPO_ROOT}/work_test5" -name "*_contam*.fastq*" -printf "%f\n" 2>/dev/null \
          | sed -E 's/.*_kneaddata_(.*)_bowtie2_.*/\1/' | sort -u | wc -l)
if [ "${mtx_dbs}" -ge 3 ]; then
    echo "[PASS] MTX KneadData used ${mtx_dbs} reference databases"
    ((pass++)) || true
else
    echo "[FAIL] MTX KneadData used ${mtx_dbs} reference database(s), expected 3"
    ((fail++)) || true
fi

# ── Test 6: mgx_mtx input guard — both folders are required ───────────────────
echo ""
echo "=== Test 6: mgx_mtx guard — missing --input_metatranscriptome must error ==="
"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx_mtx \
    --input_metagenome "${REPO_ROOT}/test/rawfastq" \
    --outdir "${RESULTS_BASE}/test6_guard" \
    --log_versions false \
    2>&1 | grep "ERROR" | head -3 | tee "${RESULTS_BASE}/test6.log" || true

if grep -q "input_metatranscriptome" "${RESULTS_BASE}/test6.log"; then
    echo "[PASS] mgx_mtx two-input guard fires correctly"
    ((pass++)) || true
else
    echo "[FAIL] mgx_mtx two-input guard did NOT fire"
    ((fail++)) || true
fi

# ── Test 7: mgx_mtx — mapped profiles, separate output folders ────────────────
echo ""
echo "=== Test 7: mgx_mtx paired halves with a mapping file (SLURM) ==="
OUT="${RESULTS_BASE}/test7_mgx_mtx"
FIX="${RESULTS_BASE}/test7_fixture"
rm -rf "${OUT}" "${FIX}" "${REPO_ROOT}/work_test7"
mkdir -p "${FIX}/mtx_in"
# Stand-in RNA samples: the same reads under distinct names, which is enough to
# exercise the mapping-driven pairing even though the ratios are meaningless.
for r in 1 2; do
    cp "${REPO_ROOT}/test/rawfastq/FG00004_S26_R${r}.fastq.gz" "${FIX}/mtx_in/RNA_FG00004_S26_R${r}.fastq.gz"
done
printf '#rna\tdna\nRNA_FG00004_S26\tFG00004_S26\n' > "${FIX}/mapping.tsv"

"${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
    --workflow mgx_mtx \
    --input_metagenome "${REPO_ROOT}/test/rawfastq" \
    --input_metatranscriptome "${FIX}/mtx_in" \
    --input_mapping "${FIX}/mapping.tsv" \
    --outdir "${OUT}" \
    --run_functional_profiling false \
    ${COMMON} \
    --log_versions false \
    -work-dir "${REPO_ROOT}/work_test7" \
    2>&1 | tee "${RESULTS_BASE}/test7.log" | grep -E "Launching|Completed|ERROR|WARN" || true

check "mgx_mtx publishes the two halves separately" "${OUT}" \
    "whole_metagenome_shotgun" "whole_metatranscriptome_shotgun"

# With a mapping file the RNA samples borrow the DNA sample's taxonomic profile,
# so the metatranscriptome half must not have run MetaPhlAn of its own.
if grep -q "TAX_MTX:metaphlan" "${RESULTS_BASE}/test7.log"; then
    echo "[FAIL] TAX_MTX:metaphlan ran despite a mapping file being given"
    ((fail++)) || true
else
    echo "[PASS] mapped run reuses the metagenome taxonomic profiles"
    ((pass++)) || true
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
echo "  Results: ${pass} passed, ${fail} failed"
echo "  Logs:    ${RESULTS_BASE}/"
echo "========================================"
[ "${fail}" -eq 0 ]
