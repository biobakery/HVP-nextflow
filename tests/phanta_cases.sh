# Phanta integration cases, sourced by tests/run_tests.sh.
#
# Not a standalone script: it runs in run_tests.sh's shell and uses its
# nf_run / check_run / check_files helpers, its ${RESULTS_BASE} and
# ${WORK_BASE}, and its pass/fail counters, so there is still one suite and one
# summary. Sourced rather than inlined only so that adding Phanta did not mean
# re-indenting four hundred lines of working bioBakery cases.
#
# Run them:
#   PHANTA_ONLY=1 NF=nextflow PROFILE="-profile local_docker" bash tests/run_tests.sh
#   PHANTA=1      bash tests/run_tests.sh      # after the bioBakery cases
#
# Prerequisites are deliberately small: a container runtime and the committed
# fixtures. No cluster, no 25 GB database, no module environment. On a laptop
# that means `docker build --platform linux/amd64 -t phanta:0.0.5-dev
# containers/phanta` and nothing else.

PHANTA_FIXTURES="${REPO_ROOT}/tests/fixtures/phanta"
PHANTA_DB="${PHANTA_DB:-${PHANTA_FIXTURES}/tiny_db/db}"
PHANTA_BASELINE="${PHANTA_BASELINE:-${PHANTA_FIXTURES}/baseline}"
# The tiny database ships database100mers.kmer_distrib and no other, so the
# read length is not a free choice here. Case 25 asserts that mismatching it
# fails with a message that says so.
PHANTA_READ_LENGTH="${PHANTA_READ_LENGTH:-100}"
# Identity of the tiny database, from tests/fixtures/phanta/tiny_db/README.md.
# Passing it on every run means these cases also exercise the fail-fast check
# on the happy path, not only in case 22.
PHANTA_DB_SHA="${PHANTA_DB_SHA:-554c323f5357f20cfd5372350676d1f55e3fc5a209015c0310c2180e7e5e299e}"

PHANTA_COMMON="--phanta_db ${PHANTA_DB}
               --phanta_read_length ${PHANTA_READ_LENGTH}
               --phanta_db_expected_inspect_sha256 ${PHANTA_DB_SHA}
               --log_versions false"

# ── Samplesheets ─────────────────────────────────────────────────────────────
# The committed sheets carry the token __TINY_INPUTS__ where the reads
# directory belongs, because a samplesheet holds absolute paths and nf_run
# launches each driver from its own directory, so a relative path would resolve
# against the wrong place. See tests/fixtures/phanta/samplesheets/README.md.
PHANTA_SHEETS="${RESULTS_BASE}/phanta_samplesheets"
rm -rf "${PHANTA_SHEETS}"; mkdir -p "${PHANTA_SHEETS}"
for sheet in "${PHANTA_FIXTURES}"/samplesheets/*.csv; do
    sed "s|__TINY_INPUTS__|${PHANTA_FIXTURES}/tiny_inputs|g" "${sheet}" \
        > "${PHANTA_SHEETS}/$(basename "${sheet}")"
done

# ── Phanta-specific helpers ──────────────────────────────────────────────────

# Compare a run's published tree against a committed baseline, byte for byte.
# bin/phanta_compare_outputs.sh does the work and excludes *.krak and
# pipeline_info/; see the migration spec, "Comparison model".
check_identical() {
    local name="$1" baseline="$2"
    local actual="${RESULTS_BASE}/${name}/phanta"
    if [ ! -d "${actual}" ]; then
        echo "[FAIL] ${name}: nothing published at ${actual}"
        ((fail++)); return
    fi
    if "${REPO_ROOT}/bin/phanta_compare_outputs.sh" "${baseline}" "${actual}" \
         > "${RESULTS_BASE}/${name}.compare.log" 2>&1; then
        echo "[PASS] ${name}: output byte-identical to $(basename "${baseline}") baseline"
        ((pass++))
    else
        echo "[FAIL] ${name}: output differs from $(basename "${baseline}") baseline"
        sed -n '1,30p' "${RESULTS_BASE}/${name}.compare.log" >&2
        ((fail++))
    fi
}

# Nextflow's own trace is the authority on what ran and what was cached. Its
# columns are task_id, hash, name, status, ... and a resumed run appends a
# second trace file to the same folder, so take the newest.
phanta_trace() {
    ls -t "${RESULTS_BASE}/$1/pipeline_info/execution_trace_"*.txt 2>/dev/null | head -1
}

# How many tasks of one process are in the trace, by exact process name.
# Barrier processes are untagged, so their name has no "(sample)" suffix.
phanta_task_count() {
    local trace="$1" procname="$2"
    [ -n "${trace}" ] && [ -f "${trace}" ] || { echo 0; return; }
    awk -F'\t' -v p="${procname}" 'NR>1 && $3 == p' "${trace}" | wc -l | tr -d ' '
}

# Did any Kraken2 task actually start? Scanning the work directory for a staged
# command is version-proof, where the wording of Nextflow's console output is
# not -- and these guard cases exist precisely to prove nothing expensive ran.
phanta_kraken2_tasks() {
    find "${WORK_BASE}/$1" -name '.command.sh' -exec grep -l '^kraken2 ' {} + 2>/dev/null \
        | wc -l | tr -d ' '
}

# One fail-fast case: run, expect a nonzero exit, expect a message, and expect
# no Kraken2 task to have been started.
phanta_guard() {
    local name="$1" expected="$2"; shift 2
    local work="${WORK_BASE}/${name}" launch="${WORK_BASE}/${name}_launch"
    rm -rf "${RESULTS_BASE}/${name}" "${work}" "${launch}"; mkdir -p "${launch}"
    local status
    ( cd "${launch}" || exit 1
      "${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
          --workflow phanta \
          --outdir "${RESULTS_BASE}/${name}" \
          -work-dir "${work}" \
          "$@" ) > "${RESULTS_BASE}/${name}.log" 2>&1
    status=$?

    local kraken2_tasks; kraken2_tasks="$(phanta_kraken2_tasks "${name}")"
    if [ "${status}" = "0" ]; then
        echo "[FAIL] ${name}: exited 0; the guard did not fire"
        ((fail++))
    # -- before the pattern: several of these messages start with "--phanta_",
    # which grep would otherwise read as an option of its own.
    elif ! grep -qF -- "${expected}" "${RESULTS_BASE}/${name}.log"; then
        echo "[FAIL] ${name}: failed, but without the expected message"
        echo "         expected to find: ${expected}" >&2
        grep -E "ERROR" -A4 "${RESULTS_BASE}/${name}.log" | head -12 >&2
        ((fail++))
    elif [ "${kraken2_tasks}" != "0" ]; then
        echo "[FAIL] ${name}: guard fired, but ${kraken2_tasks} Kraken2 task(s) had already run"
        ((fail++))
    else
        echo "[PASS] ${name}: fails fast with a clear message, no Kraken2 task"
        ((pass++))
    fi
}

note "Phanta: database ${PHANTA_DB}"
echo "  read length ${PHANTA_READ_LENGTH}, profile '${PROFILE}'"

# ── The runs, in parallel ────────────────────────────────────────────────────
note "Launching the Phanta runs"

# 17. The main cohort: three paired-end samples, in the samplesheet's own
#     deliberately non-alphabetical order, one of them Bracken-ineligible.
#     This one run is also what cases 19, 21, 23 and 24 are checked against.
nf_run test17_phanta_pe \
    --workflow phanta \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
    ${PHANTA_COMMON} &

# 18. The same three samples, R1 only. Not a copy of 17: single-end
#     classification leaves far more reads unassigned, so the baseline differs.
nf_run test18_phanta_se \
    --workflow phanta \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_se.csv" \
    ${PHANTA_COMMON} &

# 20. The same reads as 17 in alphabetical order. Same results, different
#     column order -- which is the only way to tell "columns follow the
#     samplesheet" from "columns happen to be sorted".
nf_run test20_phanta_ordering \
    --workflow phanta \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_alphabetical.csv" \
    ${PHANTA_COMMON} &

wait

# 24. Resume. Same work directory, same outdir, so the second run has something
#     to resume from; run after 17 rather than beside it.
note "Launching the Phanta resume run"
(
    cd "${WORK_BASE}/test17_phanta_pe_launch" || exit 1
    "${NF}" run "${REPO_ROOT}/main.nf" ${PROFILE} \
        --workflow phanta \
        --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
        ${PHANTA_COMMON} \
        --outdir "${RESULTS_BASE}/test17_phanta_pe" \
        -work-dir "${WORK_BASE}/test17_phanta_pe" \
        -resume
) > "${RESULTS_BASE}/test24_phanta_resume.log" 2>&1
echo "$?" > "${RESULTS_BASE}/test24_phanta_resume.status"

# ── Checks ───────────────────────────────────────────────────────────────────

note "Test 17: Phanta paired-end, tiny database"
check_run       test17_phanta_pe
check_files     test17_phanta_pe "counts.txt" "relative_read_abundance.txt" \
                                 "relative_taxonomic_abundance.txt" "total_reads.tsv" \
                                 "samples_that_failed_bracken.txt" \
                                 "*.krak.report.filtering_decisions.txt"
check_identical test17_phanta_pe "${PHANTA_BASELINE}/tiny_pe"

# The two tool contracts from the migration spec, "Two load-bearing upstream
# mechanisms". Both are asserted here because both are silent when broken: a
# six-column report produces a confusing pandas error much later, and a renamed
# Bracken input makes the second output vanish under a name nothing looks for.
kraken_report="$(find "${WORK_BASE}/test17_phanta_pe" -name 'sample_a.krak.report' | head -1)"
columns="$(awk -F'\t' 'NR==1{print NF; exit}' "${kraken_report}" 2>/dev/null || echo 0)"
if [ "${columns}" = "8" ]; then
    echo "[PASS] test17_phanta_pe: Kraken2 report has 8 columns (--report-minimizer-data)"
    ((pass++))
else
    echo "[FAIL] test17_phanta_pe: Kraken2 report has ${columns} columns, expected 8"
    ((fail++))
fi

# Bracken derives its second output's name from its INPUT filename, so this
# file existing is what proves the input was staged as {samp}.krak.report.filtered.
if [ -f "${RESULTS_BASE}/test17_phanta_pe/phanta/classification/sample_a.krak.report_bracken_species.filtered" ]; then
    echo "[PASS] test17_phanta_pe: Bracken's implicit second output is under its upstream name"
    ((pass++))
else
    echo "[FAIL] test17_phanta_pe: no sample_a.krak.report_bracken_species.filtered"
    ((fail++))
fi

note "Test 18: Phanta single-end, tiny database"
check_run       test18_phanta_se
check_identical test18_phanta_se "${PHANTA_BASELINE}/tiny_se"

note "Test 19: Bracken-ineligible sample (sample_c, within the test 17 cohort)"
# The contract is zero-byte products plus a listing, NOT a failed run: a sample
# with too few reads for any species line to reach filter_thresh is a
# biological outcome, so the cohort completes without it.
pe_class="${RESULTS_BASE}/test17_phanta_pe/phanta/classification"
ineligible_ok=true
for f in sample_c.krak.report.filtered.bracken.scaled sample_c.krak.report_bracken_species.filtered; do
    if [ ! -f "${pe_class}/${f}" ]; then
        echo "  MISSING: ${f}" >&2; ineligible_ok=false
    elif [ -s "${pe_class}/${f}" ]; then
        echo "  NOT EMPTY (should be a zero-byte substitute): ${f}" >&2; ineligible_ok=false
    fi
done
grep -qx 'sample_c' "${pe_class}/samples_that_failed_bracken.txt" 2>/dev/null \
    || { echo "  sample_c not listed in samples_that_failed_bracken.txt" >&2; ineligible_ok=false; }
# ...and it must not appear as a column in the merged tables.
if head -1 "${RESULTS_BASE}/test17_phanta_pe/phanta/final_merged_outputs/counts.txt" \
     | grep -q 'sample_c'; then
    echo "  sample_c appears as a column in counts.txt" >&2; ineligible_ok=false
fi
if $ineligible_ok; then
    echo "[PASS] test19_phanta_ineligible: zero-byte products, listed as failed, no column"
    ((pass++))
else
    echo "[FAIL] test19_phanta_ineligible: ineligible-sample contract broken"
    ((fail++))
fi

note "Test 20: output ordering follows the samplesheet, not the alphabet"
check_run test20_phanta_ordering
# Test 17's sheet is sample_b, sample_a, sample_c; test 20's is a, b, c. Same
# reads either way, so the columns must swap while the rows stay sorted.
hdr17="$(head -1 "${RESULTS_BASE}/test17_phanta_pe/phanta/final_merged_outputs/counts.txt" | cut -f3-)"
hdr20="$(head -1 "${RESULTS_BASE}/test20_phanta_ordering/phanta/final_merged_outputs/counts.txt" | cut -f3-)"
rows20="$(tail -n +2 "${RESULTS_BASE}/test20_phanta_ordering/phanta/final_merged_outputs/total_reads.tsv" | cut -f1)"
sorted20="$(printf '%s\n' "${rows20}" | sort)"
if [ "${hdr17}" = "$(printf 'sample_b\tsample_a')" ] \
   && [ "${hdr20}" = "$(printf 'sample_a\tsample_b')" ] \
   && [ "${rows20}" = "${sorted20}" ]; then
    echo "[PASS] test20_phanta_ordering: columns follow each samplesheet; total_reads rows sorted"
    ((pass++))
else
    echo "[FAIL] test20_phanta_ordering: ordering contract broken"
    echo "         test17 columns: ${hdr17}" >&2
    echo "         test20 columns: ${hdr20}" >&2
    echo "         test20 rows:    ${rows20}" >&2
    ((fail++))
fi

note "Test 21: database manifest"
manifest="${RESULTS_BASE}/test17_phanta_pe/pipeline_info/phanta_db_manifest.txt"
manifest_ok=true
for field in "inspect.out.sha256	${PHANTA_DB_SHA}" "kraken2.version" "bracken.version" \
             "pandas.version" "container" "param.phanta_read_length	${PHANTA_READ_LENGTH}"; do
    grep -qF "${field}" "${manifest}" 2>/dev/null \
        || { echo "  missing from the manifest: ${field}" >&2; manifest_ok=false; }
done
if $manifest_ok; then
    echo "[PASS] test21_phanta_manifest: database identity, tool versions and params recorded"
    ((pass++))
else
    echo "[FAIL] test21_phanta_manifest: manifest incomplete"
    ((fail++))
fi

note "Test 22: wrong database — fails before any Kraken2 task"
# The whole point of the manifest running first. The expected checksum is all
# digits on purpose: an earlier version of this check was silently skipped
# because Nextflow coerced such a value to a number and Groovy's ?: treats 0 as
# false. See the migration spec, Appendix F.
phanta_guard test22_phanta_wrong_db \
    "this is not the database this run was configured for" \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
    --phanta_db "${PHANTA_DB}" \
    --phanta_read_length "${PHANTA_READ_LENGTH}" \
    --phanta_db_expected_inspect_sha256 0000000000000000000000000000000000000000000000000000000000000000 \
    --log_versions false

note "Test 23: cohort barriers run exactly once each"
# The dataflow fans out, back in, out, and in again. A barrier that runs per
# sample instead of once would still produce output -- wrong output -- so the
# cardinality is asserted rather than inferred.
trace17="$(phanta_trace test17_phanta_pe)"
failed_tasks="$(phanta_task_count "${trace17}" 'PHANTA:phanta_failed_samples')"
totals_tasks="$(phanta_task_count "${trace17}" 'PHANTA:phanta_total_reads')"
if [ "${failed_tasks}" = "1" ] && [ "${totals_tasks}" = "1" ]; then
    echo "[PASS] test23_phanta_barriers: each cohort barrier ran exactly once"
    ((pass++))
else
    echo "[FAIL] test23_phanta_barriers: phanta_failed_samples ran ${failed_tasks} time(s), phanta_total_reads ${totals_tasks}"
    ((fail++))
fi

note "Test 24: resume"
resume_status="$(cat "${RESULTS_BASE}/test24_phanta_resume.status" 2>/dev/null || echo missing)"
trace_resume="$(phanta_trace test17_phanta_pe)"
cached=0; ran=0
if [ -n "${trace_resume}" ]; then
    cached="$(awk -F'\t' 'NR>1 && $4 == "CACHED"' "${trace_resume}" | wc -l | tr -d ' ')"
    ran="$(awk -F'\t' 'NR>1 && $4 != "CACHED"' "${trace_resume}" | wc -l | tr -d ' ')"
fi
if [ "${resume_status}" = "0" ] && [ "${ran}" = "0" ] && [ "${cached}" -gt 0 ]; then
    echo "[PASS] test24_phanta_resume: ${cached} task(s) cached, none re-run"
    ((pass++))
else
    echo "[FAIL] test24_phanta_resume: exit ${resume_status}, ${cached} cached, ${ran} re-run"
    # A launcher running in Nextflow's `mode=dev` (a jar started with java -cp)
    # does not persist the cache at all, and reports zero cached for every
    # pipeline. Check that before blaming the workflow.
    echo "         If this says 0 cached, confirm '${NF}' is an installed launcher," >&2
    echo "         not a jar run via java -cp -- mode=dev keeps no cache." >&2
    ((fail++))
fi

# ── Guards, which fail fast and need no cluster ──────────────────────────────
note "Test 25: input guards"

phanta_guard test25a_phanta_duplicate_ids \
    "duplicate sample 'sample_a'" \
    --phanta_samplesheet "${PHANTA_SHEETS}/bad_duplicate_ids.csv" ${PHANTA_COMMON}

phanta_guard test25b_phanta_mixed_layout \
    "mixes paired-end and single-end samples" \
    --phanta_samplesheet "${PHANTA_SHEETS}/bad_mixed_layout.csv" ${PHANTA_COMMON}

phanta_guard test25c_phanta_bad_header \
    "must start with the header" \
    --phanta_samplesheet "${PHANTA_SHEETS}/bad_header.csv" ${PHANTA_COMMON}

phanta_guard test25d_phanta_no_db \
    "--phanta_db is required" \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
    --phanta_read_length "${PHANTA_READ_LENGTH}" --log_versions false

# The parent of the database directory: a real path, with none of the files in
# it. Checks the manifest's per-file validation rather than a missing folder.
phanta_guard test25e_phanta_incomplete_db \
    "required database file missing or empty: hash.k2d" \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
    --phanta_db "${PHANTA_FIXTURES}/tiny_db" \
    --phanta_read_length "${PHANTA_READ_LENGTH}" --log_versions false

# Bracken's distribution is built per read length, and the wrong one is simply
# absent -- so the message has to name the parameter that must change.
phanta_guard test25f_phanta_read_length \
    "Set --phanta_read_length to one the database ships" \
    --phanta_samplesheet "${PHANTA_SHEETS}/tiny_pe.csv" \
    --phanta_db "${PHANTA_DB}" \
    --phanta_read_length 150 --log_versions false

phanta_guard test25g_phanta_no_samplesheet \
    "--phanta_samplesheet is required" \
    --phanta_db "${PHANTA_DB}" \
    --phanta_read_length "${PHANTA_READ_LENGTH}" --log_versions false
