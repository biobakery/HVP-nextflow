# Phanta Migration: Implementation Spec

What to build. Rationale, alternatives considered, and the analysis behind these
choices are in [phanta_migration_decisions.md](phanta_migration_decisions.md) —
read that only if you want to know *why*.

## Scope

Add Phanta taxonomic profiling as a selectable workflow in HVP-nextflow. All
code, tests, adapted scripts, config, and fixtures live in this repository. The
Phanta checkout is read-only: never write to it.

**In scope** — the path from FASTQ reads to four cohort tables:
`total_reads.tsv`, `counts.txt`, `relative_read_abundance.txt`,
`relative_taxonomic_abundance.txt`.

**Out of scope for v1** — Phanta post-processing scripts, the `single_end_krak`
prophage branch, the R merge variants, a second cloud provider, unifying the
samplesheet and `--readsdir` input models. Ordinary single-end samples are
supported through the main path.

## Pinned artifacts

| Item | Value |
| --- | --- |
| HVP-nextflow checkout | `/Users/jway/source/HVP-nextflow`, branch `jw_phanta` |
| Phanta checkout (read-only) | `/Users/jway/source/phanta` |
| Phanta commit | `ecbfcf9` on `main` |
| Phanta database | Select in phase 1; recommended `unmasked_db_v1.01.tar.gz`. Pin by filename + SHA256 computed at download. 20-25 GB |
| Phanta test reads | `http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/test_dataset.tar.gz` — four `.fq.gz`, **not** in the repo |
| Kraken2 | >= 2.1.2 (needs `--report-minimizer-data`); pin 2.1.6 |
| Bracken, Python, pandas, numpy | Pin exact versions in the Phanta container; record what the baseline resolved to |
| Upstream fixtures | `<phanta>/testing/{classification,final_merged_outputs}/` |

All database artifacts live at
`http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/`. Eight other builds exist
(masked, GTDB, UHGGv2, UHGV variants, MAGIC); v1 uses one, and `phanta_db`
makes switching a config change. The full catalogue is in the decisions doc.

### Required database files

Validate every one of these before any expensive work:

| File | Read by |
| --- | --- |
| `hash.k2d`, `taxo.k2d`, `opts.k2d`, `seqid2taxid.map` | Kraken2 |
| `databaseNmers.kmer_distrib` | Bracken — **specific to the configured read length** |
| `inspect.out` | `prep_to_merge_counts.py` |
| `taxonomy/nodes.dmp` | `filter_kraken_reports.py`, `process_filtered_kraken.py`, `prep_to_merge_counts.py` |
| `taxonomy/names.dmp` | taxonomy lookups |
| `library/species_genome_size.txt` | `scale_bracken.py` |

`host_prediction_to_genus.tsv` and `species_name_to_vir_score.txt` are
post-processing only — do not validate them.

On first download, confirm `library/species_genome_size.txt` ships in the main
tarball. The `library/` FASTAs are a separate download
(`database_V1.library.tar.gz`); if the size file only comes with that add-on, it
becomes mandatory and the staging footprint grows.

## Repository conventions

Three rules, each of which has bitten this repo before:

1. **Use `modules/` only.** The repo also has `processes/` (older style:
   separate `val(sample)` and `path(...)` inputs). Do not add to it.
2. **`meta.paired_end`, never `meta.single_end`.** Subworkflows carry
   `[ [id: sample, paired_end: bool], files ]` and unwrap to a plain sample
   string at the module boundary — see
   [subworkflows/quality_control.nf](../subworkflows/quality_control.nf) lines
   31-43. An inverted boolean here silently yields an empty channel.
3. **Phanta takes a samplesheet, not `--readsdir`.** This is a deliberate
   second input model for this repo. Cohort aggregation needs an explicit
   ordered sample list, which filename globbing cannot provide.

## Target layout

```text
HVP-nextflow/
├── main.nf                                 # Add phanta route; adjust --readsdir guard
├── nextflow.config                         # Phanta params (Appendix B)
├── template-params.yaml                    # Document Phanta params
├── workflows/phanta.nf
├── subworkflows/
│   ├── phanta_input.nf                     # Samplesheet -> ordered channel
│   └── phanta_profiling.nf                 # Per-sample chain + cohort barriers
├── modules/
│   ├── kraken2/main.nf
│   ├── bracken/main.nf
│   ├── taxonomy/
│   │   ├── phanta_filter/main.nf           # filter + Bracken eligibility
│   │   └── phanta_scale/main.nf
│   └── utils/
│       ├── phanta_db_manifest/main.nf      # Records + verifies database identity
│       ├── phanta_failed_samples/main.nf
│       ├── phanta_prepare_counts/main.nf   # counts + total reads
│       ├── phanta_normalize_counts/main.nf
│       └── phanta_merge_tables/main.nf
├── bin/
│   ├── phanta_convert_samplesheet.py       # Legacy TSV -> CSV
│   ├── phanta_compare_outputs.sh           # diff -r wrapper
│   └── phanta_*.py|sh                      # Adapted CLIs (Appendix A.2)
├── conf/
│   ├── base.config                         # Add Phanta withName blocks
│   ├── phanta.config                       # Pinned container, included by every profile
│   └── databases/<site>.config             # phanta_db + expected checksum per site
├── tests/
│   ├── phanta_workflow.nf.test
│   ├── phanta_modules.nf.test
│   └── fixtures/phanta/{baseline,samplesheets,tiny_db,tiny_inputs}/
└── docs/
    ├── phanta_migration_spec.md
    ├── phanta_migration_decisions.md
    └── phanta_testing.md                   # Baseline refresh + full-DB procedure
```

`conf/base.config` already exists with `withName` blocks and `check_max` —
extend it. `conf/databases/` already holds per-site database paths included by
profile; `phanta_db` belongs there alongside `metaphlan_db`.

## Process map

| Process | Upstream rule(s) | Scope | Outputs |
| --- | --- | --- | --- |
| `phanta_input` (subworkflow) | sample-file parsing | cohort | Ordered sample metadata + ordered sample-ID channel |
| `phanta_db_manifest` | — (new) | cohort, **before Kraken2** | `pipeline_info/phanta_db_manifest.txt`; fails on checksum mismatch |
| `kraken2` | `rule kraken` | per sample | Minimizer report; `.krak` (unpublished) |
| `phanta_filter` | `filter_kraken` + `process_filtered_kraken` | per sample | Filtered report; species table; filtering decisions; eligibility status; zero-byte substitutes when ineligible |
| `bracken` | `rule bracken` | per sample | Bracken table; Kraken-style species report (implicitly named) |
| `phanta_failed_samples` | `samples_failed_bracken` | **cohort barrier** | `samples_that_failed_bracken.txt`, sorted |
| `phanta_scale` | `scale_bracken` | per sample | Genome-length-scaled table |
| `phanta_prepare_counts` | `prepare_to_merge_counts` + `tot_reads_per_sample` | per sample, **needs failed list** | `.to_merge` table; per-sample count row |
| `phanta_total_reads` | `tot_reads_all` | **cohort barrier** | `total_reads.tsv`, sorted rows |
| `phanta_normalize_counts` | `prepare_to_merge_normed` | per sample, **needs `total_reads.tsv`** | `.norm_brack` per sample |
| `phanta_merge_tables` | `merge_counts_normed` + `merge_scaled_bracken` | cohort | Three merged tables, columns in samplesheet order |

Do not split the folded pairs (filter+eligibility, prepare-counts+total-reads)
back into one process per upstream rule.

### Cohort barriers

The dataflow is **not** a straight per-sample chain. Two cohort-level files feed
back into per-sample steps:

1. `phanta_prepare_counts` needs the cohort `samples_that_failed_bracken.txt`.
2. `phanta_normalize_counts` needs the cohort `total_reads.tsv`.

Shape: fan-out → fan-in → fan-out → fan-in → fan-out → fan-in. Implement each
barrier in `subworkflows/phanta_profiling.nf` as an explicit `collect()` then
`combine` back onto the per-sample channel. Assert in nf-test that each barrier
process runs exactly once per run.

Carry `[ [id:, paired_end:], files... ]` between per-sample steps. Never
recover sample IDs by parsing filenames. Pass an explicit ordered sample list
into `phanta_merge_tables`.

### Two load-bearing upstream mechanisms

**Bracken eligibility is implemented as substitute outputs, not a status flag.**
`process_filtered_kraken.py` decides Bracken will fail when no line of the
filtered report has `rank == 'S'` and `assigned_reads >= filter_thresh`. It then
writes three zero-byte files:

```
{samp}.krak.report.filtered.bracken.temp
{samp}.krak.report_bracken_species.filtered.temp
{samp}.krak.report.filtered.bracken.scaled.temp
```

`rule bracken` and `rule scale_bracken` `cp` these into place *instead of
running the tool*, and the failed-sample list is built from the existence of the
third. The observable contract for an ineligible sample: Bracken table,
Kraken-style species report, scaled table, and `.norm_brack` are all **empty
files**. Replace the canary with an explicit status/reason record, but still emit
byte-identical zero-length products — the fixtures contain them.

**Bracken's second output filename is implicit.** `rule bracken` passes only
`-o`; the Kraken-style report appears because Bracken derives its name from the
*input* filename:

```
input    {samp}.krak.report.filtered
-o       {samp}.krak.report.filtered.bracken
implied  {samp}.krak.report_bracken_species.filtered
```

Stage the filtered report under exactly `{samp}.krak.report.filtered`. Any
wrapper that renames its input breaks this.

## Output ordering contract

Nextflow channel arrival order is nondeterministic. Force order at emit time,
per file:

| Output | Order |
| --- | --- |
| `counts.txt`, `relative_read_abundance.txt`, `relative_taxonomic_abundance.txt` | Sample columns in samplesheet order |
| `total_reads.tsv` | Sample rows sorted by name |
| `samples_that_failed_bracken.txt` | Sorted by name |
| Taxon rows in merged tables | As `pd.concat(axis=1)` produces them from inputs supplied in samplesheet order |

Note the two different orders: `total_reads.tsv` is sorted while the merged
tables follow the samplesheet. This matches upstream. There is no
canonicalization step — ordering is a property of the emitting process.

## Provenance

Every published run must record the database that produced it. Upstream does
not, which is why its own default-database ambiguity cannot be resolved from
existing outputs.

**The manifest is a sidecar.** Do not add a provenance header, comment, or
column to the four tables — they must stay byte-identical to the baseline, and
per-run provenance would break that.

`modules/utils/phanta_db_manifest` writes
`${params.outdir}/pipeline_info/phanta_db_manifest.txt`, one `key\tvalue` per
line in fixed order, recording:

| Field | Notes |
| --- | --- |
| `phanta_db` path | As given plus `realpath` |
| `hash.k2d` | Size always; SHA256 only when `phanta_db_checksum_full` |
| `taxo.k2d`, `opts.k2d`, `seqid2taxid.map` | Size + SHA256 |
| `inspect.out` | Size + SHA256 — **primary identity field** |
| `library/species_genome_size.txt` | Size + SHA256 |
| `databaseNmers.kmer_distrib` | Filename (records read length) + SHA256 |
| Kraken2, Bracken versions | From `--version` |
| Python, pandas, numpy versions | Byte-identity relevant |
| Container image digest | The digest, not the tag |
| Phanta params | Every Appendix B value used |

`inspect.out` is the fingerprint that matters: a text dump of every taxon with
clade k-mer counts, so it is content-derived and cheap to hash. `hash.k2d` is
20+ GB, so its full checksum is gated.

**Fail fast.** When `phanta_db_expected_inspect_sha256` is set, a mismatch
fails the run **before any Kraken2 task**. Set it in `conf/databases/` next to
`phanta_db`.

The manifest contains a path and timestamp, so it is not byte-stable. It lives
in `pipeline_info/`, which the comparison excludes; assert its checksum fields
in nf-test instead.

## Publishing contract

Replaces upstream's `deal_with_intermediate` / `delete_intermediate`.

**Always publish:** the four tables, the failed-sample report,
`pipeline_info/phanta_db_manifest.txt`, `pipeline_log.txt`, filtering
decisions, Kraken-style Bracken species reports, scaled Bracken tables. The
manifest is never optional.

**Never publish by default:** raw `.krak`, raw or filtered Kraken reports,
standard Bracken intermediates. Retain in the work directory for `-resume` per
storage lifecycle policy.

`--save_intermediate` publishes stable intermediates without moving or deleting
task files. Document that deleting the work directory disables resume.

Decide explicitly whether to emit `pipeline_completed.txt` — it is in upstream's
`rule all`, so a directory comparison will flag it either way.

## Comparison model

Three distinct operations:

1. **Baseline refresh** (rare maintenance) — run upstream Phanta at the pinned
   commit, in the pinned container, against the pinned database, writing outside
   the Phanta checkout. Review and update `tests/fixtures/phanta/baseline/`. Only
   when the commit, database, tool versions, inputs, params, or ordering
   contract intentionally change.
2. **Routine fixture regression** (CI/development) — HVP-nextflow only,
   `diff -r` against the committed baseline. Never reruns upstream.
3. **Paired-run comparison** (validating a new dataset) — upstream and
   HVP-nextflow from identical inputs, params, database, tool versions, and
   container into two directories, then `diff -r`. Do not commit the outputs.

`bin/phanta_compare_outputs.sh` wraps `diff -r`: examines the whole tree rather
than stopping at the first mismatch, lists every missing/unexpected/differing
file, excludes `.krak` and `pipeline_info/`, exits nonzero on any difference.
Investigate every mismatch; no numeric tolerances.

## Work plan

**16-22 engineer-days over 4-5 weeks.** Owners are roles until assigned.

### Phase 0 — Prerequisites: container and test database (3-4 days)

*Owner: pipeline engineer; container reviewed by platform engineer.*

Phases 1 and 3 cannot run without these.

1. Build and publish the pinned Phanta container: Kraken2 2.1.6, Bracken,
   exactly pinned Python/pandas/numpy. Reference it from `conf/phanta.config`,
   included by every profile.
2. Build a tiny Kraken2 + Bracken database for CI per [Appendix D](#appendix-d-building-the-tiny-ci-database).
   The full database is 20-25 GB and GitHub runners have ~14 GB of disk.
3. Obtain Phanta's test reads (separate download — record URL and checksum).
   `testing/samp_file.txt` references paths absent from a fresh checkout.
4. Establish whether CI currently passes.
   [.github/workflows/ci-tests.yml](../.github/workflows/ci-tests.yml) runs
   `nf-test test --ci`, and [tests/main.nf.test](../tests/main.nf.test) loads
   `tuftshpc-params.yaml`, which needs site databases a runner lacks. If CI is
   already red, scope the Phanta tests to run standalone.

**Deliverables:** published pinned container; tiny database with its build
script and checksum; recorded test-read URL and checksum; known-good CI starting
state.

**Gate:** `kraken2` and `bracken` complete against the tiny database, in the
pinned container, on a GitHub runner within its limits.

### Phase 1 — Baseline and deployment target (2-3 days)

*Owner: pipeline engineer; baseline reviewed by domain scientist; target
confirmed with platform/data owner.*

1. Decide which profiles Phanta supports at release — at minimum one HPC profile
   plus `amazon`. The repo already ships `standard`/`local`, `tufts_hpc`,
   `harvard_rc`, `amazon` (AWS Batch, S3 work dir), and legacy `engaging`. Base
   the choice on where the FASTQs and the 25 GB database already live.
2. Select and pin the database first. Choose one artifact, have the domain
   scientist confirm it suits HVP's sample types, compute SHA256 at download,
   record filename/checksum/date in `docs/phanta_testing.md`, and set
   `phanta_db_expected_inspect_sha256` per site.
3. Run upstream Phanta at commit `ecbfcf9` on the test dataset, in the pinned
   container, against the pinned database. Record commands, params, and resolved
   tool/pandas/numpy versions. Procedure in [Appendix C](#appendix-c-running-upstream-phanta).
4. Run it again independently and confirm byte-identical output. If anything but
   `.krak` differs, stop and re-derive the ordering contract.
5. Commit the reviewed outputs to `tests/fixtures/phanta/baseline/` with
   checksums. Record which order each file carries.
6. Write `bin/phanta_compare_outputs.sh` with unit tests. Small script, not a
   framework.
7. Produce and commit a tiny-database baseline so the routine regression can run
   without the full database.

**Deliverables:** confirmed profile targets; pinned database filename, SHA256,
and selection rationale; upstream revision and tool manifest including resolved
pandas/numpy; full-database and tiny-database baselines with checksums;
per-file ordering contract recorded; comparison wrapper with unit tests;
baseline-refresh instructions in `docs/phanta_testing.md`.

**Gate:** database pinned by checksum; two independent upstream runs
byte-identical apart from `.krak`; a clean checkout compares fixtures without
rerunning Phanta.

### Phase 2 — Port and database staging (7-9 days)

*Owner: Nextflow engineer; platform engineer owns staging; domain scientist
reviews Phanta-specific transformations.*

1. Add `phanta` to the router in [main.nf](../main.nf); remove `phanta` from the
   `--readsdir` guard and add a `--phanta_samplesheet` guard. Implement
   `workflows/phanta.nf`, `subworkflows/phanta_input.nf`,
   `subworkflows/phanta_profiling.nf`. No second root pipeline.
2. Add every parameter in [Appendix B](#appendix-b-parameter-mapping) to
   [nextflow.config](../nextflow.config) and
   [template-params.yaml](../template-params.yaml); add `phanta_db` and
   `phanta_db_expected_inspect_sha256` to each `conf/databases/<site>.config`.
   Note in `template-params.yaml` that `paired_end`/`single_end`/`filepattern`
   do not apply to Phanta.
3. Validate before expensive work: duplicate sample IDs, mixed library layout,
   missing read files, parameter ranges, and every file under "Required database
   files".
4. Implement `phanta_db_manifest` and run it as part of validation, ahead of the
   first Kraken2 task. Give it a `withName` block with generous time —
   `phanta_db_checksum_full` reads 20+ GB.
5. Spend at most one engineer-day on an nf-core module fit check:
   - `KRAKEN2_KRAKEN2`: accepts `--report-minimizer-data` via `ext.args`,
     retains the eight-column report, suppresses per-read output by default,
     reports >= 2.1.2.
   - `BRACKEN_BRACKEN`: reproduces `-l S`, configured `-r` and `-t`, emits both
     outputs, and **stages its input as `{samp}.krak.report.filtered` without
     renaming** (see "Two load-bearing upstream mechanisms").
   - nf-core modules carry their own container directives, which conflict with
     the single pinned image. If reconciling costs more than the check saves,
     write local modules under `modules/kraken2/` and `modules/bracken/`.
     Upstream reuse is not a release blocker.
   Capture both contracts in nf-tests either way.
6. Establish full-database staging per target profile. HPC: shared filesystem
   plus node-local scratch. AWS Batch: compare ordinary task localization
   against download-once reuse on node-local NVMe over a few Kraken2 tasks.
   Record staging time and bytes transferred; pick the simplest adequate option.
7. Copy the scientific scripts into flat `bin/phanta_*` files and adapt them
   there — never edit upstream copies. Signatures and current behavior in
   [Appendix A.2](#a2-script-signatures). Every script takes explicit input and
   output paths. Specifically:
   - Filtering takes explicit species, filtered-report, and decision paths, and
     emits a status/reason file instead of `.temp` canaries.
   - Count normalization processes one sample per invocation instead of scanning
     a directory.
   - Total-read calculation stays shell but takes declared Kraken and Bracken
     file paths instead of a directory. Preserve `exit 64` on no reads.
   - Scaling, count preparation, and merging take explicit output paths and
     ordered input manifests.
   - Port `merging_scaled_bracken.py` only; drop the two R variants (neither is
     wired into the `Snakefile`).
8. Use isolated task directories. Do not set global `stageInMode 'copy'`. No
   process may scan a shared classification directory or derive an output beside
   a staged input.
9. Implement both cohort barriers explicitly.
10. Route Bracken-ineligible samples: zero-byte products plus a status/reason
    record, and a sorted `samples_that_failed_bracken.txt`. Tool and
    infrastructure failures stay fatal.
11. Implement the publishing contract above.
12. Do not port `single_end_krak` or post-processing.

**Deliverables:** routed `phanta` workflow; validated samplesheet input channel;
Kraken2 and Bracken modules sharing the pinned container; database manifest with
fail-fast checking; task-isolated scripts and modules; both cohort barriers;
deterministic aggregation; publishing contract implemented; one working
database-staging path per target profile, with recorded staging time and bytes
transferred.

**Gate:** end-to-end run against the tiny database produces all four tables;
every process declares inputs and outputs; no process scans a shared results
directory; module tests pass.

### Phase 3 — Validate (3-4 days)

*Owner: pipeline engineer; acceptance by domain scientist; cluster and cloud
runs by platform engineer.*

1. Add nf-tests under [tests/](../tests/) per
   [Appendix E.2](#e2-test-cases). There is no existing per-module test
   precedent — `tests/` holds two pipeline-level tests — so this establishes
   the convention. Keep it consistent with `nf-test.config`
   (`testsDir "tests"`).
2. Run the routine fixture regression against the tiny database.
3. Verify `-resume` locally: run the same test twice, confirm scientific tasks
   are cached.
4. Run the full-database workflow on each target profile. Require byte-identical
   output against the full-database baseline.
5. Paired-run comparison on at least one representative real cohort. Retain the
   report; obtain domain-scientist signoff.
6. Keep the Phanta checkout runnable as the comparator through this milestone.

**Deliverables:** passing module and workflow nf-tests; tiny-database regression
report; local resume evidence; full-database comparison report per target
profile; real-cohort paired-run report and domain-scientist signoff.

**Gate:** tiny-database regression, full-database runs, and the real-cohort
comparison all byte-identical apart from excluded files.

### Phase 4 — Tune (1-2 days)

*Owner: platform engineer with pipeline engineer.*

1. Targeted refinements to the chosen staging approach only.
2. Benchmark two or three Kraken2 CPU/memory shapes; set one default plus
   bounded OOM retry escalation in [conf/base.config](../conf/base.config),
   reusing the existing `check_max`/`task.attempt` and `errorStrategy`
   exit-status patterns.
3. One Spot/preemptible interruption test on AWS; verify `-resume` restarts the
   task without corrupting or duplicating published output.
4. Measure merge peak memory and runtime once on the representative cohort.
   Record the ceiling; build sharded aggregation only if it proves a real limit.
5. Add a README section and a `docs/` runbook: launch, the samplesheet model and
   how it differs from `--readsdir`, database staging, publishing,
   work-directory retention, resume, outputs, measured defaults.

**Deliverables:** resource defaults in `conf/base.config`; staging decision and
measurements; interruption/resume result; merge memory measurement; README
section and user runbook.

**Gate:** each target profile has documented measured defaults; AWS recovers
from one forced interruption.

## Acceptance criteria

1. No migration files, tracked or untracked, in the Phanta repository.
2. One database artifact selected, pinned by filename and SHA256, confirmed
   suitable by the domain scientist, used for both baseline and HVP runs.
3. Every published run records its database: `phanta_db_manifest.txt` with the
   `inspect.out` checksum, tool versions, and container digest. A run against a
   database whose checksum does not match
   `phanta_db_expected_inspect_sha256` fails before any Kraken2 task.
4. Two independent upstream baseline runs byte-identical apart from `.krak`;
   routine regression uses committed fixtures without rerunning Phanta.
5. Tiny database and pinned container exist; fixture regression runs on a
   GitHub runner without the 25 GB database.
6. All Phanta processes run in the same pinned image on every profile.
7. Module-fit tests prove Kraken2 minimizer compatibility and both Bracken
   outputs, or equivalent tests cover local fallbacks.
8. Both cohort barriers implemented explicitly, covered by cardinality tests.
9. No module writes outside its task directory; all outputs and
   biological-empty behavior explicit.
10. All four tables byte-identical to baseline against both tiny and full
    databases, on every target profile.
11. A representative real cohort passes comparison with domain-scientist
    signoff.
12. `-resume` passes locally and on cluster, including one AWS Spot
    interruption.
13. Database staging, Kraken2 sizing, and one aggregation-memory measurement
    documented without speculative scale tiers.

## Decisions already made

Do not relitigate these during implementation. Rationale for each is in
[phanta_migration_decisions.md](phanta_migration_decisions.md).

- All implementation lives in HVP-nextflow; Phanta stays read-only.
- Integrate through the existing router and `modules/`; never add to
  `processes/`; no separate nf-core repository.
- Prefer pinned nf-core Kraken2 and Bracken modules with a one-day fit limit and
  local fallbacks; the shared container overrides any module's own container
  directive.
- Phanta takes a samplesheet; the other workflows keep `--readsdir`. A legacy
  TSV converter ships alongside.
- Fix output ordering per file at emit time. No canonicalization layer; compare
  with `diff -r` excluding `.krak` and `pipeline_info/`.
- The pinned container and tiny test database are prerequisites, not deferred
  work.
- Keep `calc_total_reads.sh` in shell for v1.
- Port `merging_scaled_bracken.py`; drop the R merge variants.
- Fold filter with eligibility, and prepare-counts with total-reads.
- Replace destructive intermediate cleanup with the explicit publishing
  contract.
- Every run records its database identity; a mismatch fails before Kraken2.

Deferred, not denied: Phanta post-processing, `single_end_krak`, speculative
scale tiers, sharded aggregation, a second cloud provider, and unifying the
samplesheet and `--readsdir` input models.

## Open questions

Resolve with the named owner; record the answer here.

1. **Which database, and is it suitable?** `README.md` and `databases.md` name
   different artifacts (`database_V1.tar.gz` vs `unmasked_db_v1.01.tar.gz`) and
   upstream publishes no checksums. The latter is recommended. The open part is
   scientific: Phanta's default database is human-gut specific, so confirm it
   fits HVP's sample types. *Domain scientist confirms; engineer records
   filename, SHA256, date. Phase 1.*
2. **Does CI currently pass?** See phase 0 step 4. *Pipeline engineer, phase 0.*
3. **Which profiles ship at release?** *Platform/data owner, phase 1.*
4. **Tiny database committed or downloaded?** Depends on achieved size.
   *Pipeline engineer, phase 0.*
5. **Which container registry?** `aws.config` uses `public.ecr.aws`; HPC
   profiles mostly use environment modules. The pinned image needs a home
   reachable from every target executor. *Platform engineer, phase 0.*
6. **Emit `pipeline_completed.txt` or exclude it?** *Pipeline engineer, phase 2.*
7. **Which real cohort is the acceptance dataset?** Needed before phase 3 step
   5. *Domain scientist.*

---

# Appendix A: upstream contracts

From `Snakefile` at commit `ecbfcf9`. Every upstream script takes **bare
positional `sys.argv`** and **derives its output paths from its input path** —
no script accepts an output path. Adapting them to explicit arguments is a real
code change, and it is what makes task isolation possible.

## A.1 Tool command lines

Kraken2 (`rule kraken`), threads from `class_threads`:

```bash
kraken2 --db $DB --threads $T --output {samp}.krak \
  --report {samp}.krak.report --report-minimizer-data \
  [--paired] [--gzip-compressed] $READS --confidence $CONFIDENCE
```

`--paired` only for paired cohorts; `--gzip-compressed` only when `gzipped:
True`. Both are cohort-wide upstream, not per sample.

Bracken (`rule bracken`), 1 thread, 8 GB, 1 h upstream:

```bash
bracken -d $DB -i {samp}.krak.report.filtered \
  -o {samp}.krak.report.filtered.bracken -r $READ_LENGTH -l 'S' -t $FILTER_THRESH
```

`filter_thresh` does double duty: Bracken's `-t` **and** the eligibility
predicate in `process_filtered_kraken.py`. One param must drive both.

## A.2 Script signatures

Positional, in exactly this order.

| Upstream script | Positional arguments | Writes (derived from input) |
| --- | --- | --- |
| `filter_kraken_reports.py` | `krak_report`, `db`, `cov_thresh_bacterial`, `cov_thresh_viral`, `minimizer_thresh_bacterial`, `minimizer_thresh_viral`, `cov_thresh_arc`, `cov_thresh_euk`, `minimizer_thresh_arc`, `minimizer_thresh_euk` | `<krak_report>.species`, `.filtered`, `.filtering_decisions.txt` |
| `process_filtered_kraken.py` | `krak_report_filtered`, `filter_thresh`, `db` | up to three zero-byte `.temp` files |
| `scale_bracken.py` | `db`, `bracken_report`, `filtering_decisions`, `read_length`, `paired` | `<bracken_report>.scaled` |
| `prep_to_merge_counts.py` | `brack_report`, `db`, `failed_file`, `samp` | `<brack_report>.to_merge` |
| `prep_to_merge_normed.py` | `tot_reads_file`, `classification_dir`, `failed_file` | `<counts_table>.norm_brack` for **every** sample (scans the directory) |
| `calc_total_reads.sh` | `samp`, `krak_brack_dir`, `out_file` | the given out_file |
| `merge.py` | `--input-file A B C...`, `--output-file F`, `--threads N` | `F` |
| `merging_scaled_bracken.py` | `--db D`, `--input-file A B C...`, `--output-file F`, `--threads N` | `F` |

`filter_kraken_reports.py`'s thresholds are **not** grouped by domain —
bacterial and viral coverage, then bacterial and viral minimizers, then archaeal
and eukaryotic coverage, then archaeal and eukaryotic minimizers. Transcribing
in the intuitive order silently swaps thresholds.

`prep_to_merge_normed.py` is the one script whose contract genuinely changes. It
currently reads sample names from `total_reads.tsv` and builds
`<indir>/<samp>.krak.report_bracken_species.filtered.to_merge` for each. The
port processes one sample per invocation, taking that sample's counts table, its
Bracken-classified read count (field 4, `Assigned_Step_Three`) from the cohort
`total_reads.tsv`, and the failed-sample list. Preserve two edge cases: an empty
read count is treated as `0`, and a failed sample yields an empty output.
Normalization is `round(reads / bracken_classified_reads, 10)` formatted with
plain `str()` — not pandas — so reproduce it exactly.

## A.3 Filenames

Load-bearing: fixtures, Bracken's implicit naming, and the failed-sample glob
all depend on them.

```
classification/{samp}.krak                                        # raw, excluded
classification/{samp}.krak.report                                 # minimizer report
classification/{samp}.krak.report.species                         # temp upstream
classification/{samp}.krak.report.filtered
classification/{samp}.krak.report.filtering_decisions.txt
classification/{samp}.krak.report.filtered.bracken                # Bracken -o
classification/{samp}.krak.report_bracken_species.filtered         # Bracken implicit
classification/{samp}.krak.report.filtered.bracken.scaled
classification/{samp}.krak.report_bracken_species.filtered.to_merge
classification/{samp}.krak.report_bracken_species.filtered.to_merge.norm_brack
classification/{samp}_total_reads.txt
classification/samples_that_failed_bracken.txt
processed_filtered_kraken/{samp}.txt                              # dependency marker
final_merged_outputs/{total_reads.tsv,counts.txt,relative_read_abundance.txt,relative_taxonomic_abundance.txt}
pipeline_completed.txt
```

## A.4 Headers and formats

`total_reads.tsv`, written by `rule tot_reads_all` as a literal `echo -e` then
`cat`:

```
Samp_Name	Tot_Samp_Reads	Unassigned_Step_One	Assigned_Step_Three	Unassigned_Step_Three
```

`Tot_Samp_Reads` = Kraken root + unclassified. `Assigned_Step_Three` = sum of
column 6 of the filtered Bracken table. `Unassigned_Step_Three` = the
difference. `calc_total_reads.sh` exits **64** with "No reads in sample" when
both root and unclassified are absent.

`{samp}.krak.report.filtering_decisions.txt`:

```
species_taxid	superkingdom	max_cov	max_uniq_minimizers	kept
```

The three merged tables:

```
Taxon_Lineage_with_Names	Taxon_Lineage_with_IDs	<sample columns...>
```

Counts appear as pandas floats (`6564345.0`) because `merge.py` does
`pd.concat(...).fillna(0)` then `to_csv`. Do not "fix" this to integers.

Zero-byte files for an ineligible sample: the `.bracken`,
`_bracken_species.filtered`, `.scaled`, and `.norm_brack` outputs.

---

# Appendix B: parameter mapping

Upstream `config.yaml` key → HVP param → default. Defaults are upstream's.

| Upstream key | HVP param | Default |
| --- | --- | --- |
| `sample_file` | `phanta_samplesheet` | `null` (required) |
| `database` | `phanta_db` | `null`, per site in `conf/databases/` |
| `outdir` | `outdir` (existing) | `results` |
| `pipeline_directory` | — | dropped; scripts live in `bin/` |
| `confidence_threshold` | `phanta_confidence` | `0.1` |
| `gzipped` | `phanta_gzipped` | `true` |
| `class_threads` | `withName: kraken2` cpus | `16` |
| `class_mem_mb` | `withName: kraken2` memory | `32768` → `32.GB` |
| `single_end_krak` | — | not ported |
| `cov_thresh_bacterial` | `phanta_cov_thresh_bacterial` | `0.01` |
| `cov_thresh_viral` | `phanta_cov_thresh_viral` | `0.10` |
| `cov_thresh_arc` | `phanta_cov_thresh_arc` | `0.01` |
| `cov_thresh_euk` | `phanta_cov_thresh_euk` | `0` |
| `minimizer_thresh_bacterial` | `phanta_minimizer_thresh_bacterial` | `0` |
| `minimizer_thresh_viral` | `phanta_minimizer_thresh_viral` | `0` |
| `minimizer_thresh_arc` | `phanta_minimizer_thresh_arc` | `0` |
| `minimizer_thresh_euk` | `phanta_minimizer_thresh_euk` | `0` |
| `read_length` | `phanta_read_length` | `150` |
| `filter_thresh` | `phanta_filter_thresh` | `10` |
| `delete_intermediate` | replaced by `save_intermediate` | `false` |
| — | `run_phanta` | `true` |
| — | `phanta_db_expected_inspect_sha256` | `null`; per site. Mismatch fails the run |
| — | `phanta_db_checksum_full` | `false`; `true` for baselines and release validation |

Library layout comes from the samplesheet, so Phanta ignores the existing
`paired_end`/`single_end`/`filepattern` params. Say so in
`template-params.yaml`.

---

# Appendix C: running upstream Phanta

Needed in phases 0, 1, and 3. Never write into the Phanta checkout.

```bash
# one-time env
conda env create -f /Users/jway/source/phanta/phanta_env.yaml
conda activate phanta-envs

# copy the config out of the checkout, then edit the four paths
cp /Users/jway/source/phanta/config.yaml /tmp/phanta-baseline/config.yaml
# set: pipeline_directory, sample_file, outdir (outside the checkout), database

snakemake -s /Users/jway/source/phanta/Snakefile \
  --configfile /tmp/phanta-baseline/config.yaml \
  --jobs 99 --cores 1 --max-threads 16
```

Requires at least 32 GB memory. Upstream also ships `run_phanta.py`, a wrapper
with `-i/--input-dir`, `-s/--sample-sheet`, `-o/--output-dir` and threshold
flags that generates a config and calls Snakemake. Use the explicit Snakemake
invocation for baselining so the config is a committable artifact.

`testing/samp_file.txt` contains `/full/path/to/cloned/repo` placeholders and
lists `sample19` before `sample18`. That non-alphabetical order is what produces
the two distinct output orderings — keep it when reproducing the baseline.

---

# Appendix D: building the tiny CI database

No upstream recipe exists. Target under ~500 MB.

1. Choose small reference genomes spanning **all four** threshold classes —
   bacterial, viral, archaeal, eukaryotic — so every
   `cov_thresh_*`/`minimizer_thresh_*` branch in `filter_kraken_reports.py` is
   exercised. A viral genome plus one small bacterium, one archaeon, and one
   small eukaryotic contig suffices.
2. Build the Kraken2 index:
   ```bash
   kraken2-build --download-taxonomy --db tiny_db   # or a trimmed nodes/names pair
   kraken2-build --add-to-library <each>.fna --db tiny_db
   kraken2-build --build --db tiny_db
   kraken2-inspect --db tiny_db > tiny_db/inspect.out
   ```
3. Build the Bracken distribution for the CI read length (use a short one, e.g.
   100, to keep synthetic reads small):
   ```bash
   bracken-build -d tiny_db -t 4 -k 35 -l 100
   ```
   Yields `database100mers.kmer_distrib`; `phanta_read_length` must match.
4. Generate `library/species_genome_size.txt` with upstream's
   `pipeline_scripts/calculate_genome_size.py /path/to/tiny_db`. It needs
   Biopython, which `phanta_env.yaml` does not install — add it to the build
   environment.
5. Verify `taxonomy/nodes.dmp`, `taxonomy/names.dmp`, `seqid2taxid.map`, and
   `inspect.out` are present.
6. Create synthetic reads: two eligible samples, and **one deliberately
   Bracken-ineligible sample** whose filtered report has no species line
   reaching `filter_thresh` — cheapest route is very few reads. That sample
   exercises the zero-byte substitutes and the failed-sample list.
7. Keep the build as a committed script
   (`tests/fixtures/phanta/tiny_db/build.sh`) with a checksum, so the database
   is reproducible rather than an opaque blob.

If the result exceeds what should live in git, host it and have CI download by
pinned URL plus checksum.

---

# Appendix E: deliverables and test matrix

## E.1 Files to create or modify

| Path | Action | Done when |
| --- | --- | --- |
| `main.nf` | modify | `phanta` routes; `--readsdir` guard excludes it; `--phanta_samplesheet` guard added |
| `nextflow.config` | modify | all Appendix B params present with defaults |
| `template-params.yaml` | modify | Phanta params documented, including that layout params do not apply |
| `conf/phanta.config` | create | pinned container attached to all Phanta processes; included by every profile |
| `conf/base.config` | modify | `withName` blocks for `kraken2`, `bracken`, `phanta_db_manifest`, and the small table steps |
| `conf/databases/*.config` | modify | `phanta_db` + `phanta_db_expected_inspect_sha256` per site |
| `workflows/phanta.nf` | create | validates, calls subworkflows, publishes |
| `subworkflows/phanta_input.nf` | create | emits `[[id:,paired_end:],files]` plus an ordered sample-ID channel |
| `subworkflows/phanta_profiling.nf` | create | per-sample chain with both cohort barriers explicit |
| `modules/kraken2/main.nf` | create | A.1 command; minimizer report; `.krak` unpublished |
| `modules/bracken/main.nf` | create | A.1 command; input staged as `.krak.report.filtered`; both outputs emitted |
| `modules/taxonomy/phanta_filter/main.nf` | create | filtered report, species table, decisions, eligibility status |
| `modules/taxonomy/phanta_scale/main.nf` | create | scaled table; empty passthrough when ineligible |
| `modules/utils/phanta_db_manifest/main.nf` | create | writes the manifest; fails on checksum mismatch; runs before Kraken2 |
| `modules/utils/phanta_failed_samples/main.nf` | create | sorted failed list |
| `modules/utils/phanta_prepare_counts/main.nf` | create | `.to_merge` plus the per-sample total-reads row |
| `modules/utils/phanta_normalize_counts/main.nf` | create | one sample per task |
| `modules/utils/phanta_merge_tables/main.nf` | create | three tables, columns in samplesheet order |
| `bin/phanta_convert_samplesheet.py` | create | legacy TSV → CSV |
| `bin/phanta_*.py`, `bin/phanta_calc_total_reads.sh` | create | explicit I/O paths; no derived output names |
| `bin/phanta_compare_outputs.sh` | create | `diff -r` wrapper; excludes `.krak` and `pipeline_info/`; nonzero on any diff |
| `tests/phanta_modules.nf.test` | create | E.2 module cases pass |
| `tests/phanta_workflow.nf.test` | create | E.2 workflow cases pass |
| `tests/fixtures/phanta/**` | create | tiny DB, inputs, samplesheets, baselines |
| `docs/phanta_testing.md` | create | baseline refresh, full-DB procedure, resolved versions, checksums |
| `README.md` | modify | Phanta section |

## E.2 Test cases

| Case | Level | Asserts |
| --- | --- | --- |
| Paired-end cohort, tiny DB | workflow | four tables byte-identical to fixture |
| Ordinary single-end cohort | workflow | same, single-end channel path |
| Bracken-ineligible sample | workflow | zero-byte products; sample in `samples_that_failed_bracken.txt`; run succeeds |
| Duplicate sample IDs | workflow | fails before any Kraken2 task |
| Mixed paired/single cohort | workflow | fails with a clear message |
| Missing database file | workflow | fails before any Kraken2 task, naming the file |
| Read length with no `kmer_distrib` | workflow | fails before any Kraken2 task |
| Barrier cardinality | workflow | `phanta_failed_samples` and `phanta_total_reads` each run exactly once |
| Column ordering | workflow | non-alphabetical samplesheet still yields samplesheet-ordered columns and sorted `total_reads.tsv` rows |
| Database manifest written | workflow | manifest records the tiny DB's `inspect.out` SHA256, tool versions, container digest |
| Wrong database rejected | workflow | with a wrong `phanta_db_expected_inspect_sha256`, fails before Kraken2 and names the mismatch |
| Manifest excluded from diff | workflow | fixture comparison passes despite the manifest's varying path and timestamp |
| `-resume` | workflow | second identical run caches all scientific tasks |
| Kraken2 module fit | module | `--report-minimizer-data` accepted; eight-column report; version >= 2.1.2 |
| Bracken module fit | module | both outputs present with upstream names |
| Each `phanta_*` module | module | explicit outputs; nothing written outside the task dir |
