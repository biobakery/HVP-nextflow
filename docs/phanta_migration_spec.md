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
prophage branch, the R merge variants, Tufts HPC, a second cloud provider,
unifying the samplesheet and `--readsdir` input models. Ordinary single-end
samples are supported through the main path.

**Targets, in order** — AWS Batch (`amazon`) first: it is where the baseline is
generated. Harvard FASRC (`harvard_rc`) second, and required: it must
reproduce the AWS baseline byte-for-byte from the same image digest. The order
is a design input, not a schedule — every "which site first" choice below
follows from it. AWS support is Phanta-only; the `amazon` profile is
documented as incomplete for every other workflow, and this work does not fix
that.

## Pinned artifacts

| Item | Value |
| --- | --- |
| HVP-nextflow checkout | `/Users/jway/source/HVP-nextflow`, branch `jw_phanta` |
| Phanta checkout (read-only) | `/Users/jway/source/phanta` |
| Phanta commit | `ecbfcf9` on `main` |
| Phanta database | Select in phase 1; recommended `unmasked_db_v1.01.tar.gz`. Pin by filename + SHA256 computed at download. 20-25 GB |
| Phanta test reads | `http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/test_dataset.tar.gz` — four `.fq.gz`, **not** in the repo |
| Kraken2 | >= 2.1.2 (needs `--report-minimizer-data`); pin 2.1.6 |
| Bracken, Python, pandas, numpy | Pin exact versions in the Phanta container; record what the baseline resolved to. **Resolved:** bracken 3.1, python 3.10.22, pandas 2.3.3, numpy 2.2.6 (also biopython 1.88, snakemake 7.15.2). Locked by URL + checksum in `containers/phanta/conda-lock.txt` |
| Phanta image | One digest in `public.ecr.aws`; every profile references it (convention 4). **Built** as `containers/phanta/` (`linux/amd64`, 1.73 GB); **not yet pushed**, so no digest exists yet |
| Tiny test database | `tests/fixtures/phanta/tiny_db/db/`, built and committed. `inspect.out` sha256 `554c323f5357f20cfd5372350676d1f55e3fc5a209015c0310c2180e7e5e299e` — the value `phanta_db_expected_inspect_sha256` takes for tiny-database runs |
| Nextflow | `24.10.4` — what `tests/run_tests.sh` pins; `manifest.nextflowVersion` is `>= 24.04.0` |
| Java | 21 |
| `manifest.version` | currently `0.0.4`; bump to `0.0.5` with the Phanta release |
| Upstream fixtures | `<phanta>/testing/{classification,final_merged_outputs}/` |

The two checkout paths are this machine's. Substitute a clone URL plus the
commits above when handing off elsewhere.

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

Five rules. The first three have bitten this repo before; the last two are
what CI enforces.

1. **Every `process` lives in `modules/**/main.nf` and nowhere else.**
   `subworkflows/` is for stages shared by more than one workflow, or a Groovy
   helper that is a top-level `def` (like `read_input.nf`). Phanta has one
   consumer, so its samplesheet parsing, per-sample chain, and cohort barriers
   go directly in `workflows/phanta.nf`; split a subworkflow out only when the
   mgx wrapper is built. Processes are lowercase snake_case, workflows
   UPPER_SNAKE; per-sample processes carry `tag "$sample"`; publish with
   `publishDir path: { "${params.outdir}/…" }, mode: 'copy'`.
2. **`meta.paired_end`, never `meta.single_end`.** Subworkflows carry
   `[ [id: sample, paired_end: bool], files ]` and unwrap to a plain sample
   string at the module boundary — see
   [subworkflows/quality_control.nf](../subworkflows/quality_control.nf) lines
   31-43. An inverted boolean here silently yields an empty channel.
3. **Phanta takes a samplesheet, not `--readsdir`.** This is the first
   samplesheet in the repo: there is no `splitCsv` anywhere, and
   `input_mapping` (mgx_mtx's TSV manifest) is the nearest precedent. Cohort
   aggregation needs an explicit ordered sample list, which filename globbing
   cannot provide.
4. **Software comes from the profile, never from a shared include.**
   [docs/dependencies.md](dependencies.md): a process runs whatever is on
   `PATH`, and what puts it there "lives in the profile … and nowhere else".
   One image digest, referenced three ways:
   - `amazon`: `container = 'public.ecr.aws/<repo>@sha256:<digest>'` in each
     Phanta `withName` block.
   - `harvard_rc`: a `def env_phanta` at the top of the profile, following the
     existing `def env_<stack>` pattern, whose `beforeScript` pulls that same
     digest with Apptainer.
   - `tufts_hpc` (deferred): `apptainer` is already enabled there; a
     `container` directive suffices.
   There is no `conf/phanta.config`; the config layering in
   [nextflow.config](../nextflow.config) has no slot for one.
5. **Every `params.phanta_*` referenced in any `.nf` has a default in
   `conf/params.config`.** Strict mode is off, so a missing default resolves
   to `null` and only fails when dereferenced — which is exactly when
   `-preview` breaks CI.

### CI obligations

CI is tool-free by policy ([tests/README.md](../tests/README.md), "What CI
checks"). It runs three things, and a Phanta change can fail each:

- **`nextflow config -profile <p>` for every profile** — any config syntax
  error, including in the new `conf/databases/aws.config`.
- **`bin/check_profile_resources.py`** — fails if a profile reuses a selector
  string that `conf/base.config` also uses and omits any of cpus/memory/time.
  So env-only blocks in `harvard_rc.config` use a selector base does not (one
  process per line, as the existing KneadData blocks do), and any real
  override restates all three directives.
- **`bin/make_diagrams.py --check`** — regenerates
  `docs/diagrams/architecture.mmd`, `docs/diagrams/<workflow>.mmd`, and
  `docs/workflow_reference.md` with `nextflow run -preview` and fails on any
  byte difference. The moment `main.nf` gains a router case,
  `architecture.mmd` changes. To pass:
  - `workflows/phanta.nf` declares `workflow PHANTA {` at column 0 on one
    line, with a single-line `//` comment directly above — it becomes the
    diagram box text.
  - `main.nf` has `case 'phanta':` with `PHANTA(` on the following line,
    exactly that shape, and `phanta` is added to the hand-maintained
    `default:` error string.
  - Each process: `^process name {` at column 0, a `//` comment above it, and
    a `publishDir` within 12 lines whose path is a double-quoted string
    containing `${params.outdir}`.
  - A `("phanta", "Phanta taxonomic profiling", [...])` tuple in `WORKFLOWS`
    in `bin/make_diagrams.py`. `make_stubs()` only writes a stub reads dir and
    a report folder — **add a `{samplesheet}` stub**, a two-row CSV pointing
    at the stub reads, and plumb it through `render_dag`. The preview must
    exit 0 under `-profile local` with `STUB_DBS` plus `--phanta_db
    /stub/phanta`. Therefore the samplesheet guard checks only that the file
    exists and parses; FASTQ and database file existence is verified inside
    tasks, which `-preview` never runs.
  - The three generated files are regenerated and committed in the same
    change.

## Target layout

```text
HVP-nextflow/
├── main.nf                                 # include PHANTA; router case; own input guard; default: string
├── nextflow.config                         # amazon profile gains conf/databases/aws.config
├── workflows/phanta.nf                     # Samplesheet parsing, per-sample chain, cohort barriers
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
│       ├── phanta_merge_tables/main.nf
│       └── version_log/main.nf             # modify: phanta_db line + params block
├── bin/
│   ├── make_diagrams.py                    # modify: WORKFLOWS tuple + {samplesheet} stub
│   ├── phanta_convert_samplesheet.py       # Legacy TSV -> CSV
│   ├── phanta_compare_outputs.sh           # diff -r wrapper
│   └── phanta_*.py|sh                      # Adapted CLIs (Appendix A.2)
├── conf/
│   ├── params.config                       # Phanta options section; phanta_db; phanta_options
│   ├── base.config                         # Phanta withName blocks
│   ├── harvard_rc.yaml                     # commented ### Phanta workflow section
│   ├── tufts_hpc.yaml                      # same, plus "not set yet" note
│   ├── profiles/
│   │   ├── aws.config                      # Phanta blocks: container digest + queue + resources
│   │   └── harvard_rc.config               # def env_phanta; env-only withName blocks
│   └── databases/
│       ├── aws.config                      # create: phanta_db + expected checksum as s3://
│       └── harvard_rc.config               # phanta_db + expected checksum
├── tests/
│   ├── run_tests.sh                        # modify: NF/PROFILE overrides, PHANTA_ONLY, Phanta cases
│   └── fixtures/phanta/{baseline,samplesheets,tiny_db,tiny_inputs}/
└── docs/
    ├── diagrams/phanta.mmd                 # generated
    ├── workflow_reference.md               # generated
    ├── dependencies.md                     # hand-maintained: env_phanta rows
    ├── phanta_migration_spec.md
    ├── phanta_migration_decisions.md
    └── phanta_testing.md                   # Baseline refresh + full-DB procedure
```

`conf/base.config` already exists with `withName` blocks and `check_max` —
extend it. `conf/params.config` is the single place a parameter is declared;
`conf/databases/` holds per-site database paths included by profile. `README.md`
and `docs/architecture.md` are hand-maintained and need a Phanta entry;
`architecture.md` is already stale (it says `version_log` writes `versions.txt`
— it writes `pipeline_log.txt`), so fix that while there.

## Process map

| Process | Upstream rule(s) | Scope | Outputs |
| --- | --- | --- | --- |
| samplesheet parsing (in `workflows/phanta.nf`) | sample-file parsing | cohort | Ordered sample metadata + ordered sample-ID channel |
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
barrier in `workflows/phanta.nf` as an explicit `collect()` then `combine` back
onto the per-sample channel. The `run_tests.sh` barrier-cardinality case
asserts each barrier process ran exactly once (one task directory per barrier
process under the case's work dir).

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
in `pipeline_info/`, which the comparison excludes — Nextflow's own
timeline/report/trace/dag files land there with launch timestamps in their
names ([conf/reports.config](../conf/reports.config)). The `run_tests.sh`
manifest case asserts the checksum fields with `grep` instead.

`modules/utils/version_log` also publishes to `pipeline_info/`. Extend it with a
`phanta_db:` line under its "Database Paths" banner and a Phanta params block,
following its indent convention; otherwise Phanta settings never appear in the
run log.

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
2. **Routine fixture regression** (the `run_tests.sh` Phanta cases) —
   HVP-nextflow only, `diff -r` against the committed baseline. Never reruns
   upstream. Not in CI: CI runs no tools.
3. **Paired-run comparison** (validating a new dataset) — upstream and
   HVP-nextflow from identical inputs, params, database, tool versions, and
   container into two directories, then `diff -r`. Do not commit the outputs.

`bin/phanta_compare_outputs.sh` wraps `diff -r`: examines the whole tree rather
than stopping at the first mismatch, lists every missing/unexpected/differing
file, excludes `.krak` and `pipeline_info/`, exits nonzero on any difference.
Investigate every mismatch; no numeric tolerances.

## Work plan

**21-30 engineer-days over 5-6 weeks.** Owners are roles until assigned.

### Phase 0 — Prerequisites: container, test database, AWS groundwork (5-7 days)

*Owner: pipeline engineer; container and AWS wiring reviewed by platform
engineer.*

Nothing downstream runs without these, and AWS is the critical path.

1. **DONE (except the push).** Build the Phanta image: Kraken2 2.1.6, Bracken,
   exactly pinned Python/pandas/numpy. Push to `public.ecr.aws` and record the
   **digest**. Every profile references that digest (convention 4).

   Built as `containers/phanta/`, `linux/amd64`, 1.73 GB. Everything installs
   from `conda-lock.txt`, a 235-package `@EXPLICIT` list pinned by URL and
   checksum — version specs alone would re-solve differently later and break
   byte-identity. Resolved: kraken2 2.1.6, bracken 3.1, python 3.10.22, pandas
   2.3.3, numpy 2.2.6, biopython 1.88, snakemake 7.15.2.

   Two notes. The Python pin is not free choice: an unconstrained solve gives
   Python 3.15 / pandas 3.0.6, and upstream's own `snakemake=7.15.2` — which the
   image needs anyway, since the baseline is generated by running the Snakefile
   *in this image* — constrains it to 3.10.22 / 2.3.3. And Bracken 3.x was
   checked against `scale_bracken.py`'s by-name column reads: `est_abundance.py`
   emits the same seven-column header as 2.x, so no downgrade is needed.

   **Still open:** the push and the digest, blocked on AWS credentials and an
   ECR repository. Until a digest exists, no profile can reference the image, so
   the `container` directives in `conf/profiles/aws.config` and the Apptainer
   pull in `def env_phanta` cannot be written against a real value.
2. **DONE.** Build a tiny Kraken2 + Bracken database per
   [Appendix D](#appendix-d-building-the-tiny-test-database). Its purpose: fast,
   deterministic test runs that exercise the Bracken-ineligible path without
   staging 25 GB. Commit it if small enough, otherwise host it and have
   `run_tests.sh` fetch it by pinned URL plus checksum.

   Committed at `tests/fixtures/phanta/tiny_db/` (recipe plus the built 5.4 MB
   database) and `tests/fixtures/phanta/tiny_inputs/` (reads plus generator),
   5.8 MB in total — see open question 3. Verified end to end by running
   upstream's own scripts against it: eight-column minimizer report, all four
   threshold classes exercised with distinct coverages, both Bracken outputs
   under their upstream names, and `sample_c` Bracken-ineligible with zero-byte
   products.

   **Appendix D needed three corrections to work at all; see the amendments
   recorded there.** Each failure was silent.
3. Obtain Phanta's test reads (separate download — record URL and checksum).
   `testing/samp_file.txt` references paths absent from a fresh checkout.
4. AWS groundwork, against the infrastructure handed over per "Inputs a fresh
   agent cannot produce":
   - Create `conf/databases/aws.config` (`phanta_db` and
     `phanta_db_expected_inspect_sha256` as `s3://` paths) and add it to the
     `amazon` block in [nextflow.config](../nextflow.config), matching how
     `harvard_rc` pairs profile and databases.
   - Add Phanta `withName` blocks to
     [conf/profiles/aws.config](../conf/profiles/aws.config) with `container`
     (the ECR digest), `queue`, and cpus/memory/time restated in full — the
     selector matches base's, so `check_profile_resources.py` requires it.
   - Replace the hardcoded `aws.batch.cliPath` with a param, or document it as
     a site input; document that `-w s3://...` must be passed.
   - Upload the tiny database to S3 and run one Kraken2 task against it.
5. Database staging decision for AWS: compare ordinary S3 task localization
   against download-once reuse on node-local NVMe, over a few Kraken2 tasks
   with the full database. Record staging time and bytes transferred; pick the
   simplest adequate option. Decided here because the phase 1 baseline runs on
   whatever is chosen.

**Deliverables:** ~~ECR image digest~~ (image built, push blocked);
~~tiny database with build script and checksum~~ **done**; recorded test-read
URL and checksum; `conf/databases/aws.config`; Phanta blocks in `aws.config`;
one Kraken2 task completed on Batch; staging decision with measurements.

**Gate:** `kraken2` and `bracken` complete against the tiny database on AWS
Batch in the pinned image, and `nextflow config -profile amazon` plus
`bin/check_profile_resources.py` pass.

**Gate status:** kraken2 and bracken complete against the tiny database in the
pinned image — but locally, not on Batch. The AWS half of the gate, and steps
3-5, remain blocked on the infrastructure handover ("Inputs a fresh agent cannot
produce", item 1).

### Phase 1 — Baseline on AWS (2-3 days)

*Owner: pipeline engineer; baseline reviewed by domain scientist.*

1. Targets are fixed: AWS first, FASRC second, Tufts deferred. No decision to
   make here — the baseline below is generated on AWS and FASRC reproduces it
   in phase 3.
2. Select and pin the database first. Choose one artifact, have the domain
   scientist confirm it suits HVP's sample types, compute SHA256 at download,
   record filename/checksum/date in `docs/phanta_testing.md`, and set
   `phanta_db_expected_inspect_sha256` per site.
3. Run upstream Phanta at commit `ecbfcf9` on the test dataset, inside the
   pinned image on an AWS launch host, against the pinned database. Record
   commands, params, and resolved tool/pandas/numpy versions. Procedure in
   [Appendix C](#appendix-c-running-upstream-phanta).
4. Run it again independently and confirm byte-identical output. If anything but
   `.krak` differs, stop and re-derive the ordering contract.
5. Commit the reviewed outputs to `tests/fixtures/phanta/baseline/` with
   checksums. Record which order each file carries.
6. Write `bin/phanta_compare_outputs.sh` with unit tests. Small script, not a
   framework.
7. Produce and commit a tiny-database baseline so the routine regression can run
   without the full database.

**Deliverables:** pinned database filename, SHA256, and selection rationale;
upstream revision and tool manifest including resolved pandas/numpy;
full-database and tiny-database baselines with checksums, both generated on
AWS; per-file ordering contract recorded; comparison wrapper with unit tests;
baseline-refresh instructions in `docs/phanta_testing.md`.

**Gate:** database pinned by checksum; two independent upstream runs
byte-identical apart from `.krak`; a clean checkout compares fixtures without
rerunning Phanta.

### Phase 2 — Port (8-10 days)

*Owner: Nextflow engineer; domain scientist reviews Phanta-specific
transformations.*

1. Add `include { PHANTA }` and a `case 'phanta':` / `PHANTA(` router case to
   [main.nf](../main.nf), in exactly the shape "CI obligations" requires.
   Leave `phanta` out of the `--readsdir` guard list; add a
   `--phanta_samplesheet` existence-and-parse guard as a `def` at the bottom
   of `main.nf`, beside `report_input_dir`. Add `phanta` to the `default:`
   error string. Implement everything in `workflows/phanta.nf`: samplesheet
   parsing emitting `[[id:,paired_end:],files]` plus an ordered sample-ID
   channel, the per-sample chain, and both cohort barriers. No subworkflow
   files, no second root pipeline.
2. Parameters, per [Appendix B](#appendix-b-parameter-mapping):
   - [conf/params.config](../conf/params.config): a `// ── Phanta options ──`
     section; `phanta_db = null` under *Database paths*; `phanta_options = ''`
     under *Per-tool extra arguments*; `phanta` added to the token comment on
     the `workflow` line.
   - [conf/harvard_rc.yaml](../conf/harvard_rc.yaml) and
     [conf/tufts_hpc.yaml](../conf/tufts_hpc.yaml): a commented
     `### Phanta workflow` section in the existing menu style, `phanta` added
     to each `### Which workflow` token list, and a line under Tufts'
     `### Databases this site has not set yet`. State there that
     `paired_end`/`single_end`/`filepattern` do not apply to Phanta.
   - `conf/databases/harvard_rc.config`: `phanta_db` and
     `phanta_db_expected_inspect_sha256` (the AWS file was created in phase 0).
3. Validate before expensive work: duplicate sample IDs, mixed library layout,
   missing read files, parameter ranges, and every file under "Required database
   files". These checks run in a task, not in the workflow body, so `-preview`
   is unaffected.
4. Implement `phanta_db_manifest` and run it as part of validation, ahead of the
   first Kraken2 task. Give it a `withName` block with generous time —
   `phanta_db_checksum_full` reads 20+ GB. Extend `version_log` as described
   under "Provenance".
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
   Either way, both contracts are asserted by the `run_tests.sh` Phanta cases
   (the report has eight columns; both Bracken outputs exist under their
   upstream names).
6. CI obligations: add the `WORKFLOWS` tuple and the `{samplesheet}` stub to
   `bin/make_diagrams.py`; add `withName` blocks to `conf/base.config` for
   every new process; add `def env_phanta` and env-only `withName` blocks to
   `conf/profiles/harvard_rc.config`; regenerate `docs/diagrams/` and
   `docs/workflow_reference.md`; add the `env_phanta` rows to
   `docs/dependencies.md`. Run all three CI checks locally before pushing.
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
Kraken2 and Bracken modules running in the pinned image; database manifest with
fail-fast checking; task-isolated scripts and modules; both cohort barriers;
deterministic aggregation; publishing contract implemented; `version_log`
extended; all three CI checks green with the regenerated diagrams committed.

**Gate:** end-to-end run on AWS against the tiny database produces all four
tables; every process declares inputs and outputs; no process scans a shared
results directory; `make_diagrams.py --check` and `check_profile_resources.py`
pass.

### Phase 3 — Validate on AWS, then reproduce on FASRC (4-6 days)

*Owner: pipeline engineer; acceptance by domain scientist; FASRC run by
platform engineer.*

1. Parameterise [tests/run_tests.sh](../tests/run_tests.sh) so one suite runs
   at both sites: `NF="${NF:-/n/lab_storage/.../nextflow}"`,
   `PROFILE="${PROFILE:--profile harvard_rc}"`, and a `PHANTA_ONLY=1` guard
   that skips the bioBakery cases (which have no tool stack on AWS). Then add
   the Phanta cases per [Appendix E.2](#e2-test-cases) using the existing
   `nf_run`/`check_run`/`check_files` helpers plus a new
   `check_identical <name> <baseline_dir>` that calls
   `bin/phanta_compare_outputs.sh`. On AWS:
   `PROFILE="-profile amazon" NF=nextflow PHANTA_ONLY=1 bash tests/run_tests.sh`.
2. Run the Phanta cases on AWS against the tiny database.
3. Verify `-resume` on AWS: run the same case twice, confirm scientific tasks
   are cached.
4. Run the full-database workflow on AWS. Require byte-identical output against
   the full-database baseline.
5. **FASRC parity.** Resolve open question 2 first. Pull the ECR digest with
   Apptainer via `def env_phanta`; stage the database on shared storage with
   node-local scratch; add `phanta_db` and the expected checksum to
   `conf/databases/harvard_rc.config`; run the same Phanta cases with the
   default `PROFILE`; require byte-identity against the **AWS-generated**
   baselines, tiny and full.
6. Paired-run comparison on at least one representative real cohort. Retain the
   report; obtain domain-scientist signoff.
7. Keep the Phanta checkout runnable as the comparator through this milestone.

**Deliverables:** parameterised `run_tests.sh` with passing Phanta cases;
tiny-database regression report; AWS resume evidence; full-database comparison
report on AWS; FASRC parity report against the AWS baselines; real-cohort
paired-run report and domain-scientist signoff.

**Gate:** tiny-database regression, full-database runs, and the real-cohort
comparison all byte-identical apart from excluded files — on AWS, and then on
FASRC against the same baselines.

### Phase 4 — Tune (1-2 days)

*Owner: platform engineer with pipeline engineer.*

1. Targeted refinements to the AWS staging approach chosen in phase 0 only.
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

**Deliverables:** resource defaults in `conf/base.config`; staging refinements
with measurements; interruption/resume result; merge memory measurement; README
section, `docs/architecture.md` entry, and user runbook; `manifest.version`
bumped to `0.0.5`.

**Gate:** AWS and FASRC both have documented measured defaults; AWS recovers
from one forced interruption.

## Acceptance criteria

1. No migration files, tracked or untracked, in the Phanta repository.
2. One database artifact selected, pinned by filename and SHA256, confirmed
   suitable by the domain scientist, used for both baseline and HVP runs.
3. Every published run records its database: `phanta_db_manifest.txt` with the
   `inspect.out` checksum, tool versions, and container digest. A run against a
   database whose checksum does not match
   `phanta_db_expected_inspect_sha256` fails before any Kraken2 task.
4. Two independent upstream baseline runs on AWS byte-identical apart from
   `.krak`; routine regression uses committed fixtures without rerunning
   Phanta.
5. Tiny database and ECR image exist; the fixture regression runs as
   `run_tests.sh` cases on AWS and on FASRC without the 25 GB database.
6. All Phanta processes run from the same image digest on `amazon` and
   `harvard_rc`, each through that profile's own mechanism.
7. The `run_tests.sh` Phanta cases assert Kraken2's eight-column minimizer
   report and both Bracken outputs under their upstream names.
8. Both cohort barriers implemented explicitly; the barrier-cardinality case
   passes.
9. No module writes outside its task directory; all outputs and
   biological-empty behavior explicit.
10. All four tables byte-identical to the AWS-generated baseline against both
    tiny and full databases, on AWS and then on FASRC.
11. A representative real cohort passes comparison with domain-scientist
    signoff.
12. `-resume` passes on AWS and on FASRC, including one AWS Spot interruption.
13. All three CI checks pass with the regenerated `docs/diagrams/` and
    `docs/workflow_reference.md` committed.
14. Database staging, Kraken2 sizing, and one aggregation-memory measurement
    documented without speculative scale tiers.

## Decisions already made

Do not relitigate these during implementation. Rationale for each is in
[phanta_migration_decisions.md](phanta_migration_decisions.md).

- All implementation lives in HVP-nextflow; Phanta stays read-only.
- AWS Batch first and generates the baseline; Harvard FASRC second and must
  reproduce it; Tufts deferred. AWS support is Phanta-only.
- Integrate through the existing router; every process in `modules/**/main.nf`;
  the whole workflow in `workflows/phanta.nf`, no subworkflow files until a
  second consumer exists; no separate nf-core repository.
- One image digest in ECR, referenced per profile through that profile's own
  mechanism. No shared `conf/phanta.config`.
- Tests are `run_tests.sh` cases, parameterised by site. CI stays tool-free;
  Phanta's CI contribution is the three existing checks passing.
- Prefer pinned nf-core Kraken2 and Bracken modules with a one-day fit limit and
  local fallbacks; the per-profile container overrides any module's own
  container directive.
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

## Inputs a fresh agent cannot produce

Hand these over before starting, so the implementer neither guesses nor stalls.
Each is tagged with the phase it blocks.

1. **AWS infrastructure** (blocks phase 0): a Batch compute environment with
   NVMe-backed instance types for Kraken2 (≥ 64 GB RAM to hold the 25 GB
   index), job queue names, an S3 bucket for the database and the work dir, an
   ECR repository to push to, and the IAM role the launch host assumes.
2. **The database decision** (blocks phase 1, not phase 0): open question 1 —
   or permission to proceed with `unmasked_db_v1.01.tar.gz` pending signoff.
   Phase 0 and the tiny-DB work need no answer.
3. **FASRC container policy** (blocks phase 3 step 5): open question 2, plus a
   shared-storage path for the 25 GB database.
4. **The acceptance cohort** (blocks phase 3 step 6): open question 4.
5. **Repo access** if not on this machine: a clone URL plus the commits in
   "Pinned artifacts".

The tiny-database genomes are *not* an input — Appendix D pins them.

## Open questions

Resolve with the named owner; record the answer here.

1. **Which database, and is it suitable?** `README.md` and `databases.md` name
   different artifacts (`database_V1.tar.gz` vs `unmasked_db_v1.01.tar.gz`) and
   upstream publishes no checksums. The latter is recommended. The open part is
   scientific: Phanta's default database is human-gut specific, so confirm it
   fits HVP's sample types. *Domain scientist confirms; engineer records
   filename, SHA256, date. Gates phase 1.*
2. **Does FASRC permit Apptainer for this stack?** If yes, `def env_phanta`
   pulls the ECR digest and byte-identity with AWS holds. If FASRC requires a
   hutlab module instead, byte-identity across sites is **not** guaranteed and
   FASRC gets its own baseline — record that explicitly. *Platform engineer.
   Gates phase 3 step 5, not AWS work.*
3. ~~**Tiny database committed or downloaded?** Depends on achieved size.~~
   **ANSWERED — committed.** The built database is 5.4 MB and the synthetic
   reads 0.4 MB, 5.8 MB in total, so hosting plus a checksummed fetch in
   `run_tests.sh` is not worth its weight. Committed at
   `tests/fixtures/phanta/tiny_db/db/` alongside the recipe that produced it, so
   a test run needs neither Docker nor network access. Only `hash.k2d` (5.6 MB)
   is binary-large; the other eight files are small text that will diff
   readably. *Pipeline engineer, phase 0.*
4. **Which real cohort is the acceptance dataset?** *Domain scientist. Gates
   phase 3 step 6.*
5. **Emit `pipeline_completed.txt` or exclude it?** *Pipeline engineer, phase 2.*

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
| — | `phanta_options` | `''` — appended verbatim to the Kraken2 command, per the repo's `<tool>_options` convention |
| — | `phanta_db_expected_inspect_sha256` | `null`; per site. Mismatch fails the run |
| — | `phanta_db_checksum_full` | `false`; `true` for baselines and release validation |

There is no `run_phanta`: `run_*` are stage toggles inside mgx/mtx, and Phanta
is selected with `--workflow phanta`. Library layout comes from the
samplesheet, so Phanta ignores the existing
`paired_end`/`single_end`/`filepattern` params; say so in the
`### Phanta workflow` section of both `conf/*.yaml` templates.

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

# Appendix D: building the tiny test database

**BUILT.** Implemented as `tests/fixtures/phanta/tiny_db/build.sh`; the built
database, its recipe, and the synthetic reads are committed. Achieved 5.4 MB,
well under the target.

**Three amendments, recorded because following the text below as originally
written produces a database that is broken in ways nothing reports.** The
original wording is kept; the amendments are marked inline.

1. *Step 2 — do not download the taxonomy.* `filter_kraken_reports.py` looks up
   the literal rank string `superkingdom`, and NCBI has since renamed that rank
   to `domain` (and `Viruses` to `acellular root`), which was confirmed against
   live NCBI. With a freshly downloaded taxonomy every superkingdom lookup
   fails, falls through the catch-all `else` at `filter_kraken_reports.py:248`
   that **keeps** the species, and filtering silently does nothing — while the
   output still looks entirely normal. A download would also make the database
   differ run to run, since NCBI's dump is a moving target. The committed
   `taxonomy/nodes.dmp` and `names.dmp` are hand-authored and use
   `superkingdom`, matching what the real Phanta database ships. This does
   **not** affect the full database: the taxonomy is read from inside the
   database directory, and each Phanta tarball carries its own frozen snapshot.
2. *Step 1 — four genomes are not enough; there is a fifth.* With only genomes
   that resolve cleanly to a species, every row of
   `library/species_genome_size.txt` is numeric, pandas types `species_level_taxa`
   as `int64`, and `scale_bracken.py` dies on `ValueError: You are trying to
   merge on object and int64 columns` because it casts the Bracken side to
   `str`. The real database never hits this — it holds many genomes classified
   no further than genus, for which `calculate_genome_size.py` writes the literal
   `error - taxid above species`, forcing the column to object and the merge to
   string-to-string. This is **not** a pandas-version issue: pandas 1.5.3 raises
   the identical error, which was checked directly. A fifth genome
   (`NC_000912.1`, tagged at genus `2995234`) makes the fixture representative.
3. *Step 5 — NCBI dump line endings are load-bearing.* Dump lines end `\t|`, not
   `\t|\t`. Kraken2 discards the final two characters before reading the last
   field, so a line ending `\t|\t` leaves a stray `\t|` glued to the name class,
   the `scientific name` test never matches, and **every name in `inspect.out`
   comes out empty** — silently. Phanta reads taxon names from `inspect.out`, so
   the merged tables would carry blank lineage names. Add this to step 5's
   presence checks: verify `inspect.out` has a populated name column, not just
   that the file exists.

Accessions are pinned by *versioned* accession and SHA256 in
`tests/fixtures/phanta/tiny_db/accessions.tsv`; the table below gives
unversioned ones. Full detail is in
[tests/fixtures/phanta/tiny_db/README.md](../tests/fixtures/phanta/tiny_db/README.md).

No upstream recipe exists. Target under ~500 MB; the genomes below total about
2.5 Mb of sequence, so the Kraken2 index will be tens of MB.

1. Use these four genomes — one per threshold class, so every
   `cov_thresh_*`/`minimizer_thresh_*` branch in `filter_kraken_reports.py` is
   exercised. Fetch by RefSeq accession; do not substitute without updating the
   baseline.

   | Class | Organism | Accession | Size |
   | --- | --- | --- | --- |
   | viral | phiX174 | `NC_001422` | 5 kb |
   | bacterial | *Mycoplasma genitalium* G37 | `NC_000908` | 580 kb |
   | archaeal | *Methanocaldococcus jannaschii* | `NC_000909` | 1.7 Mb |
   | eukaryotic | *S. cerevisiae* chromosome I | `NC_001133` | 230 kb |
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

If the result exceeds what should live in git, host it (S3 for AWS, shared
storage for FASRC) and have `run_tests.sh` fetch it by pinned URL plus checksum.
It did not — see open question 3. Only the nine files in "Required database
files" are committed; the `library/` source FASTAs are dropped, which also
mirrors how upstream ships the real database.

---

# Appendix E: deliverables and test matrix

## E.1 Files to create or modify

| Path | Action | Done when |
| --- | --- | --- |
| `main.nf` | modify | `include { PHANTA }`; `case 'phanta':` / `PHANTA(` in the exact shape; excluded from the `--readsdir` guard; `--phanta_samplesheet` guard as a `def`; `phanta` in the `default:` string |
| `nextflow.config` | modify | `amazon` profile also includes `conf/databases/aws.config` |
| `conf/params.config` | modify | `// ── Phanta options ──` section with every Appendix B default; `phanta_db = null`; `phanta_options = ''`; `phanta` in the workflow-token comment |
| `conf/harvard_rc.yaml`, `conf/tufts_hpc.yaml` | modify | commented `### Phanta workflow` section; `phanta` in `### Which workflow`; Tufts "not set yet" line |
| `conf/base.config` | modify | `withName` blocks for `kraken2`, `bracken`, `phanta_db_manifest`, and the small table steps |
| `conf/profiles/aws.config` | modify | Phanta blocks with `container` digest, `queue`, cpus/memory/time restated; `cliPath` no longer hardcoded |
| `conf/profiles/harvard_rc.config` | modify | `def env_phanta` (Apptainer pull of the digest); env-only `withName` blocks on selectors base does not use |
| `conf/databases/aws.config` | create | `phanta_db` + `phanta_db_expected_inspect_sha256` as `s3://` paths |
| `conf/databases/harvard_rc.config` | modify | `phanta_db` + `phanta_db_expected_inspect_sha256` |
| `workflows/phanta.nf` | create | `workflow PHANTA {` at col 0 with a `//` line above; samplesheet parsing emitting `[[id:,paired_end:],files]` + ordered sample-ID channel; per-sample chain; both barriers explicit; publishes |
| `modules/kraken2/main.nf` | create | A.1 command; minimizer report; `.krak` unpublished |
| `modules/bracken/main.nf` | create | A.1 command; input staged as `.krak.report.filtered`; both outputs emitted |
| `modules/taxonomy/phanta_filter/main.nf` | create | filtered report, species table, decisions, eligibility status |
| `modules/taxonomy/phanta_scale/main.nf` | create | scaled table; empty passthrough when ineligible |
| `modules/utils/phanta_db_manifest/main.nf` | create | writes the manifest; fails on checksum mismatch; runs before Kraken2 |
| `modules/utils/phanta_failed_samples/main.nf` | create | sorted failed list |
| `modules/utils/phanta_prepare_counts/main.nf` | create | `.to_merge` plus the per-sample total-reads row |
| `modules/utils/phanta_normalize_counts/main.nf` | create | one sample per task |
| `modules/utils/phanta_merge_tables/main.nf` | create | three tables, columns in samplesheet order |
| `modules/utils/version_log/main.nf` | modify | `phanta_db:` under "Database Paths"; Phanta params block |
| `bin/make_diagrams.py` | modify | `phanta` tuple in `WORKFLOWS`; `{samplesheet}` stub in `make_stubs()` plumbed through `render_dag` |
| `bin/phanta_convert_samplesheet.py` | create | legacy TSV → CSV |
| `bin/phanta_*.py`, `bin/phanta_calc_total_reads.sh` | create | explicit I/O paths; no derived output names |
| `bin/phanta_compare_outputs.sh` | create | `diff -r` wrapper; excludes `.krak` and `pipeline_info/`; nonzero on any diff |
| `tests/run_tests.sh` | modify | `NF`/`PROFILE` env overrides; `PHANTA_ONLY` guard; `check_identical` helper; E.2 cases |
| `tests/fixtures/phanta/**` | create | tiny DB (or fetch script + checksum), inputs, samplesheets, baselines — **tiny DB and inputs done** (committed, with build script and pinned accessions); samplesheets and baselines still to come |
| `docs/diagrams/architecture.mmd`, `docs/diagrams/phanta.mmd`, `docs/workflow_reference.md` | generate | `bin/make_diagrams.py` output committed; `--check` passes |
| `docs/dependencies.md` | modify | `env_phanta` row in the stack table; `amazon` row updated |
| `docs/architecture.md` | modify | Phanta entry; stale `version_log` claim fixed |
| `docs/phanta_testing.md` | create | baseline refresh, full-DB procedure, resolved versions, checksums |
| `README.md` | modify | Phanta in Workflow & I/O, Project Structure, Project Status |

## E.2 Test cases

All are `tests/run_tests.sh` cases under `PHANTA_ONLY=1`, run first on AWS
(`PROFILE="-profile amazon"`) and then on FASRC (default). Numbering continues
from the existing 16.

| Case | Asserts |
| --- | --- |
| `test17_phanta_pe` — paired-end cohort, tiny DB | `check_run`; `check_identical` against the tiny-DB baseline: four tables byte-identical |
| `test18_phanta_se` — ordinary single-end cohort | same, single-end path |
| `test19_phanta_ineligible` — one Bracken-ineligible sample | `check_run`; zero-byte products present; sample listed in `samples_that_failed_bracken.txt` |
| `test20_phanta_ordering` — non-alphabetical samplesheet | merged-table header columns in samplesheet order; `total_reads.tsv` rows sorted |
| `test21_phanta_manifest` | `pipeline_info/phanta_db_manifest.txt` exists; `grep` finds the tiny DB's `inspect.out` SHA256, Kraken2 version, and the image digest |
| `test22_phanta_wrong_db` | with a wrong `phanta_db_expected_inspect_sha256`: nonzero exit; log names the mismatch; no `kraken2` task directory exists under the case's work dir |
| `test23_phanta_barriers` | exactly one task directory each for `phanta_failed_samples` and `phanta_total_reads` under the case's work dir |
| `test24_phanta_resume` | second identical `nf_run` with `-resume`: all scientific tasks reported as cached in the log |
| `test25_phanta_guards` — fast, no cluster | duplicate sample IDs, mixed layout, missing database file, and read length with no `kmer_distrib` each fail before any Kraken2 task with a clear message (same shape as the existing tests 15-16) |
| Tool contracts (asserted inside `test17`) | Kraken2 report has eight columns; both Bracken outputs exist under their upstream names; no file written outside any task dir |
