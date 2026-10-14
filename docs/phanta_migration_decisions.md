# Phanta Migration: Decisions and Rationale

Why the migration is specified the way it is. The instructions themselves are in
[phanta_migration_spec.md](phanta_migration_spec.md) — this document exists for
whoever asks, months from now, why a particular choice was made, and to record
what the planning got wrong along the way.

Each section states a decision, the evidence behind it, and what was rejected.

---

## 1. No canonicalization layer; ordering is fixed at emit time

**Decision:** HVP-nextflow forces output order per file as it writes, and
comparison is a `diff -r` wrapper. There is no canonicalization utility and no
bespoke comparison tool.

An earlier draft assumed upstream Phanta suffered from glob-order
nondeterminism, and built a substantial apparatus around it: a
`phanta_canonicalize_outputs.py`, a `phanta_compare_outputs.py` with cell-level
diff artifacts, bounded previews and binary checksums, and a baseline-refresh
wrapper chaining them. Reading the code does not support the premise.

**Upstream is deterministic for a fixed sample file and fixed inputs:**

- `merge.py:35` and `merging_scaled_bracken.py:33` use
  `concurrent.futures.Executor.map`, which yields results **in input order**
  regardless of thread count. Sample columns come out in `expand(...)` order,
  i.e. sample-file order.
- `rule tot_reads_all` builds `total_reads.tsv` with
  `cat {classdir}/*total_reads.txt`; `rule samples_failed_bracken` uses
  `ls {classdir}/*.temp`. Shell globs are **sorted** — deterministic, just
  alphabetical rather than sample-file order.
- The only thread-order-dependent artifact is the raw `.krak` per-read
  assignment file, excluded from equivalence and from publishing.

**A single canonicalization rule would have been actively wrong.** Two
different orders coexist in one upstream run. From Phanta's committed fixtures,
generated from a sample file listing `sample19` before `sample18`:

```
total_reads.tsv     rows:    sample18, sample19    (alphabetical, from the glob)
counts.txt          columns: sample19, sample18    (sample-file order)
```

"Canonicalize everything to samplesheet order" would reorder `total_reads.tsv`
and destroy the match. The non-alphabetical upstream sample file is what exposes
this, which is why the spec says to keep it when reproducing the baseline.

**The real risk is on the HVP side.** Nextflow channel arrival order *is*
nondeterministic. So the problem is real, but it belongs to the emitting
process, not to a post-hoc cleanup pass. Preventing disorder is cheaper than
absorbing it, and it removes two scripts and a framework from the deliverables.

**Rejected:** a canonicalizer (wrong premise, and its single rule would corrupt
`total_reads.tsv`); a custom comparison tool with diff artifacts and bounded
previews (over-built for comparing four TSVs — `git diff --no-index` gives
readable tabular diffs).

---

## 2. One pinned container is a prerequisite, not deferred work

**Decision:** build one Phanta image before baselining, push it to ECR, and
reference that single digest from every profile through that profile's own
mechanism — `container` in `amazon`'s `withName` blocks, an Apptainer pull in a
`def env_phanta` `beforeScript` on `harvard_rc`.

Byte-identity is a property of the software environment, not of the port.
`filter_kraken_reports.py`, `scale_bracken.py`, `merge.py`, and
`merging_scaled_bracken.py` all write through `pandas.to_csv`, whose float
formatting is version-sensitive. The evidence is visible in the fixtures:
`counts.txt` values appear as `6564345.0`, floats produced by
`pd.concat(...).fillna(0)` then `to_csv`.

Upstream does not pin these. `phanta_env.yaml` pins `kraken2=2.1.6` and
`snakemake=7.15.2` but leaves `bracken`, `pandas`, and `r-base` unpinned. So the
baseline's exact pandas and numpy versions must be captured at generation time
and frozen in the image.

This repository delivers software per profile, and says so as a rule —
[docs/dependencies.md](dependencies.md): a process runs whatever is on `PATH`,
and what puts it there "lives in the profile … and nowhere else". `amazon` uses
per-process ECR containers, `tufts_hpc` has Apptainer enabled, and `harvard_rc`
loads hutlab environment modules through a `def env_<stack>` per software stack.
An earlier draft proposed a shared `conf/phanta.config` carrying one container
directive for all profiles; the config layering in `nextflow.config` has no slot
for such a file, and it would have been the only place in the repo where a tool
environment lived outside a profile.

So the mechanism follows the repo, and the *digest* is what is shared. Because
AWS is the first target, the ECR image is canonical and the baseline is
generated there. FASRC then reproduces it by pulling the same digest with
Apptainer — which `tufts_hpc` already does for its own stack, so there is
precedent.

**The caveat that matters:** byte-identity across sites holds only if FASRC
runs the container. If FASRC policy requires a hutlab module for this stack
instead, pandas/numpy versions are no longer guaranteed to match the image, and
FASRC must generate its own baseline. That is open question 2 in the spec; it
gates FASRC work, not AWS work.

**Rejected:** a shared `conf/phanta.config` (layering violation, and only
"shared" in name — each executor still needs its own delivery mechanism);
deferring containers to a separate "custom container / SBOM program" (as an
earlier draft did) while simultaneously requiring byte-identical output. Those
two positions are incompatible.

---

## 3. Keep `calc_total_reads.sh` in shell for v1

**Decision:** wrap it with explicit paths; do not rewrite it in Python yet.

An earlier draft called for porting it to Python to remove the `bc` dependency
and the directory glob. Removing the glob is necessary for task isolation; the
language change is not, and it carries risk against a byte-identity bar: bash
integer arithmetic and Python can differ in output formatting. The script also
has a behavior that is easy to lose in translation — `exit 64` with "No reads in
sample" when both the Kraken root and unclassified counts are absent.

Rewriting working code before a baseline exists inverts the order of operations.
Revisit after the comparison passes.

---

## 4. A samplesheet, accepted as a second input model

**Decision:** Phanta takes a `sample,fastq_1,fastq_2` CSV. The other workflows
keep `--readsdir`.

An earlier draft claimed this "follows existing conventions". It does not. This
repository has exactly one input model: `--readsdir` plus filepattern globbing
with sample IDs derived from filenames
([subworkflows/read_input.nf](../subworkflows/read_input.nf)), and
[main.nf](../main.nf) hard-errors when `--readsdir` is unset for read-based
workflows. A CSV is a genuinely new second model — and with the legacy TSV
converter, a third user-facing format.

It is still the right call. Cohort aggregation needs an explicit, ordered,
user-supplied sample list; filename globbing cannot provide one, and the merged
tables' column order depends on having it. But the cost is real and should be
documented rather than disguised: two input conventions in one repository is a
papercut for users, and the earlier draft's own instruction to "not maintain two
internal channel models" was being violated at the repository level while being
honored at the channel level.

Mitigations in the spec: one internal channel model
(`[[id:,paired_end:],files]`), a converter for the legacy format, and explicit
documentation of the difference. Unifying the two input models is deferred, not
denied.

**Reconsidered after merging `main` (v0.0.4), and kept.** The merge removed
what little precedent there was: `conf/params.config` now declares itself "the
single place a parameter is declared", there is no `splitCsv` anywhere in the
tree, and `bin/make_diagrams.py` can only stub a reads directory — so a
samplesheet needs its own stub for the CI graph check. The alternative was put
to the user plainly: make Phanta an mgx stage toggle (`--run_phanta true`),
consuming the existing `read_input.nf` channel after KneadData, with no CSV at
all. Its costs: sample column order would have to become alphabetical (and the
upstream baseline regenerated to match), every Phanta test would need
KneadData's host database, and the relationship between `run_phanta` and
`run_taxonomic_profiling` (MetaPhlAn) would need defining. The user chose the
standalone workflow. A thin mgx wrapper that feeds KneadData output into the
same chain remains a possible later addition — which is also why the chain
lives in `workflows/phanta.nf` rather than a subworkflow until that second
consumer exists.

---

## 5. Fold upstream rules rather than one process per rule

**Decision:** `phanta_filter` covers `filter_kraken` + `process_filtered_kraken`;
`phanta_prepare_counts` covers `prepare_to_merge_counts` +
`tot_reads_per_sample`. Eight processes, not ten.

Six of the ten upstream per-sample rules are small text transformations. One
Nextflow process per rule means six extra container starts plus stage/unstage
rounds per sample — material on AWS Batch for a large cohort, and noise in the
trace on Slurm. Each folded pair reads the same inputs at the same point in the
dataflow with no barrier between them, so folding costs nothing in
parallelism.

Resource blocks follow the existing "small table and folder assembly steps"
pattern already in [conf/base.config](../conf/base.config): 1 cpu, 8 GB, 1 h.

**Rejected:** one-process-per-rule symmetry. It reads tidier against the
`Snakefile` and costs real money.

---

## 6. Provenance capture, because upstream has none

**Decision:** every run writes `pipeline_info/phanta_db_manifest.txt`, and a
checksum mismatch against `phanta_db_expected_inspect_sha256` fails the run
before any Kraken2 task.

Upstream leaves no trace of which database produced a result. The database
arrives as a bare filesystem path — `config['database']`, consumed at eight
sites in the `Snakefile` — filled in by the operator from a placeholder in
`config.yaml`. No version string, no checksum, no manifest is written anywhere.
`run_phanta.py`'s `getDefaults()` compounds it by validating a candidate
database directory through the presence of `species_name_to_vir_score.txt`, a
*post-processing* file the core pipeline never reads — so a directory can pass
that check while missing `inspect.out` or `library/species_genome_size.txt`.

The consequence: given a completed Phanta output directory, there is no way to
determine which database build produced it. Results cannot be audited or
reproduced without out-of-band notes about what someone typed into a config
file.

Two design points worth preserving:

**The manifest must be a sidecar.** Embedding provenance as a header, comment,
or column in the four tables would break byte-identity with the baseline
immediately, since provenance varies per run. This is the one place where the
obvious implementation is wrong.

**`inspect.out` is the fingerprint, not `hash.k2d`.** `inspect.out` is a text
dump of every taxon with clade k-mer counts — content-derived, so two different
builds cannot share its checksum, and small enough to hash every run.
`hash.k2d` is 20+ GB, so a full SHA256 of it costs a complete read and is gated
behind `phanta_db_checksum_full`, required only for baselines and release
validation.

The manifest carries a path and timestamp, so it is not byte-stable across runs.
It lives in `pipeline_info/`, which the comparison excludes — Nextflow's own
timeline/report/trace/dag files there are already timestamped — and its checksum
fields are asserted by a `run_tests.sh` case with `grep` instead.

---

## 7. The database is a fresh choice, pinned by checksum

**Decision:** select one artifact in phase 1 (recommended
`unmasked_db_v1.01.tar.gz`), pin it by filename and SHA256 computed at download,
and use it for both the upstream baseline and all HVP runs.

Upstream contradicts itself about which database is "the default":

| Upstream source | Artifact |
| --- | --- |
| `README.md` line 76, the install step everyone follows | `database_V1.tar.gz` |
| `databases.md` line 5, the row labelled "Default database" | `unmasked_db_v1.01.tar.gz` |

`V1` versus `v1.01`, both presented as the default. Not a clean alias either:
the optional genome-FASTA add-on is `database_V1.library.tar.gz`, following the
README's naming lineage rather than the table's. Two people following two
different upstream docs can profile against different databases and get
different counts for identical reads.

`unmasked_db_v1.01.tar.gz` is recommended because `databases.md` is the
maintained catalogue — it carries recent commit activity, including at the
pinned commit `ecbfcf9`, while the README's link has not been updated alongside
it.

A URL is not a version pin. Upstream publishes no checksums, so the tarball can
be replaced in place; hence pinning by SHA256 computed at download, and treating
the catalogue as a moving target to re-verify at each baseline refresh.

**Determining whether the two tarballs are the same build is optional.** It
costs a second 25 GB transfer and changes nothing, because we control which
database both sides of the comparison use.

The remaining question is scientific, not technical: Phanta's default database
was built from human-gut viral and bacterial genomes, so its suitability for
HVP's sample types needs domain-scientist confirmation. Upstream's
`databases.md` says as much and points non-gut users toward custom databases.

### Full database catalogue

Recorded so a future switch is a documented decision rather than a
rediscovery. All at `http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/` unless
noted; each 20-25 GB.

| Database | Artifact | Prokaryotic | Viral | Taxonomy |
| --- | --- | --- | --- | --- |
| Default (per `databases.md`) | `unmasked_db_v1.01.tar.gz` | HumGut | MGV + RefSeq viral | NCBI |
| Default (per `README.md`) | `database_V1.tar.gz` | HumGut | MGV + RefSeq viral | NCBI |
| Prophage-masked default | `masked_db_v1.tar.gz` | HumGut | MGV + RefSeq viral | NCBI |
| Default, GTDB taxonomy | `unmasked_db_v1_gtdb.tar.gz` | HumGut | MGV + RefSeq viral | GTDB |
| UHGGv2 + MGV | `uhgg2_mgv_v1.tar.gz` | UHGGv2 | MGV + RefSeq viral | GTDB |
| HumGut + UHGV "MQ+" | `humgut_uhgv_mqplus_v1.tar.gz` | HumGut | UHGV >= medium quality | NCBI |
| HumGut + UHGV "HQ+" | `humgut_uhgv_hqplus_v1.tar.gz` | HumGut | UHGV >= high quality | NCBI |
| UHGGv2 + UHGV "MQ+" | `uhgg2_uhgv_v2.tar.gz` | UHGGv2 | UHGV >= medium quality | GTDB |
| MAGIC (early-life gut) | Zenodo `10.5281/zenodo.10369093` | MAGIC_pMAGs | MAGIC_vMAGs | GTDB |

Only the masked database is referenced by upstream's own docs as a deliberate
alternative; the rest are community or extension builds.

---

## 8. A tiny test database is mandatory, and had to be invented

**Decision:** build a small Kraken2 + Bracken database as a phase 0
deliverable, with a committed build script and pinned genome accessions.

An earlier draft required a "routine fixture regression" but never scheduled
the artifact that makes it possible. Phanta's committed `testing/` directory
holds expected *outputs* but no database and no reads — the test reads are a
separate `test_dataset.tar.gz` download, and `testing/samp_file.txt` references
paths absent from a fresh checkout.

The original justification was CI: the full database is 20-25 GB and a
GitHub-hosted runner has about 14 GB of disk. That justification evaporated
with the merge of `main` — the repo's CI is tool-free by policy and runs no
pipeline at all (see §10). The tiny database still earns its place on three
grounds: a regression run completes in minutes rather than hours; it needs no
25 GB staged per test run on either site; and it is the only practical way to
include a deliberately Bracken-ineligible sample, which exercises the zero-byte
substitute outputs that the real test reads never trigger.

No upstream recipe exists for a minimal database, so the spec carries one
(Appendix D), with four named RefSeq accessions spanning the four threshold
classes — bacterial, viral, archaeal, eukaryotic — so every filtering branch
runs. The accessions are pinned so a fresh implementer does not choose, and so
a baseline generated from them stays valid.

This is also why the effort estimate moved from 10-15 to 16-22 engineer-days:
the container and the tiny database are several days of work that the earlier
plan assumed away. (It moved again, to 21-30, for the AWS groundwork and FASRC
parity — see §9.)

---

## 9. AWS first and generates the baseline; FASRC second and reproduces it

**Decision:** AWS Batch (`amazon`) is the first target and the site where the
baseline is generated. Harvard FASRC (`harvard_rc`) is the second, required
target: it must reproduce the AWS baseline byte-for-byte. Tufts is deferred. A
second cloud provider stays deferred.

An earlier draft opened with "select the initial cloud on day one … implement
one of AWS Batch/S3 or Google Batch/GCS", and a later one treated the Slurm
profiles as primary with AWS as an optional add-on. The user settled it: AWS
first, FASRC second, both required.

The `amazon` profile is explicitly incomplete. Its header, post-merge: "AWS
Batch. INCOMPLETE, and unmaintained since v0.0.2 … configures containers and
queues for the six read-based core processes only … Treat this as a starting
point for an AWS port, not as a working profile." It has no
`conf/databases/aws.config`, a hardcoded `aws.batch.cliPath`, and requires `-w
s3://...` on every launch. Phanta is therefore the first workflow to make
`amazon` actually work. That scope is deliberately narrow: the profile stays
broken for the other fifty-odd processes, and fixing them is not Phanta's job.

**Why the order is a design input, not a schedule.** Three things follow from
"AWS first":

- *Baseline provenance.* The byte-identity bar needs one canonical environment.
  AWS is container-native, so the ECR digest is the canonical environment with
  no translation step, and the baseline carries that digest in its manifest.
  FASRC's job becomes "match this", which is a cleaner acceptance test than two
  sites each generating a baseline and hoping they agree.
- *Where the groundwork goes.* `conf/databases/aws.config`, the Phanta blocks in
  `aws.config`, and the S3-vs-NVMe staging decision all moved into phase 0 as
  the critical path, rather than sitting in phase 2 behind the port.
- *What the test suite must do.* `tests/run_tests.sh` is bound to FASRC — it
  hardcodes the Nextflow path, `-profile harvard_rc`, and FASRC fixture
  locations. Phanta's tests have to run on AWS first, so the suite is
  parameterised (§12) rather than waiting for FASRC.

**Rejected:** FASRC-first (would generate the baseline under a module stack
whose pandas/numpy are not pinned to the image — see §2); AWS-only (the HPC
cluster is where the lab actually runs); dropping AWS from v1 (the user's call
was to keep it).

---

## 10. Tests are `run_tests.sh` cases; CI stays tool-free

**Decision:** Phanta's tests are cases in `tests/run_tests.sh`, the repo's one
integration suite. CI does not run them. Phanta's CI obligation is that the
three existing tool-free checks pass.

Earlier drafts planned an nf-test suite under `tests/` and a "routine fixture
regression" running in GitHub Actions. Both targets were removed by the merge of
`main`: the v0.0.4 cleanup deleted nf-test entirely because its cases "asserted
a four-process pipeline, loaded params files that were not in the repo, and so
had not run in a long time" (`tests/run_tests.sh` header), and CI was rewritten
to run only what needs no tools — profile parsing, `check_profile_resources.py`,
and `make_diagrams.py --check` (`tests/README.md`, "What CI checks"). Those are
stated policy, not accident.

So the tests follow the suite that exists: bash cases using `nf_run`,
`check_run`, `check_files`, parallel drivers, `.status` files, plus one new
helper for byte-identity. The assertions from the earlier nf-test matrix are all
preserved; only the harness changed.

**Rejected:** a Phanta-only tool-running CI job. The tiny database and the
container would make it feasible on a GitHub runner, so it is not a technical
impossibility — it contradicts the repo's stated CI policy. If that policy
changes, Phanta is the easiest workflow to bring under it. Also rejected:
re-introducing nf-test for Phanta alone (two test frameworks in one repo, one of
them just deleted on purpose).

---

## 11. The generated docs are part of the deliverable

**Decision:** `docs/diagrams/architecture.mmd`, `docs/diagrams/phanta.mmd`, and
`docs/workflow_reference.md` are regenerated with `bin/make_diagrams.py` and
committed in the same change as the workflow.

This is not optional tidiness. `make_diagrams.py --check` runs in CI and fails
on any byte difference between the committed files and a fresh regeneration.
The architecture diagram is parsed from `main.nf`'s router and the `include`
statements, so `architecture.mmd` changes the moment `case 'phanta':` is added
— before any `WORKFLOWS` entry exists. A Phanta pull request that does not
regenerate the diagrams is red on arrival.

Two consequences the spec now carries. The script's regexes impose formatting
on the code: `workflow PHANTA {` at column 0 with a `//` comment directly
above, `case 'phanta':` with `PHANTA(` on the next line, every process with a
comment above and a `publishDir` within twelve lines. And the script's stub
generator only knows how to fake a reads directory or a report folder — a
samplesheet is new, so `make_stubs()` needs a CSV stub and the input guard must
accept a file that exists but points at empty FASTQs, because `-preview` builds
the graph without running a task.

---

## 12. One test suite, parameterised by site

**Decision:** `tests/run_tests.sh` gets `NF` and `PROFILE` environment
overrides and a `PHANTA_ONLY=1` guard, so the same Phanta cases run on AWS
(`PROFILE="-profile amazon" NF=nextflow`) and on FASRC (defaults unchanged).

The suite hardcodes FASRC: the Nextflow binary under `/n/lab_storage`,
`-profile harvard_rc`, and fixture paths under `$HOME`. Phanta must be tested on
AWS first (§9), where none of that exists. The bioBakery cases cannot run on AWS
either — there is no KneadData/MetaPhlAn/HUMAnN stack there — hence the guard
that skips them.

**Rejected:** a separate `tests/run_phanta_tests.sh`. Two suites with copied
helpers drift; the existing one's comments explain at length why it is "the
only test suite", and that is worth preserving. The overrides are three lines.

---

## Corrections made during planning

Recorded because each one changed the specification, and because a reader
comparing this to an earlier draft will otherwise wonder what happened.

**The eligibility mechanism was misread.** The first review covered `Snakefile`
lines 1-120 and 200-364, skipping 120-200 — which contains `rule bracken` and
`rule process_filtered_kraken`. The `.temp` files are not status canaries, as
initially described; they are **substitute outputs** that `rule bracken` and
`rule scale_bracken` copy into place instead of running the tool. Consequence
for the spec: an ineligible sample's observable contract is byte-identical
zero-length products, which the fixtures contain.

**Bracken's second output filename was missed.** `rule bracken` declares two
outputs but passes only `-o`. The Kraken-style species report appears because
Bracken derives its name from the *input* filename. Consequence: the filtered
report must be staged as exactly `{samp}.krak.report.filtered`, and this became
a hard constraint on the nf-core module fit check.

**A premise about existing consortium results was fabricated.** An earlier draft
instructed the implementer to "confirm with the domain scientist which database
the consortium's existing Phanta results were generated against", and suggested
checksumming the database used for those results. There are no such results:
HVP-nextflow contains no reference to Phanta outside these documents, and Phanta
has never been run here. The assumption was never verified before being written
down. Correcting it simplified the plan — the database becomes a fresh choice
rather than an archaeology problem, since byte-identity only requires both sides
of the comparison to agree.

**The database was not pinned at all.** An earlier draft cited the database as a
directory URL (`.../Phanta_DBs/`) rather than an artifact, which hid the
`database_V1` / `unmasked_db_v1.01` contradiction entirely and left the
byte-identity requirement resting on an unspecified input.

**Both documents were written against the pre-v0.0.4 tree.** They were committed
as `5793365` and `main` was merged in immediately after as `2a5a9df`, bringing
the v0.0.4 restructure. The spec was then reconciled against the merged tree.
What had gone stale: the nf-test suite and `nf-test.config` (deleted);
`template-params.yaml` (replaced by `conf/harvard_rc.yaml` and
`conf/tufts_hpc.yaml`); the `processes/` directory (deleted); the `engaging`
profile (deleted); params declared in `nextflow.config` (moved to
`conf/params.config`); a `run_phanta` stage toggle (a category error — Phanta
is a `--workflow` token); a shared `conf/phanta.config` (a layering violation);
"GitHub runner" gates (CI runs no tools); "HPC as the primary path" (the user
chose AWS first); and two subworkflow files for a single-consumer workflow
(against the repo's layering rule). The merge also added obligations the spec
had not known about: `make_diagrams.py --check`, `check_profile_resources.py`,
and the `def env_<stack>` convention. §§9-12 record the resulting decisions.
