# Testing Phanta

How to run the Phanta cases, what the committed fixtures are, and how to
regenerate a baseline when one legitimately has to change.

The design behind all of this is in
[phanta_migration_spec.md](phanta_migration_spec.md) — "Comparison model" for
why there are three distinct comparison operations, and Appendix F for the
traps that cost time during the port. This document is the procedure.

---

## The bar: byte-identity

Phanta's four published tables must be **byte-identical** to the baseline, not
approximately equal. There are no numeric tolerances anywhere in the
comparison, and adding one would defeat its purpose: every script in the chain
writes through `pandas.to_csv`, whose float formatting varies between pandas
versions, so a difference in the last decimal place is evidence that something
about the environment changed. That is exactly what the comparison exists to
catch.

A mismatch is therefore always investigated, never tolerated. The three usual
causes, in order of likelihood:

1. **A different database.** Compare `pipeline_info/phanta_db_manifest.txt`
   between the two runs; `inspect.out.sha256` is the field that matters.
2. **A different container.** `pandas.version` and `numpy.version` are in the
   manifest for this reason.
3. **A reordered samplesheet.** Merged-table columns follow the samplesheet;
   `total_reads.tsv` rows are sorted by name. The two orders are different on
   purpose — see "Output ordering contract" in the spec.

---

## Running the cases

```bash
# Build the image once (linux/amd64 even on an arm64 machine).
docker build --platform linux/amd64 -t phanta:0.0.5-dev containers/phanta

# Phanta only — no cluster, no module stack, no 25 GB database.
PHANTA_ONLY=1 NF=nextflow PROFILE="-profile local_docker" bash tests/run_tests.sh
```

Other ways in:

| Invocation | Runs |
| --- | --- |
| `PHANTA_ONLY=1 …` | the Phanta cases alone, then the summary |
| `PHANTA=1 bash tests/run_tests.sh` | the bioBakery cases, then the Phanta cases, one summary |
| `bash tests/run_tests.sh` | the bioBakery cases only; Phanta is skipped with a note |

Phanta is opt-in rather than on by default because the two halves of the suite
need different things. The bioBakery cases need the hutlab module environment
and a cluster; the Phanta cases need a container runtime and the committed
fixtures, and nothing else. Neither host has the other's prerequisites.

`NF` and `PROFILE` are both overridable so the same suite runs at more than one
site. On AWS, once the infrastructure exists:

```bash
PHANTA_ONLY=1 NF=nextflow PROFILE="-profile amazon" bash tests/run_tests.sh
```

### What the cases cover

The cases live in [tests/phanta_cases.sh](../tests/phanta_cases.sh), which
`run_tests.sh` sources so that there is still one suite and one summary.

| Case | Asserts |
| --- | --- |
| 17 | paired-end cohort byte-identical to the baseline; Kraken2's report has eight columns; Bracken's implicit second output is under its upstream name |
| 18 | single-end cohort byte-identical to its own baseline |
| 19 | the Bracken-ineligible sample gets zero-byte products, is listed as failed, and is absent from the merged tables |
| 20 | merged-table columns follow the samplesheet while `total_reads.tsv` rows stay sorted |
| 21 | the manifest records the database checksum, tool versions and params |
| 22 | a wrong expected checksum fails the run with **no Kraken2 task started** |
| 23 | each cohort barrier ran exactly once |
| 24 | `-resume` caches every task and re-runs none |
| 25a–g | seven input guards: duplicate IDs, mixed layout, bad header, no database, incomplete database, a read length the database has no distribution for, no samplesheet |

Three of these are worth knowing the reasoning for.

**Case 20 needs two samplesheets.** On an alphabetical cohort, "columns follow
the samplesheet" and "columns happen to be sorted" are indistinguishable. So
`tiny_pe.csv` lists `sample_b` before `sample_a`, `tiny_alphabetical.csv` lists
the same reads in order, and the case asserts the columns differ between the
two runs while the `total_reads.tsv` rows stay sorted in both. Upstream's own
`testing/samp_file.txt` lists `sample19` before `sample18` for the same reason.

**Case 22's checksum is all zeroes on purpose.** An earlier version of the
fail-fast guard was silently skipped, because Nextflow coerces a numeric-looking
CLI value to a number and Groovy's `?:` treats `0` as false. An all-digit value
is the shape that broke it.

**Case 25's guards each assert that no Kraken2 task started**, not merely that
the run failed. Failing is easy; failing *before* an hour of classification is
the actual contract.

---

## The fixtures

Everything needed to run the cases is committed. No network access, no Docker
pull beyond the image, no database download.

| Path | What | Size |
| --- | --- | --- |
| `tests/fixtures/phanta/tiny_db/` | the tiny Kraken2 + Bracken database, its build script, and pinned accessions | 5.4 MB |
| `tests/fixtures/phanta/tiny_inputs/` | three synthetic paired-end samples and the generator | 0.4 MB |
| `tests/fixtures/phanta/samplesheets/` | six samplesheets, three valid and three deliberately broken | 2 KB |
| `tests/fixtures/phanta/baseline/` | the expected output trees, one per layout | 104 KB |

### The tiny database

Five genomes, one per threshold class in `filter_kraken_reports.py`, plus a
fifth classified only to genus. Full detail and the reasoning are in
[tests/fixtures/phanta/tiny_db/README.md](../tests/fixtures/phanta/tiny_db/README.md)
and in Appendix D of the spec, which records three amendments without which the
build produces a database that is broken in ways nothing reports.

Its identity, which the cases pass as `--phanta_db_expected_inspect_sha256`:

```
inspect.out  sha256  554c323f5357f20cfd5372350676d1f55e3fc5a209015c0310c2180e7e5e299e
```

It ships `database100mers.kmer_distrib` and no other, so
`--phanta_read_length 100` is mandatory. Case 25f asserts that getting this
wrong produces a message naming the parameter rather than a tool error.

`sample_c` is Bracken-ineligible by construction: too few reads for any species
line to reach `filter_thresh`. It is in every cohort, so the zero-byte
substitute path is exercised by the ordinary cases rather than only by a
special one.

### What produced the committed baselines

Recorded per baseline in `tests/fixtures/phanta/baseline/<name>.provenance.txt`:

| | |
| --- | --- |
| Image | `phanta:0.0.5-dev`, built from `containers/phanta/` |
| Kraken2 | 2.1.6 |
| Bracken | conda package **3.1**; `bracken -v` self-reports **v3.0.1** |
| Python / pandas / numpy | 3.10.22 / 2.3.3 / 2.2.6 |
| Database | the committed tiny database, checksum above |
| Read length | 100 |

The Bracken discrepancy is the tool's own doing: the bioconda package is
`bracken-3.1-h9948957_0` and its `-v` string was never bumped. Both numbers are
correct, for different things. The manifest records what `-v` says, so the
manifest will say 3.0.1.

**These baselines were generated locally, not on AWS.** The spec's acceptance
criterion 10 calls for AWS-generated baselines, and that remains the bar for
the full database. Generating the tiny-database baseline locally was a
deliberate call: it comes out of the same pinned image, and leaving the entire
regression suite blocked on infrastructure that has not been handed over buys
nothing. Re-verify on Batch when AWS exists; if it differs, that is a finding
about the environment and belongs in the manifest comparison, not in a
tolerance.

Two independent runs produced byte-identical trees before the baseline was
committed, which is the spec's own two-run requirement applied to the tiny
database.

---

## Regenerating a baseline

**Only when something intentionally changed**: the container, the database, the
tool versions, the inputs, the params, or the ordering contract. A failing
comparison is not by itself a reason to regenerate — it is a reason to find out
why, and regenerating to make a test pass destroys the evidence that the test
existed to preserve.

```bash
docker build --platform linux/amd64 -t phanta:0.0.5-dev containers/phanta

OUT=$(mktemp -d)
nextflow run main.nf -profile local_docker --workflow phanta \
    --phanta_samplesheet <materialised sheet> \
    --phanta_db        tests/fixtures/phanta/tiny_db/db \
    --phanta_read_length 100 \
    --phanta_db_expected_inspect_sha256 554c323f5357f20cfd5372350676d1f55e3fc5a209015c0310c2180e7e5e299e \
    --outdir "$OUT" --log_versions false
```

The samplesheets are committed with a `__TINY_INPUTS__` placeholder rather than
real paths — see
[tests/fixtures/phanta/samplesheets/README.md](../tests/fixtures/phanta/samplesheets/README.md).
Materialise one the way `tests/phanta_cases.sh` does:

```bash
sed "s|__TINY_INPUTS__|$PWD/tests/fixtures/phanta/tiny_inputs|g" \
    tests/fixtures/phanta/samplesheets/tiny_pe.csv > /tmp/sheet.csv
```

Then:

1. **Run it twice**, into two different output and work directories, and
   confirm the two agree. A baseline that does not reproduce itself is not a
   baseline.
2. **Diff the new tree against the old one** with
   `bin/phanta_compare_outputs.sh` and account for every difference in the
   commit message. "Regenerated" is not an explanation.
3. Copy `$OUT/phanta` over `tests/fixtures/phanta/baseline/<name>` and refresh
   `<name>.provenance.txt` from `$OUT/pipeline_info/phanta_db_manifest.txt`.

---

## The full database

Not yet done — it needs the database pinned (open question 1) and the AWS
infrastructure. The procedure when it is:

1. Pin one artifact from `http://ab_phanta.os.scg.stanford.edu/Phanta_DBs/`,
   compute its SHA256 **at download**, and record filename, checksum and date
   here. Upstream publishes no checksums, so this is the only record.
2. Confirm `library/species_genome_size.txt` ships in the main tarball. If it
   only comes with `database_V1.library.tar.gz`, that add-on becomes mandatory
   and the staging footprint grows.
3. Set `phanta_db` and `phanta_db_expected_inspect_sha256` in the site's
   `conf/databases/` file.
4. Run upstream Phanta at commit `ecbfcf9` in the pinned image, per Appendix C
   of the spec. Never write into the Phanta checkout.
5. Run it a second time independently; require byte-identity apart from
   `.krak`.
6. Commit to `tests/fixtures/phanta/baseline/full/` with checksums, and add a
   case that uses it.

Run with `--phanta_db_checksum_full` for a baseline: it hashes the 20+ GB
`hash.k2d` as well as `inspect.out`. Ordinary runs leave it off, because that
is a complete read of the index every time.

---

## Paired-run comparison against upstream

For validating a new dataset, rather than regression. Run upstream and this
pipeline from identical inputs, params, database, tool versions and container
into two directories, then:

```bash
bin/phanta_compare_outputs.sh <upstream_out> <our_out>/phanta
```

It will always report a handful of **structural** differences, which are
publishing-contract choices rather than content. The script prints the list
itself; the table is in Appendix F of the spec. In short: upstream publishes
three zero-byte canary files and a `pipeline_completed.txt` that this pipeline
does not, and it moves the filtering decisions into `classification/intermediate/`
where this pipeline publishes them in place.

Everything both sides publish in the same place is byte-identical. To see only
that:

```bash
diff -r -x '*.krak' -x pipeline_info -x '*.temp' \
     -x pipeline_completed.txt -x intermediate <upstream_out> <our_out>/phanta
```

Do not commit the outputs of a paired run.

---

## What is not covered yet

Honest list, so nobody mistakes a green suite for a finished migration.

- **Nothing has run on AWS Batch.** Every result above is from
  `-profile local_docker` on one machine.
- **Nothing has run at FASRC.** `conf/profiles/harvard_rc.config` carries
  `phanta_image = null`, because there is no ECR digest to point it at yet.
- **No full-database run, and no real cohort.** Both need the database pinned
  and the infrastructure handed over.
- **The case results above are from Nextflow 26.10.4 with
  `NXF_SYNTAX_PARSER=v1`**, not the pinned 24.10.4, because Nextflow 26 cannot
  parse `conf/base.config` at all and the 24.10.4 build available here is a
  bare jar. Re-confirm on a properly installed 24.10.4 — particularly the
  resume case, since a jar started with `java -cp` puts Nextflow in `mode=dev`,
  which keeps no resume cache for *any* pipeline.

### A warning about Nextflow 26 and the diagrams

Under Nextflow 26, `bin/make_diagrams.py --check` reports every diagram stale —
including workflows nobody has touched — and regenerating under 26 **corrupts
them**: node descriptions are dropped and identifiers are mangled
(`_kraken2_("<b>"kraken2"</b>")`). The committed diagrams are correct and the
check passes on the pinned 24.10.4. Do not commit diagrams regenerated under
Nextflow 26, however loudly the check asks for it.
