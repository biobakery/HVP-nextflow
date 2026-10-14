# Phanta tool image

One image, one digest, referenced by every profile through that profile's own
mechanism. It carries Kraken2, Bracken, and the Python stack the Phanta scripts
run on. It carries no databases — `phanta_db` is staged per site from
`conf/databases/<site>.config`.

Why an image at all, and why pinned this hard, is
[docs/phanta_migration_decisions.md](../../docs/phanta_migration_decisions.md) §2.
The short version: the four published tables must be byte-identical between the
upstream Phanta baseline and this repository's port, every Phanta script writes
through `pandas.to_csv`, and `to_csv` float formatting is version-sensitive. So
both sides of that comparison have to run the same pandas. Upstream pins only
`kraken2` and `snakemake`; the rest is pinned here.

## What is in it

| Package | Version | Why this one |
| --- | --- | --- |
| `kraken2` | 2.1.6 | Upstream's own pin in `phanta_env.yaml`. Needs >= 2.1.2 for `--report-minimizer-data`. |
| `bracken` | 3.1 (reports `v3.0.1`) | Unpinned upstream; this is what the solve resolves to. Output header verified compatible — see below. The conda package version and the string `bracken -v` prints disagree; the package is the authority. |
| `python` | 3.10.22 | Forced by `snakemake=7.15.2`, not chosen freely. |
| `pandas` | 2.3.3 | Byte-identity critical. |
| `numpy` | 2.2.6 | Byte-identity critical. |
| `snakemake` | 7.15.2 | Upstream's own pin. Present so the baseline can be generated **in this image**. |
| `biopython` | 1.88 | Needed by `calculate_genome_size.py` when building a database; `phanta_env.yaml` omits it. |
| `bc` (apt) | distro | `calc_total_reads.sh` refuses to run without it. |
| `procps` (apt) | distro | Nextflow's task metrics shell out to `ps`. |

### How the Python pin was chosen

An unconstrained solve of `kraken2=2.1.6 bracken pandas biopython` today
resolves to **Python 3.15 / pandas 3.0.6 / numpy 2.5.4**. That is not a safe
foundation: pandas 3.x postdates every line of Phanta's code by years, and
changes the default string dtype and the `fillna` downcasting rules that
`merge.py` and `merging_scaled_bracken.py` depend on.

Adding upstream's own `snakemake=7.15.2` pin — which the image needs anyway, to
run the Snakefile for the baseline — constrains the solve to **Python 3.10.22 /
pandas 2.3.3 / numpy 2.2.6**. That is the pin in use. It is derived from an
upstream constraint rather than invented, which is the reason to prefer it over
an arbitrary ceiling.

The alternative considered was forcing `pandas<2`, which yields pandas 1.5.3 /
numpy 1.26.4 — closer to what Phanta's authors actually ran in 2023. It was not
chosen because it requires inventing a constraint upstream never stated. The
risk either way is bounded: both the upstream baseline and the port run in
*this* image, so byte-identity holds regardless of which is picked. What a wrong
pick would cost is divergence from what Phanta's authors intended, not a broken
comparison.

The pipeline scripts were then run end to end against the tiny test database in
this image. They use only `read_csv`, `concat(axis=1).fillna(0)`, `merge`,
`groupby`, and `to_csv`; every `.append` in them is a Python list append, not the
`DataFrame.append` that pandas 2.0 removed.

Two behaviour notes, both verified rather than assumed:

- `calculate_genome_size.py`'s `groupby('Taxonomy').sum()` concatenates a string
  column under pandas 2.x. It only affects `library/taxa_genome_size.txt`, which
  nothing in the pipeline reads — `species_genome_size.txt` is built from
  explicit column selections — so it is recorded rather than patched.
- `scale_bracken.py` raises `ValueError: You are trying to merge on object and
  int64 columns` when `library/species_genome_size.txt` has an all-numeric
  `species_level_taxa` column, because it casts the Bracken side to `str` first.
  **This is not a pandas-2 regression** — pandas 1.5.3 raises the identical
  error, which was checked directly. It is a property of the database, not of
  the image: a real Phanta database holds genomes that resolve no further than
  genus, so that column is object dtype and the merge is string-to-string. It is
  why the tiny test database carries a fifth, genus-tagged genome. See
  [tests/fixtures/phanta/tiny_db/README.md](../../tests/fixtures/phanta/tiny_db/README.md).

## The lock file

`conda-lock.txt` is an `@EXPLICIT` list: 235 packages, each pinned by URL and
checksum. The Dockerfile installs from it and from nothing else.

Version specs alone would not be enough. A conda solve is a moving target, so a
rebuild months from now would quietly drift pandas or numpy and break
byte-identity against a baseline generated today — the exact failure mode the
image exists to prevent.

**Regenerating the lock is a baseline-refresh trigger.** Do not do it to pick up
an unrelated package; do it deliberately, and expect to regenerate the baselines
in `tests/fixtures/phanta/baseline/` and re-record the resolved versions in
`docs/phanta_testing.md`.

```bash
docker run --rm --platform linux/amd64 mambaorg/micromamba:1.5.8 bash -c "
  micromamba create -y -n phanta -c conda-forge -c bioconda \
    python=3.10.22 pandas=2.3.3 numpy=2.2.6 \
    kraken2=2.1.6 bracken snakemake=7.15.2 biopython >/dev/null 2>&1 &&
  micromamba env export -n phanta --explicit
" > /tmp/lock.txt
# then replace everything below the header comment in conda-lock.txt
```

## Build

Always `linux/amd64`, including on an arm64 Mac: AWS Batch is amd64, and
bioconda publishes no arm64 build of kraken2. On Apple silicon this runs under
emulation and takes a few minutes.

```bash
docker build --platform linux/amd64 -t phanta:0.0.5-dev containers/phanta
```

The build ends with a smoke test — `kraken2 --version`, `bracken -v`, a pandas /
numpy / Bio import, and a check that `kraken2-build`, `bracken-build`, and `bc`
are on `PATH`. A build that cannot start its own tools fails rather than ships.

## Push, and the digest

**Not yet done — this needs AWS credentials and an ECR repository, neither of
which exists yet.** See "Inputs a fresh agent cannot produce" item 1 in the
spec. When the repository exists:

```bash
ECR=public.ecr.aws/<registry-alias>/phanta
aws ecr-public get-login-password --region us-east-1 \
  | docker login --username AWS --password-stdin public.ecr.aws

docker tag phanta:0.0.5-dev "$ECR:0.0.5"
docker push "$ECR:0.0.5"

# The digest is the deliverable. Record THIS, not the tag.
docker inspect --format='{{index .RepoDigests 0}}' "$ECR:0.0.5"
```

Then put that digest — never the tag — in three places, per convention 4:

- `conf/profiles/aws.config`: `container = 'public.ecr.aws/…@sha256:…'` in each
  Phanta `withName` block.
- `conf/profiles/harvard_rc.config`: the Apptainer pull inside `def env_phanta`.
- `docs/phanta_testing.md`, alongside the resolved tool versions, as part of the
  baseline's provenance.

A tag is mutable and two sites resolving the same tag at different times can get
different images, which is precisely the failure byte-identity would then hide.

## Verified in this image

Resolved versions, read out of the built image. These are what
`docs/phanta_testing.md` must record against any baseline generated here:

```
kraken2    2.1.6
bracken    3.1  (conda package; `bracken -v` prints "Bracken v3.0.1")
python     3.10.22
pandas     2.3.3
numpy      2.2.6
biopython  1.88
snakemake  7.15.2
```

**Bracken 3.x output header — checked, compatible.** `scale_bracken.py` reads
the `.bracken` table by column *name*, and Phanta was written against Bracken
2.x, so a rename would have broken the port silently. `est_abundance.py:437` in
this image writes:

```
name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads  new_est_reads  fraction_total_reads
```

— identical to 2.x, and it covers every column `scale_bracken.py` names
(`taxonomy_id`, `new_est_reads`, `fraction_total_reads`). No Bracken downgrade
needed. `fraction_total_reads` is emitted at `%0.5f`, but `scale_bracken.py`
recomputes that column before writing, so the printed precision never reaches a
published table.

## Open items

- **Whether to keep pandas 2.3.3 or drop to 1.5.3 remains genuinely open.** The
  one concrete failure found so far is not pandas-version dependent (above), so
  there is no evidence forcing the change. 1.5.3 is still the more conservative
  choice, being contemporaneous with Phanta's code. Decide before the phase 1
  baseline is generated, because changing it afterwards invalidates the baseline.
- **ECR push and digest capture**, above. Blocked on AWS credentials and an ECR
  repository.
- **The image is built but unpushed.** Until it has a digest, no profile can
  reference it, so `conf/profiles/aws.config` and `def env_phanta` in
  `conf/profiles/harvard_rc.config` cannot be written against a real value.
