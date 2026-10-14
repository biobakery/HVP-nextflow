# Phanta tiny test database

A complete, working Phanta database of five genomes — about 9 MB including the
synthetic reads, against the real database's 20-25 GB. Built by
[build.sh](build.sh); the built copy under `db/` is committed so that a test
run needs neither Docker nor network access.

Why it exists is [decisions §8](../../../docs/phanta_migration_decisions.md):
a regression run finishes in seconds instead of hours, neither site has to stage
25 GB per test run, and — the part the real database cannot do at all — it
contains a sample that Bracken refuses to process, which is the only way to
exercise the zero-byte-substitute path that the port depends on.

**Open question 3 in the spec is answered: committed, not downloaded.** At 9 MB
the hosting-plus-checksum-fetch machinery is not worth it.

## What is in it

| Class | Organism | Accession | Tagged taxid |
| --- | --- | --- | --- |
| viral | *Sinsheimervirus* phiX174 | `NC_001422.1` | 10847 (species) |
| bacterial | *Mycoplasmoides genitalium* G37 | `NC_000908.2` | 243273 (strain) |
| archaeal | *Methanocaldococcus jannaschii* DSM 2661 | `NC_000909.1` | 243232 (strain) |
| eukaryotic | *S. cerevisiae* S288C chromosome I | `NC_001133.9` | 559292 (strain) |
| genus-only | *Mycoplasmoides pneumoniae* M129 | `NC_000912.1` | 2995234 (genus) |

The first four are one per threshold class, so every
`cov_thresh_*`/`minimizer_thresh_*` branch in `filter_kraken_reports.py` runs
against real data. Three are tagged at strain level and one at species level,
which exercises both arms of the "lowest rank for this species" logic — an S
line followed by S1 is handled differently from an S line standing alone.

The fifth is there for a less obvious reason, below.

## Three traps this fixture exists to avoid

Each of these was hit while building it, and each failed *silently* or
misleadingly. They are the reason the taxonomy is committed rather than
downloaded, and the reason there are five genomes rather than four.

**1. `superkingdom` no longer exists in NCBI's taxonomy.**
`filter_kraken_reports.py` looks up the literal rank string `superkingdom`.
NCBI has since renamed that rank to `domain` (and `Viruses` to
`acellular root`). Appendix D of the spec says to run
`kraken2-build --download-taxonomy`; **do not**. With a freshly downloaded
taxonomy every superkingdom lookup fails, falls through the catch-all `else` at
`filter_kraken_reports.py:248` that *keeps* the species, and filtering silently
does nothing — while the output still looks completely normal. Phanta's entire
purpose is that filtering step. `taxonomy/nodes.dmp` and `names.dmp` are
therefore hand-authored, committed, and use `superkingdom`, matching what the
real Phanta database ships. A download would also make the database differ run
to run, since NCBI's dump is a moving target.

**2. NCBI dump lines end with `\t|`, not `\t|\t`.**
Kraken2 discards the final two characters of each line before reading the last
field. A line ending `\t|\t` leaves a stray `\t|` glued to the name class, so
the `scientific name` test never matches and **every name in `inspect.out`
comes out empty** — with no error anywhere. Phanta reads taxon names from
`inspect.out`, so the merged tables would carry empty lineage names.

**3. A database where every genome resolves to a species is unrepresentative.**
This is what the fifth genome is for. `calculate_genome_size.py` writes the
literal string `error - taxid above species` into `species_level_taxa` for any
genome that does not resolve to a species. With only the first four genomes
every row is numeric, pandas types that column `int64`, and `scale_bracken.py`
dies on:

```
ValueError: You are trying to merge on object and int64 columns for key 'taxonomy_id'
```

because it casts the Bracken side to `str` first. The real database never hits
this — it holds many genomes classified no further than genus, so the column is
object dtype and the merge is string-to-string. **This is not a pandas-version
problem:** pandas 1.5.3 raises the identical error. Tagging
`NC_000912.1` at genus makes the fixture behave like the real database.

## Rebuilding

Run inside the pinned image, so kraken2, bracken, dustmasker, pandas and
biopython are exactly the versions the baselines were generated with:

```bash
docker build --platform linux/amd64 -t phanta:0.0.5-dev containers/phanta

docker run --rm --platform linux/amd64 -v "$PWD":/w -w /w phanta:0.0.5-dev \
  bash tests/fixtures/phanta/tiny_db/build.sh -o /w/build/tiny
```

Then copy `build/tiny/reads/*.fq.gz` over `../tiny_inputs/`, and copy from
`build/tiny/tiny_db/` into `db/` exactly the nine files the spec's "Required
database files" table lists:

```
hash.k2d  taxo.k2d  opts.k2d  seqid2taxid.map  database100mers.kmer_distrib
inspect.out  taxonomy/nodes.dmp  taxonomy/names.dmp  library/species_genome_size.txt
```

Everything else a build leaves behind is dropped: `database.kraken` and
`database100mers.kraken` are Bracken build intermediates; `library/added/` holds
the source FASTAs, which only `calculate_genome_size.py` reads at *build* time
and which kraken2-build names randomly, so committing them would add churn and
half the fixture's size for nothing; `taxa_genome_size.txt` is read by nothing at
all. This also mirrors how the real database is used — upstream ships the
`library/` FASTAs as a separate download.

Rebuilding invalidates the committed baselines. Treat it the same way as
regenerating the image lock.

## Identity

```
inspect.out sha256  554c323f5357f20cfd5372350676d1f55e3fc5a209015c0310c2180e7e5e299e
```

That value is what `phanta_db_expected_inspect_sha256` must be set to for any
test run against this database, and what the `test21_phanta_manifest` case
greps for. `inspect.out` is content-derived, so it fingerprints the build.

## Samples

| Sample | Composition | Role |
| --- | --- | --- |
| `sample_a` | 3450 pairs across all five genomes | ordinary cohort member |
| `sample_b` | 2560 pairs, different proportions | ordinary cohort member |
| `sample_c` | 4 pairs | **Bracken-ineligible** |

`sample_c` has too few reads for any species to reach `filter_thresh` (10), so
`process_filtered_kraken.py` writes its zero-byte substitutes, Bracken never
runs, and the sample lands in `samples_that_failed_bracken.txt`. Reads are
generated by [../tiny_inputs/make_reads.py](../tiny_inputs/make_reads.py),
which fixes the PRNG seed and writes gzip headers with `mtime=0` so the files
are byte-identical between runs.
