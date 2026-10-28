# Phanta samplesheet fixtures

The samplesheets `tests/run_tests.sh` uses for the Phanta cases, against the
tiny database in `../tiny_db/db/` and the synthetic reads in `../tiny_inputs/`.

## Why the `__TINY_INPUTS__` placeholder

A Phanta samplesheet carries read **paths**, and a committed file cannot carry a
path that is correct on every machine. Relative paths do not solve it either:
`run_tests.sh` launches each Nextflow driver from its own directory so that
concurrent runs do not share `.nextflow/`, and Nextflow resolves a relative path
against the launch directory, not against the samplesheet.

So each sheet is committed with the literal token `__TINY_INPUTS__` where the
reads directory belongs, and `run_tests.sh` expands it with `sed` into the
results directory before the run. The committed file stays reviewable and
diffable; the materialised one is disposable.

## The sheets

| File | Cohort | Used by |
| --- | --- | --- |
| `tiny_pe.csv` | 3 paired-end samples, **deliberately not in alphabetical order** (`sample_b`, `sample_a`, `sample_c`) | the main regression case, and the baseline in `../baseline/tiny_pe/` |
| `tiny_se.csv` | the same three samples, R1 only | the single-end path, baseline in `../baseline/tiny_se/` |
| `tiny_alphabetical.csv` | the same three samples, paired, in alphabetical order | the ordering case — same reads as `tiny_pe.csv`, different column order |
| `bad_duplicate_ids.csv` | `sample_a` twice | guard: duplicate sample IDs |
| `bad_mixed_layout.csv` | one paired, one single-end | guard: mixed library layout |
| `bad_header.csv` | `sample_id,read1,read2` | guard: wrong header |

`tiny_pe.csv` is out of alphabetical order on purpose, and so is the baseline
built from it. That is the only thing that distinguishes "columns follow the
samplesheet" from "columns happen to be sorted" — the two are indistinguishable
on an alphabetical cohort, which is why upstream's own `testing/samp_file.txt`
lists `sample19` before `sample18`. The ordering case pairs it with
`tiny_alphabetical.csv`: identical reads, identical results, different column
order, while `total_reads.tsv` rows stay sorted by name in both.

`sample_c` is Bracken-ineligible by construction — too few reads for any species
line to reach `filter_thresh`. Every cohort here includes it, so the zero-byte
substitute path and the failed-sample list are exercised by the main case rather
than only by a special one.
