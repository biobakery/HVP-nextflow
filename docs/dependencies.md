# Dependencies

Where every step's software comes from, and which versions are pinned together
on purpose. The companion to [`docs/architecture.md`](architecture.md),
which says how the pipeline is put together; this page says what it needs to run.

The pipeline itself declares no dependencies: a process runs whatever is on
`PATH` when its task starts. What puts it there is site-specific, so it lives in
the profile — `conf/profiles/<site>.config` — and nowhere else.

## How a process gets its software

| Site | `-profile` | Mechanism |
|---|---|---|
| Harvard FASRC (Cannon) | `harvard_rc` | a `beforeScript` per process group, loading hutlab environment modules |
| Tufts HPC | `tufts_hpc` | tools on `PATH`, or Apptainer where a process declares a `container` |
| This machine | `local` | tools on `PATH` |
| AWS Batch | `amazon` | per-process containers — **incomplete**, six processes only |

Harvard is the only site where the whole pipeline is wired up. There, one Groovy
variable per software stack is defined at the top of the profile and referenced
by the process groups, so a version is changed in one line:

| Stack variable | hutlab module | Provides |
|---|---|---|
| `env_kneaddata` | `rocky8/kneaddata/0.12.0-devel` | `kneaddata`, `kneaddata_read_count_table` |
| `env_metaphlan` | `rocky8/metaphlan4/4.0.6_vOct22_fixed` | `metaphlan`, `merge_metaphlan_tables.py` |
| `env_humann` | `rocky8/humann4/4.0-alpha-1-final` | `humann` and its `humann_*_table` helpers; auto-loads MetaPhlAn |
| `env_baqlava` | `rocky8/baqlava/1.2.0-devel` | `baqlava`; auto-loads HUMAnN 4, anadama2, MetaPhlAn |
| `env_biobakery` | `rocky8/biobakeryworkflows/3.2` | everything below |

`rocky8/biobakeryworkflows/3.2` carries two unrelated groups of dependencies,
which is why so many processes load it:

1. **The assembly tools** — MEGAHIT, bowtie2, samtools, MetaBAT2,
   `jgi_summarize_bam_contig_depths`, CheckM2, PhyloPhlAn, Mash, FastANI,
   assembly-stats — plus the `CHECKM2DB` and `PHYLOPHLAN_PATH` variables. Taking
   them all from one module keeps the assembly component on exactly the
   stack `biobakery_workflows assembly` uses.
2. **Scripts that ship with `biobakery_workflows` rather than with a tool**:
   `count_features.py`, `get_counts_from_humann_logs.py`, `humann_join_tables`.

## A version pin that is not a preference

**The MetaPhlAn / HUMAnN / BAQLaVa triple.** MetaPhlAn `4.0.6_vOct22_fixed` with
`mpa_vOct22_CHOCOPhlAnSGB_202403` is the only combination HUMAnN 4 and BAQLaVa
1.2.0 can both read. HUMAnN 4 cannot yet use the vJan25 database, and
`metaphlan_4/4.2` is built against python3.12, which conflicts with the
python3.10 stack the rest of the profile loads. Change one and all three have to
move. The database halves of the pin are in `conf/databases/harvard_rc.config`;
the module halves are in the profile.

## There are two CheckMs, and they do different jobs

This is the easiest thing in the repo to get wrong, so it is worth stating
plainly.

| | CheckM2 | CheckM 1.2.0 |
|---|---|---|
| Where from | `rocky8/biobakeryworkflows/3.2`, on `PATH` as `checkm2` | **vendored** in `bin/scripts/checkm/`, run as `bin/scripts/checkm.py` |
| Used by | `checkm2` (MAG completeness and contamination) | `abundance` (per-sample MAG abundance) |
| Subcommands used | `predict` | `coverage`, `profile` |
| Extra arguments param | `--checkm_predict_options` | `--checkm_coverage_options` |
| Database | `CHECKM2DB`, exported by the module | none — these two subcommands need no reference data |

CheckM2 is a rewrite, not a new version: it dropped the `coverage` and `profile`
subcommands, which is the step the anadama2 assembly workflow uses to turn a
sample's BAM and bins into an abundance table. So quality assessment uses
CheckM2 and abundance uses a trimmed copy of CheckM 1. Both are needed; neither
substitutes for the other.

The vendored copy is a subset — enough for `coverage` and `profile` and nothing
else. It is CheckM 1.2.0, GPLv3; see `bin/scripts/checkm/README.md` and
`bin/scripts/checkm/LICENSE`.

## Other vendored code

| Path | What it is | Why it is here |
|---|---|---|
| `bin/scripts/checkm/` | CheckM 1.2.0, `coverage` and `profile` only | see above |

## Running the pipeline itself

Nextflow 25.10.4 or newer, as the nf-core template requires (nf-schema 2.7 needs
it), and a JDK 17 or newer; the module loads 21. The hutlab module provides a
24.10.4 launcher and the JDK. `NXF_VER` makes that launcher download and run
25.10.4:

```bash
source /n/lab_storage/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery-workflows-nextflow/0.0.4   # Nextflow launcher + JDK
export NXF_VER=25.10.4                                  # the pipeline needs >= 25.10.4
```

`bin/make_diagrams.py` needs only Python. `bin/check_profile_resources.py`
needs Nextflow and a JVM but no tools and no databases: it only resolves each
profile with `nextflow config`. CI runs both, plus `nextflow run -preview`,
which builds the whole graph without launching a task.

Running them by hand off a login node, where `hutlab load` is a shell function
and so is not available in a script:

```sh
module load jdk/21.0.2-fasrc01
export PATH=/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin:$PATH
bin/check_profile_resources.py && bin/make_diagrams.py --check
```
