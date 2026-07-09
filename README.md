# biobakery-nextflow

> Nextflow implementation of the [bioBakery](https://github.com/biobakery) metagenomics suite,
> ported from the [anadama2 biobakery-workflows](https://github.com/biobakery/biobakery_workflows)
> with multi-cluster support, nf-core-style modules, and an SGB assembly pipeline.

![Static Badge](https://img.shields.io/badge/Author-Kevin_Bonham_PhD-purple)
![Static Badge](https://img.shields.io/badge/Author-Guilherme_Fahur_Bottino_PhD-purple)
![Static Badge](https://img.shields.io/badge/Author-Danielle_Pinto-purple)

---

## Table of Contents

1. [Overview](#overview)
2. [Tools & Versions](#tools--versions)
3. [Installation & Environment Setup](#installation--environment-setup)
   - [Harvard FASRC (Cannon)](#harvard-fasrc-cannon)
   - [Tufts HPC](#tufts-hpc)
   - [AWS Batch](#aws-batch)
   - [Local / Docker](#local--docker)
4. [Quick Start](#quick-start)
5. [Use Cases](#use-cases)
6. [All Parameters](#all-parameters)
7. [Compute Profiles](#compute-profiles)
8. [Changing Databases](#changing-databases)
9. [Changing Resources](#changing-resources)
10. [Output Structure](#output-structure)
11. [Testing](#testing)

---

## Overview

Select a workflow with `--workflow`:

| `--workflow` | Status | Description |
|---|---|---|
| `mgx` *(default)* | ✅ Ready | Whole metagenome shotgun: QC + taxonomy + function |
| `sgb_pipeline` | ✅ Ready | MAG assembly → binning → SGB clustering (ported from anadama2 `feature/sgb_pipeline`) |
| `mgx_mtx` | 🚧 Partial | Paired MGX + MTX (runs MGX half) |
| `mtx` | 🔜 Stub | Metatranscriptome |
| `16s` | 🔜 Stub | 16S amplicon (DADA2 / QIIME2 planned) |
| `vis` | 🔜 Stub | Quarto report generation |
| `stats` | 🔜 Stub | MaAsLin2 / LEfSe / diversity |

---

## Tools & Versions

| Tool | Purpose | Supported versions |
|---|---|---|
| [KneadData](https://github.com/biobakery/kneaddata) | Host decontamination + QC trimming | 0.12+ |
| [MetaPhlAn](https://github.com/biobakery/MetaPhlAn) | Species-level taxonomic profiling | v3.1, v4.0.6, v4.1.1 |
| [HUMAnN](https://github.com/biobakery/humann) | Functional pathway profiling | 3.7, 4.0-alpha |
| [BAQLaVa](https://github.com/biobakery/baqlava) | Viral profiling | 0.5+ |
| [StrainPhlAn](https://github.com/biobakery/MetaPhlAn) | Strain-level profiling (SGB mode) | bundled with MetaPhlAn 4 |
| [MEGAHIT](https://github.com/voutcn/megahit) | De novo metagenome assembly | 1.2+ |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | MAG binning | 2.15+ |
| [CheckM2](https://github.com/chklovski/CheckM2) | MAG quality assessment | 1.0+ |
| [PhyloPhlAn](https://github.com/biobakery/phylophlan) | Phylogenetic MAG placement | 3.0+ |
| [Mash](https://github.com/marbl/Mash) | Pairwise genome distance (SGB clustering) | 2.0+ |

---

## Installation & Environment Setup

### Harvard FASRC (Cannon)

Tools are pre-installed via the hutlab module system. No containers needed.

```sh
# One-time setup (run once per login shell or add to ~/.bashrc)
source /n/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery-workflows-nextflow/0.0.1

# Verify Nextflow is available
nextflow -version
```

Database paths are automatically loaded when you use `-profile harvard_rc`
(from `conf/databases/harvard_rc.config`). No params file needed for defaults.

---

### Tufts HPC

Tools are provided via Apptainer containers already configured in the `tufts_hpc` profile.

```sh
# Add Nextflow and shared binaries to PATH
module load nextflow
export PATH="/cluster/tufts/bonhamlab/shared/bin:$PATH"

# Database paths are loaded from conf/databases/tufts.config
```

Create a params file for your run (see `template-params.yaml` for all options):

```yaml
# my-run-params.yaml
readsdir: /path/to/my/fastqs
outdir:   /path/to/results
paired_end: true
filepattern: "*_R{1,2}*.fastq.gz"
```

```sh
nextflow run main.nf -profile tufts_hpc -params-file my-run-params.yaml
```

---

### AWS Batch

Requires prior setup of AMIs, job queues, S3 buckets, and ECR containers.
ECR container URIs are already configured in `conf/profiles/aws.config`.

```sh
nextflow run main.nf -profile amazon -params-file params.yaml
```

See `conf/profiles/aws.config` for queue and container assignments per process.

---

### Local / Docker

Requires all tools installed in your `$PATH` (conda recommended).

```sh
# Install via conda/mamba
mamba create -n biobakery -c biobakery -c conda-forge \
  kneaddata metaphlan=4 humann nextflow
conda activate biobakery

nextflow run main.nf -profile local \
  --readsdir /path/to/fastqs \
  --outdir   /path/to/results \
  --human_genome  /path/to/kneaddata/hg38 \
  --metaphlan_db  /path/to/metaphlan_db \
  --humann_db     /path/to/humann_db
```

> **Note:** MetaPhlAn requires ≥ 15 GB RAM. A 16 GB laptop is often insufficient.
> Prefer HPC or cloud for production runs.

---

## Quick Start

### Minimal MGX run (paired-end)

```sh
nextflow run main.nf \
  --workflow     mgx \
  --readsdir     /path/to/fastqs \
  --outdir       results \
  --human_genome /path/to/kneaddata/hg38 \
  --metaphlan_db /path/to/metaphlan_db \
  --humann_db    /path/to/humann_db
```

### Single-end reads

```sh
nextflow run main.nf --workflow mgx \
  --paired_end false \
  --filepattern "*.fastq.gz" \
  --readsdir /path/to/fastqs \
  [other required params]
```

### Resume an interrupted run

```sh
nextflow run main.nf [same params as original run] -resume
```

### Using a params file (recommended for reproducibility)

```sh
# Copy template and fill in your paths
cp template-params.yaml my-params.yaml
# Edit my-params.yaml ...

nextflow run main.nf -params-file my-params.yaml
```

---

## Use Cases

### 1. Standard MGX: taxonomy + function only

Default settings. Both MetaPhlAn and HUMAnN run; regroup/rename/merge are on.

```sh
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results
```

---

### 2. Taxonomy only — with MetaPhlAn options

Skip HUMAnN and tune MetaPhlAn output format, alignment length, and database:

```sh
# Marker abundance table instead of relative abundance (default: rel_ab_w_read_stats)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_functional_profiling false \
  --metaphlan_analysis_type marker_ab_table

# Stricter minimum read alignment length (default: 70 bp)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_functional_profiling false \
  --metaphlan_read_min_len 100

# Switch to a different MetaPhlAn database index
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --metaphlan_db   /path/to/other/metaphlan_databases \
  --metaphlan_index mpa_vOct22_CHOCOPhlAnSGB_202403 \
  --metaphlan_version metaphlan_v4
```

> On Harvard FASRC hutlab builds, pass `--metaphlan_version metaphlan_v3` — the hutlab binary uses v3-style CLI flags (`--bowtie2db`) even though it is MetaPhlAn 4.

---

### 3. Functional profiling — HUMAnN options

Tune the regrouping target, skip internal searches, or disable post-processing:

```sh
# KEGG orthologs instead of MetaCyc reactions (default: uniref90_rxn)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --humann_regroup_grouping uniref90_ko

# Available groupings: uniref90_rxn | uniref90_ko | uniref90_eggnog | uniref90_pfam

# Skip MetaPhlAn prescreen inside HUMAnN (useful when species profile already exists)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --humann_bypass_prescreen true

# Skip nucleotide search — go straight to translated search (faster, lower sensitivity)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --humann_bypass_nucleotide_search true

# Disable post-processing (no regroup / rename / merge tables)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_humann_regroup false \
  --run_humann_rename  false \
  --run_humann_merge   false
```

---

### 4. Add viral profiling — BAQLaVa options

```sh
# Standard viral profiling
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_viral_profiling true

# Bypass bacterial depletion — needed for test/tiny samples with 0 prescreen species
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_viral_profiling true \
  --baqlava_bypass_depletion true

# Use a custom BAQLaVa database (default: bundled db)
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_viral_profiling true \
  --baqlava_db /path/to/custom/baqlava_db
```

---

### 5. Add strain-level profiling (StrainPhlAn)

StrainPhlAn runs after MetaPhlAn and uses the compressed SAM output.
By default it profiles all detected clades.

```sh
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_strain_profiling true
```

Target specific clades only:

```sh
  --run_strain_profiling true \
  --strainphlan_clades "s__Akkermansia_muciniphila,s__Bacteroides_fragilis"
```

---

### 6. Full MGX with all optional modules

```sh
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   results \
  --run_viral_profiling  true \
  --run_strain_profiling true
```

---

### 7. SGB pipeline (MAG assembly → binning → SGB clustering)

Requires MEGAHIT, MetaBAT2, CheckM2, PhyloPhlAn, Mash, and the PhyloPhlAn database.

```sh
nextflow run main.nf -profile harvard_rc \
  --workflow      sgb_pipeline \
  --readsdir      /path/to/kneaddata_output \
  --outdir        results \
  --run_qc        false \
  --phylophlan_path /path/to/phylophlan_databases
```

Lower quality thresholds for more permissive binning:

```sh
  --sgb_completeness  30 \
  --sgb_contamination 15
```

Cross-dataset abundance (align reads against all MAGs in the cohort):

```sh
  --sgb_abundance_type by_dataset
```

---

### 8. Skip QC (start from already-cleaned reads)

```sh
nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/kneaddata_output \
  --outdir   results \
  --run_qc   false
```

---

### 9. HUMAnN v4 (alpha)

Switch both the software module and the database:

```sh
nextflow run main.nf -profile harvard_rc \
  --humann_version humann_v4a \
  --humann_db /n/huttenhower_lab/tools/nextflow/databases/humann4 \
  --readsdir /path/to/fastqs \
  --outdir   results
```

> Also update `beforeScript` in `conf/profiles/harvard_rc.config` to load
> `rocky8/humann4/4.0-alpha-1` for the `humann` process.

---

### 10. Mouse / ribosomal RNA decontamination

Change only the KneadData reference; everything else stays the same.

```sh
# Mouse C57BL reference
nextflow run main.nf -profile harvard_rc \
  --human_genome /n/huttenhower_lab/data/kneaddata_databases/mouse_C57BL \
  --readsdir /path/to/fastqs --outdir results

# Ribosomal RNA (for rRNA depletion check)
nextflow run main.nf -profile harvard_rc \
  --human_genome /n/huttenhower_lab/data/kneaddata_databases/ribosomal_RNA/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2 \
  --readsdir /path/to/fastqs --outdir results
```

---

### 11. Demo run (bundled test sample)

```sh
# Harvard FASRC — single-end demo
nextflow run main.nf -profile harvard_rc \
  --readsdir    test/single_end_rawfastq \
  --filepattern "*.fastq.gz" \
  --paired_end  false \
  --outdir      test/results

# Same demo + viral profiling
nextflow run main.nf -profile harvard_rc \
  --readsdir    test/single_end_rawfastq \
  --filepattern "*.fastq.gz" \
  --paired_end  false \
  --outdir      test/results \
  --run_viral_profiling true \
  --baqlava_bypass_depletion true
```

---

## All Parameters

### Workflow & I/O

| Parameter | Default | Description |
|---|---|---|
| `--workflow` | `mgx` | Workflow: `mgx \| mtx \| mgx_mtx \| 16s \| vis \| stats \| sgb_pipeline` |
| `--readsdir` | *required* | Directory containing input FASTQ files |
| `--outdir` | `results` | Output directory |
| `--paired_end` | `true` | `true` for paired-end, `false` for single-end |
| `--filepattern` | `*_R{1,2}*.fastq.gz` | Glob to match reads. Paired-end must use `{1,2}`. |

### Tool toggles

| Parameter | Default | Description |
|---|---|---|
| `--run_qc` | `true` | Run KneadData QC. Set `false` to start from cleaned reads. |
| `--run_taxonomic_profiling` | `true` | Run MetaPhlAn |
| `--run_functional_profiling` | `true` | Run HUMAnN |
| `--run_viral_profiling` | `false` | Run BAQLaVa viral profiling |
| `--run_baqlava` | `false` | Alias for `--run_viral_profiling` (backward compat) |
| `--run_strain_profiling` | `false` | Run StrainPhlAn SGB mode |

### QC options (KneadData)

| Parameter | Default | Description |
|---|---|---|
| `--human_genome` | *required* | KneadData Bowtie2 reference directory (hg38, mouse, rRNA, etc.) |
| `--kneaddata_bypass_trim` | `false` | Skip Trimmomatic trimming step |
| `--kneaddata_remove_intermediate_files` | `true` | Delete intermediate files to save disk |

### Taxonomic profiling options (MetaPhlAn)

| Parameter | Default | Description |
|---|---|---|
| `--metaphlan_db` | *required* | Directory containing MetaPhlAn bowtie2 index + pkl |
| `--metaphlan_index` | `mpa_vJun23_CHOCOPhlAnSGB_202307` | Database index name (must exist inside `metaphlan_db`) |
| `--metaphlan_version` | `metaphlan_v4` | CLI flag style: `metaphlan_v3` (uses `--bowtie2db`) or `metaphlan_v4` (uses `--db_dir`). **Harvard FASRC hutlab builds use `metaphlan_v3`.** |
| `--metaphlan_analysis_type` | `rel_ab_w_read_stats` | `-t` flag: `rel_ab \| rel_ab_w_read_stats \| marker_ab_table` |
| `--metaphlan_read_min_len` | `70` | Minimum read alignment length |

### Functional profiling options (HUMAnN)

| Parameter | Default | Description |
|---|---|---|
| `--humann_db` | *required* | Parent directory; pipeline appends `/chocophlan`, `/uniref`, `/utility_mapping` |
| `--humann_version` | `humann_v37` | Output naming: `humann_v37` or `humann_v4a` |
| `--humann_bypass_prescreen` | `false` | Skip MetaPhlAn-based taxonomic prescreen |
| `--humann_bypass_nucleotide_search` | `false` | Skip nucleotide (ChocoPhlAn) search |
| `--run_humann_regroup` | `true` | Regroup gene families (`humann_regroup_table`) |
| `--humann_regroup_grouping` | `uniref90_rxn` | Target grouping: `uniref90_rxn \| uniref90_ko \| uniref90_eggnog \| uniref90_pfam` |
| `--run_humann_rename` | `true` | Add human-readable names (`humann_rename_table`) |
| `--run_humann_merge` | `true` | Merge per-sample tables into one matrix (`humann_join_tables`) |

### Viral profiling options (BAQLaVa)

| Parameter | Default | Description |
|---|---|---|
| `--run_viral_profiling` | `false` | Enable BAQLaVa |
| `--baqlava_bypass_depletion` | `false` | Bypass bacterial depletion step (needed for test samples with 0 prescreen species) |
| `--baqlava_db` | `null` | Custom BAQLaVa database path; uses bundled db when `null` |

### Strain profiling options (StrainPhlAn)

| Parameter | Default | Description |
|---|---|---|
| `--run_strain_profiling` | `false` | Enable StrainPhlAn |
| `--strainphlan_clades` | `null` | Comma-separated clade list (e.g. `s__Akkermansia_muciniphila`). `null` = all detected. |
| `--strainphlan_phylophlan_mode` | `accurate` | Phylogenetic reconstruction: `accurate \| fast` |
| `--strainphlan_marker_in_n_samples` | `2` | Minimum samples a marker must be present in |
| `--strainphlan_db` | `null` | StrainPhlAn reference db; auto-resolved from `metaphlan_db` when `null` |

### SGB pipeline options

| Parameter | Default | Description |
|---|---|---|
| `--phylophlan_path` | `null` | **Required for `sgb_pipeline`**. Path to PhyloPhlAn database directory. |
| `--megahit_min_contig_length` | `2500` | Minimum contig length for MEGAHIT and MetaBAT2 `-m` |
| `--megahit_options` | `''` | Extra MEGAHIT flags as a quoted string |
| `--metabat_options` | `''` | Extra MetaBAT2 flags |
| `--checkm_predict_options` | `''` | Extra `checkm2 predict` flags |
| `--checkm_coverage_options` | `''` | Extra `checkm coverage` flags |
| `--phylophlan_metagenomic_options` | `''` | Extra `phylophlan_metagenomic` flags |
| `--mash_sketch_options` | `''` | Extra `mash sketch` flags |
| `--sgb_completeness` | `50` | Minimum MAG completeness % for SGB inclusion |
| `--sgb_contamination` | `10` | Maximum MAG contamination % for SGB inclusion |
| `--sgb_abundance_type` | `by_sample` | Abundance estimation: `by_sample \| by_dataset` |
| `--sgb_gc_length_stats` | `false` | Calculate GC content and length statistics per bin |

### Resources & logging

| Parameter | Default | Description |
|---|---|---|
| `--max_memory` | `128.GB` | Global memory cap (used by `check_max()` in `conf/base.config`) |
| `--max_cpus` | `32` | Global CPU cap |
| `--max_time` | `240.h` | Global walltime cap |
| `--log_versions` | `true` | Write tool and database versions to `results/pipeline_info/` |

---

## Compute Profiles

Use `-profile <name>` to select the execution environment.

| Profile | Environment | Executor | Tool delivery |
|---|---|---|---|
| `standard` / `local` | Laptop / workstation | local | conda / system PATH |
| `tufts_hpc` | Tufts HPC | SLURM `batch` | Apptainer containers |
| `harvard_rc` | Harvard FASRC Cannon | SLURM `hsph` | hutlab module system |
| `amazon` | AWS | AWS Batch | ECR containers |
| `engaging` | MIT Engaging | SLURM `newnodes` | system PATH |

> The `harvard_rc` and `tufts_hpc` profiles also auto-load their respective `conf/databases/*.config`,
> so database paths are set automatically — override them on the CLI if needed.

---

## Changing Databases

Every database path is a named parameter. Override on the command line without editing any file:

```sh
nextflow run main.nf -profile harvard_rc \
  --human_genome /n/huttenhower_lab/data/kneaddata_databases/mouse_C57BL \
  --readsdir /path/to/fastqs --outdir results
```

### Available databases on Harvard FASRC

| Tool | Reference | Path |
|---|---|---|
| KneadData | Human hg38 | `/n/huttenhower_lab/data/kneaddata_databases/hg38` |
| KneadData | Mouse C57BL | `/n/huttenhower_lab/data/kneaddata_databases/mouse_C57BL` |
| KneadData | Ribosomal RNA | `/n/huttenhower_lab/data/kneaddata_databases/ribosomal_RNA/SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2` |
| MetaPhlAn | `mpa_vJun23_CHOCOPhlAnSGB_202307` *(default)* | `/n/huttenhower_lab/tools/metaphlan4/rocky8/v4.1.1/lib/python3.10/site-packages/metaphlan/metaphlan_databases` |
| MetaPhlAn | `mpa_vOct22_CHOCOPhlAnSGB_202403` | `/n/huttenhower_lab/tools/metaphlan4/rocky8/v4.0.6_vOct_fixed/lib/python3.10/site-packages/metaphlan/metaphlan_databases` |
| HUMAnN 3.9 | ChocoPhlAn + UniRef90 | `/n/huttenhower_lab/tools/nextflow/databases/humann3` |
| HUMAnN 4 | nucleotide + protein | `/n/huttenhower_lab/tools/nextflow/databases/humann4` |

### Downloading databases

**KneadData:**
```sh
kneaddata_database --download human_genome bowtie2 /path/to/kneaddata_databases
```

**MetaPhlAn** (do not put multiple indexes in the same directory):
```sh
metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202307 \
  --bowtie2db /path/to/metaphlan_databases
```

**HUMAnN:**
```sh
humann_databases --download chocophlan full      /path/to/humann_db/chocophlan
humann_databases --download uniref uniref90_diamond /path/to/humann_db/uniref
humann_databases --download utility_mapping full /path/to/humann_db/utility_mapping
```

---

## Changing Resources

Per-process defaults live in `conf/base.config`. They scale with `task.attempt` (up to 2 retries).
Override for a specific cluster in the profile config, e.g. `conf/profiles/harvard_rc.config`.

**Edit a profile config** (permanent change for that cluster):

```groovy
// conf/profiles/harvard_rc.config
withName: humann {
    memory = '64.G'   // was 32.G
    cpus   = 16       // was 8
    time   = '24.h'   // was 12.h
}
```

**One-off override on the command line** (without editing files):

```sh
nextflow run main.nf -profile harvard_rc \
  --max_memory 64.GB --max_cpus 16 \
  --readsdir /path/to/fastqs --outdir results
```

> `--max_memory` / `--max_cpus` / `--max_time` set global caps applied by `check_max()` in `conf/base.config`.
> They don't set per-process values, but they cap any automatic scaling.

---

## Output Structure

```
results/
├── pipeline_info/
│   ├── versions.txt                               # Tool versions (kneaddata, metaphlan, humann, ...)
│   ├── db_paths.txt                               # Database paths used in this run
│   ├── execution_timeline_YYYY-MM-DD_HH-mm-ss.html
│   ├── execution_report_YYYY-MM-DD_HH-mm-ss.html
│   ├── execution_trace_YYYY-MM-DD_HH-mm-ss.txt
│   └── pipeline_dag_YYYY-MM-DD_HH-mm-ss.svg
│
├── kneaddata/                                     # (run_qc=true)
│   ├── ${SAMPLE}_kneaddata.log
│   ├── ${SAMPLE}_kneaddata_paired_1.fastq.gz      # (paired-end only)
│   ├── ${SAMPLE}_kneaddata_paired_2.fastq.gz
│   ├── ${SAMPLE}_kneaddata_unmatched_1.fastq.gz
│   └── ${SAMPLE}_kneaddata_unmatched_2.fastq.gz
│
├── metaphlan/                                     # (run_taxonomic_profiling=true)
│   ├── ${INDEX}/
│   │   ├── ${SAMPLE}_profile_${INDEX}.tsv         # species relative abundances
│   │   └── ${SAMPLE}_bowtie2_${INDEX}.tsv         # bowtie2 alignment stats
│   ├── bzip/
│   │   └── ${SAMPLE}_${INDEX}.sam.bz2             # compressed SAM (used by StrainPhlAn)
│   └── merged_metaphlan_profiles.tsv              # all samples merged
│
├── humann/                                        # (run_functional_profiling=true)
│   └── ${H_VERSION}/
│       ├── ${SAMPLE}_genefamilies_${H_VERSION}.tsv
│       ├── ${SAMPLE}_pathabundance_${H_VERSION}.tsv
│       ├── ${SAMPLE}_pathcoverage_${H_VERSION}.tsv   # (humann_v37 only)
│       ├── regroup/                               # (run_humann_regroup=true)
│       │   └── ${SAMPLE}_${GROUPING}_${H_VERSION}.tsv
│       ├── rename/                                # (run_humann_rename=true)
│       │   └── ${SAMPLE}_${GROUPING}_named_${H_VERSION}.tsv
│       └── merged/                               # (run_humann_merge=true)
│           ├── merged_genefamilies_${H_VERSION}.tsv
│           ├── merged_pathabundance_${H_VERSION}.tsv
│           └── merged_pathcoverage_${H_VERSION}.tsv
│
├── baqlava/                                       # (run_viral_profiling=true)
│   └── ${SAMPLE}_baqlava/
│
├── strainphlan/                                   # (run_strain_profiling=true)
│   ├── markers/
│   │   └── ${SAMPLE}.pkl
│   └── ${CLADE}/
│       └── output/
│           ├── *.tre                              # phylogenetic tree
│           └── *.tsv                             # marker table
│
└── [sgb_pipeline workflow outputs]
    ├── assembly/main/${SAMPLE}/${SAMPLE}.final.contigs.fa
    ├── assembly/contig_depths/${SAMPLE}.contig_depths.txt
    ├── bins/${SAMPLE}/bins/*.bin.*.fa
    ├── checkm/
    │   ├── ${SAMPLE}/quality_report.tsv
    │   ├── merged_quality_report.tsv
    │   ├── n50/mags_n50.tsv
    │   └── qa/checkm_qa_and_n50.tsv
    ├── phylophlan/
    │   ├── phylophlan_out.tsv
    │   └── phylophlan_relab.tsv
    ├── sgbs/
    │   ├── mash/mash_dist_out.tsv
    │   └── sgbs/SGB_info.tsv
    └── final_profile.tsv                          # merged SGB abundance profile
```

---

## Testing

Test FASTQ files live in `test/`:

| Directory | Contents |
|---|---|
| `test/rawfastq/` | Two paired-end samples |
| `test/single_end_rawfastq/` | `HD32R1_subsample.fastq.gz` — [bioBakery tutorial](https://github.com/biobakery/biobakery_workflows/tree/master/examples/tutorial/input) demo sample |

Run with [nf-test](https://www.nf-test.com/):

```sh
# From HPC
nf-test test tests/main.nf.test         --profile tufts_hpc
nf-test test tests/validate_output.nf.test --profile tufts_hpc

# Locally (requires databases)
nextflow run main.nf -profile standard -params-file template-params.yaml
```

CI runs automatically on every push via `.github/workflows/ci-tests.yml`.

---

## Project Structure

```
biobakery-nextflow/
├── main.nf                              # Router: --workflow mgx|sgb_pipeline|...
├── workflows/
│   ├── mgx.nf                           # MGX: QC → taxonomy → function (+ viral/strain optional)
│   ├── sgb_pipeline.nf                  # MAG assembly → binning → SGB clustering
│   ├── mtx.nf                           # MTX stub
│   ├── mgx_mtx.nf                       # MGX+MTX stub
│   ├── sixteens.nf                      # 16S stub
│   ├── vis.nf                           # Visualization stub
│   └── stats.nf                         # Statistics stub
├── subworkflows/
│   ├── quality_control.nf               # KneadData (single + paired routing)
│   ├── taxonomic_profiling.nf           # MetaPhlAn + bzip + merge
│   ├── functional_profiling.nf          # HUMAnN + regroup + rename + merge
│   ├── viral_profiling.nf               # BAQLaVa
│   └── strain_profiling.nf              # StrainPhlAn (SGB mode)
├── modules/
│   ├── kneaddata/main.nf
│   ├── metaphlan/main.nf                # metaphlan + metaphlan_bzip + metaphlan_merge
│   ├── humann/main.nf
│   ├── strainphlan/main.nf              # sample2markers + strainphlan
│   ├── viral/baqlava/main.nf
│   ├── sgb_pipeline/megahit/main.nf     # MEGAHIT assembly
│   ├── binning/metabat2/main.nf
│   ├── qc/checkm2/main.nf               # checkm2 + merge + n50 + wrangling
│   ├── phylogenomics/phylophlan_metagenomic/main.nf
│   └── utils/
│       ├── align_and_depth/main.nf      # Bowtie2 + jgi_summarize_bam_contig_depths
│       ├── humann_merge/main.nf
│       ├── humann_regroup/main.nf
│       ├── humann_rename/main.nf
│       ├── mash/main.nf                 # sketch + paste + dist + sgb_cluster + merge_tax
│       └── version_log/main.nf
├── conf/
│   ├── base.config                      # Default resources + check_max()
│   ├── profiles/
│   │   ├── harvard_rc.config            # hutlab module system + SLURM hsph
│   │   ├── tufts_hpc.config             # Apptainer + SLURM batch
│   │   ├── aws.config                   # AWS Batch + ECR containers
│   │   └── local.config
│   └── databases/
│       ├── harvard_rc.config            # Default DB paths on Cannon
│       └── tufts.config                 # Default DB paths on Tufts
├── bin/
│   ├── scripts/                         # Python helpers (from anadama2 assembly_tasks/)
│   │   ├── checkm_wrangling.py
│   │   ├── mag_n50_calc.py
│   │   ├── mash_list_inputs.py
│   │   ├── phylophlan_add_tax_assignment.py
│   │   └── ...
│   └── Rscripts/
│       ├── mash_clusters.R
│       └── merge_tax_and_abundance.R
├── assets/
│   ├── diagrams/                        # draw.io source files
│   └── templates/                       # Quarto report templates (planned)
└── test/                                # Test FASTQ files + expected outputs
```
