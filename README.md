# HVP-nextflow

> Human Virome Project workflows in Nextflow: an orchestrator over pluggable
> read-based and assembly-based virome components, built on
> [biobakery-nextflow](https://github.com/biobakery/biobakery-nextflow).

This repository carries only what HVP needs from biobakery-nextflow: KneadData,
MetaPhlAn, HUMAnN, BAQLaVa and the assembly pipeline. On top of those it adds an
orchestrator that runs **components** stage by stage. Tools that are planned but
not yet implemented sit in the pipeline as placeholders, each with its slot,
contract and TODO, so whoever implements one only touches that component.

## Architecture

<!-- generated:architecture -->
<!-- written by bin/make_diagrams.py from conf/components.json; do not edit by hand -->
```mermaid
%%{init: {"flowchart": {"curve": "basis", "nodeSpacing": 24, "rankSpacing": 44, "padding": 10}}}%%
flowchart TB
    entry(["<b>main.nf</b> → orchestrator<br/><i>--components</i> · <i>--only</i>"])
    subgraph st_preprocess["<b>1 · Read QC & host removal</b> &nbsp;<i>in order</i>"]
        direction LR
        kneaddata["<b>kneaddata</b><br/><i>KneadData</i><br/>▣ local modules"]
        host_scrubber["<b>host_scrubber</b><br/><i>sra-human-scrubber</i>"]
        kneaddata --> host_scrubber
    end
    subgraph st_read_profiling["<b>2 · Read-based profiling</b> &nbsp;<i>side by side</i>"]
        direction LR
        baqlava["<b>baqlava</b><br/><i>MetaPhlAn, HUMAnN, BAQLaVa</i><br/>▣ local modules"]
        phanta["<b>phanta</b><br/><i>Phanta</i>"]
        sylph["<b>sylph</b><br/><i>Sylph</i>"]
    end
    subgraph st_assembly["<b>3 · De novo assembly</b> &nbsp;<i>side by side</i>"]
        direction LR
        assembly["<b>assembly</b><br/><i>MEGAHIT, MetaBAT2, CheckM2, PhyloPhlAn, Mash</i><br/>▣ local modules"]
        metaspades["<b>metaspades</b><br/><i>metaSPAdes</i>"]
    end
    subgraph st_virus_identification["<b>4 · Virus identification</b> &nbsp;<i>side by side</i>"]
        direction LR
        genomad["<b>genomad</b><br/><i>geNomad</i>"]
        cenote_taker3["<b>cenote_taker3</b><br/><i>Cenote-Taker 3</i>"]
    end
    subgraph st_virus_curation["<b>5 · Decontamination & curation</b> &nbsp;<i>in order</i>"]
        direction LR
        checkv_contamination["<b>checkv_contamination</b><br/><i>CheckV</i>"]
        circular_contigs["<b>circular_contigs</b><br/><i>tr-trimmer, reorient-circular-seq</i>"]
        virus_rename["<b>virus_rename</b><br/><i>to be designed</i>"]
        checkv_contamination --> circular_contigs
        circular_contigs --> virus_rename
    end
    subgraph st_characterisation["<b>6 · Completeness & taxonomy</b> &nbsp;<i>side by side</i>"]
        direction LR
        checkv_completeness["<b>checkv_completeness</b><br/><i>CheckV</i>"]
        virus_taxonomy["<b>virus_taxonomy</b><br/><i>MMseqs2 (ICTV)</i>"]
    end
    entry -- "raw reads (fastq)" --> st_preprocess
    st_preprocess -- "clean reads (fastq)" --> st_read_profiling
    st_preprocess -- "clean reads (fastq)" --> st_assembly
    st_assembly -- "contigs (fasta)" --> st_virus_identification
    st_virus_identification -- "viral contigs (fasta + scores)" --> st_virus_curation
    st_virus_curation -- "curated viral genomes (fasta + summary)" --> st_characterisation
    st_read_profiling --> out_read_profiles[("read-based profiles (tsv)")]
    st_characterisation --> out_annotations[("completeness + taxonomy (tsv)")]
    subgraph legend[" "]
        direction LR
        l_module["<b>module component</b><br/>▣ runs local or nf-core modules"]
        l_plugin{{"<b>plugin component</b><br/>⧉ built on a Nextflow plugin"}}
        l_skel["<b>placeholder</b><br/>wired in, not implemented"]
    end

    classDef ready fill:#dbeafe,stroke:#2563eb,stroke-width:1.5px,color:#0f172a
    classDef skeleton fill:#ffffff,stroke:#94a3b8,stroke-width:1px,stroke-dasharray:4 3,color:#64748b
    classDef entry fill:#e2e8f0,stroke:#475569,stroke-width:1.5px,color:#0f172a
    classDef output fill:#dcfce7,stroke:#16a34a,stroke-width:1.5px,color:#14532d
    class kneaddata,baqlava,assembly,l_module,l_plugin ready
    class host_scrubber,phanta,sylph,metaspades,genomad,cenote_taker3,checkv_contamination,circular_contigs,virus_rename,checkv_completeness,virus_taxonomy,l_skel skeleton
    class entry entry
    class out_read_profiles,out_annotations output
    style st_preprocess fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style st_read_profiling fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style st_assembly fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style st_virus_identification fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style st_virus_curation fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style st_characterisation fill:#f8fafc,stroke:#cbd5e1,stroke-width:1px,color:#334155
    style legend fill:none,stroke:none
    linkStyle default stroke:#64748b,stroke-width:1.5px
```
<!-- /generated:architecture -->

Blue boxes are ready, and dashed boxes are placeholders. Rectangles (`▣`) are
module components, built from local or nf-core modules. Hexagons (`⧉`) are
plugin components, built on a Nextflow plugin. The diagram is generated
from `conf/components.json` and the code by `bin/make_diagrams.py`. CI fails if
it is stale. [docs/architecture.md](docs/architecture.md) explains the layers,
the stage contracts and how to add a component.

## Components

<!-- generated:components -->
<!-- written by bin/make_diagrams.py from conf/components.json; do not edit by hand -->
| Stage | Component | Status | Tools | Implemented as | What it does |
|---|---|---|---|---|---|
| 1 · Read QC & host removal | `kneaddata` | ready | KneadData | ▣ local modules: `kneaddata` | Trim reads and remove host reads |
| 1 · Read QC & host removal | `host_scrubber` | placeholder | sra-human-scrubber | — | Scrub residual human reads |
| 2 · Read-based profiling | `baqlava` | ready | MetaPhlAn, HUMAnN, BAQLaVa | ▣ local modules: `baqlava`, `humann`, `humann_merge`, `humann_regroup` +3 | Bacterial depletion, then viral profiling |
| 2 · Read-based profiling | `phanta` | placeholder | Phanta | — | Virus-inclusive k-mer profiling |
| 2 · Read-based profiling | `sylph` | placeholder | Sylph | — | Containment-based profiling |
| 3 · De novo assembly | `assembly` | ready | MEGAHIT, MetaBAT2, CheckM2, PhyloPhlAn, Mash | ▣ local modules: `abundance`, `align_and_depth`, `checkm2`, `mash` +3 | Contigs, plus MAG binning and SGB clustering |
| 3 · De novo assembly | `metaspades` | placeholder | metaSPAdes | — | Alternative short-read assembler |
| 4 · Virus identification | `genomad` | placeholder | geNomad | — | Viral contigs and proviruses |
| 4 · Virus identification | `cenote_taker3` | placeholder | Cenote-Taker 3 | — | Hallmark-gene virus discovery |
| 5 · Decontamination & curation | `checkv_contamination` | placeholder | CheckV | — | Trim host regions from proviruses |
| 5 · Decontamination & curation | `circular_contigs` | placeholder | tr-trimmer, reorient-circular-seq | — | Trim DTRs and reorient circular genomes |
| 5 · Decontamination & curation | `virus_rename` | placeholder | to be designed | — | Stable IDs and a per-genome summary table |
| 6 · Completeness & taxonomy | `checkv_completeness` | placeholder | CheckV | — | Genome completeness estimates |
| 6 · Completeness & taxonomy | `virus_taxonomy` | placeholder | MMseqs2 (ICTV) | — | ICTV taxonomic assignment |
<!-- /generated:components -->

## Quick start

Harvard FASRC, everything that is ready (KneadData → BAQLaVa, and assembly).
The pipeline needs Nextflow ≥ 25.10.4, the current nf-core standard. The
hutlab launcher is 24.10.4, so set `NXF_VER`; it downloads 25.10.4 on first use:

```sh
source /n/lab_storage/huttenhower_lab/tools/hutlab/src/hutlabrc_rocky8.sh
hutlab load rocky8/biobakery-workflows-nextflow/0.0.4    # Nextflow launcher + JDK
export NXF_VER=25.10.4

nextflow run main.nf -profile harvard_rc \
  --readsdir /path/to/fastqs \
  --outdir   /path/to/output \
  -w         /path/to/output/work
```

Keep the work directory on shared storage: compute nodes cannot see a login
node's `/tmp`.

### Choosing components

```sh
# Read-based profiling only
--components kneaddata,baqlava

# Reads already clean: skip KneadData (pairs are concatenated for the profilers)
--components baqlava,assembly

# BAQLaVa without HUMAnN
--baqlava_humann false

# One component on its own. genomad needs contigs, so the
# assembly component is added automatically. genomad is still a placeholder.
--only genomad
```

### Demo (bundled single-end sample)

```sh
nextflow run main.nf -profile harvard_rc \
  --readsdir tests/data/single_end_rawfastq --single_end true \
  --components kneaddata,baqlava --baqlava_bypass_depletion true \
  --outdir tests/results/demo
```

### Help and parameters

```sh
nextflow run main.nf --help                 # every parameter, from nextflow_schema.json
nextflow run main.nf --help components      # one parameter in detail
nextflow run main.nf --version
```

Parameters are validated against `nextflow_schema.json` at launch by the
[nf-schema](https://nextflow-io.github.io/nf-schema/) plugin, through the
standard nf-core utility subworkflows. Turn that off with
`--validate_params false`. Defaults live in `conf/params.config`.

## nf-core and Nextflow plugins

- **nf-core library.** `subworkflows/nf-core/` holds `utils_nextflow_pipeline`,
  `utils_nfcore_pipeline` and `utils_nfschema_plugin`, installed from
  [nf-core/modules](https://github.com/nf-core/modules) and pinned in
  `modules.json`. Update them with `nf-core subworkflows update`. nf-core tool
  modules go in `modules/nf-core/` (`nf-core modules install <tool>`).
- **Nextflow plugins.** Every plugin the pipeline loads is declared in the
  `plugins { }` block of `nextflow.config`: nf-schema, plus any plugin a
  component is built on. A component lists the plugins it needs in
  `conf/components.json`, and CI checks they're declared. For a one-off plugin,
  use `-plugins nf-<name>@<version>`.

See [docs/architecture.md](docs/architecture.md#adding-a-component) for the
three ways to implement a component: nf-core module, local module or Nextflow
plugin.

## Layout

```
main.nf                                  entry: init → BIOBAKERY_HVP → completion
workflows/hvp.nf                         the orchestrator (six stages)
subworkflows/nf-core/                    nf-core utility subworkflows (pinned in modules.json)
subworkflows/local/
  utils_hvp_pipeline/                    init/completion on the nf-core utils; input channel; component registry
  components/<name>/main.nf              one component each: ready or placeholder
  read_input/ taxonomic_profiling/ functional_profiling/   shared building blocks
modules/local/<tool>/main.nf             one process per tool invocation
conf/components.json                     the component registry (stages, status, plugins)
modules.json                             nf-core modules/subworkflows and their pinned commits
conf/params.config                       every parameter's default
conf/base.config                         default resources per process
conf/profiles/, conf/databases/          per-site runtime and database paths
nextflow_schema.json                     parameter schema (nf-schema)
bin/make_diagrams.py                     generates the diagram and tables above
```

## Configuration

Settings are layered. Later layers win:

1. `conf/params.config`: every parameter's default
2. `conf/base.config`: default resources, capped by `--max_memory/cpus/time`
3. `conf/profiles/<site>.config`: executor, queue and tool modules (`-profile`)
4. `conf/databases/<site>.config`: database paths for that site (`-profile`)
5. `-params-file my-run.yaml`: this run (`conf/harvard_rc.yaml` is a template)
6. `--param value` on the command line

Profiles: `harvard_rc` (FASRC, SLURM + hutlab modules), `tufts_hpc` (SLURM +
Apptainer), `local`, and `amazon` (incomplete). Where each tool's software comes
from is in [docs/dependencies.md](docs/dependencies.md).

## Output

```
<outdir>/
  kneaddata/        clean reads, read count table
  metaphlan/        per-sample and merged taxonomic profiles
  humann/           gene families, pathways, ECs, merged tables
  baqlava/          per-sample viral profiles
  assembly/         contigs and contig depths
  bins/ checkm/ phylophlan/ sgbs/ abundance_*/   MAGs, quality, placement, SGB clusters
  pipeline_info/    pipeline_log.txt (versions, databases, components), timeline, report, trace
```

## Testing

```sh
bin/make_diagrams.py --check         # registry, code and docs agree (Python only)
bin/check_profile_resources.py       # no profile drops a resource request
sbatch tests/submit_tests.sh         # integration suite on FASRC (tools + databases)
```

CI runs the first two plus `nextflow run -preview` for the default, read-only,
no-QC and `--only` selections. That builds the whole graph without running a task.

## History

Forked from biobakery-nextflow 0.0.4; see [HISTORY.md](HISTORY.md).
