# vis / stats port status — v0.0.4

Working notes for the port of `biobakery_workflows vis` and `stats` (3.2, AnADAMA2)
to Nextflow. Nothing here has been pushed to GitHub and nothing has been deployed
to `lab_storage`. Everything below is local, on branch
`feature/standard-biobakery-workflow`.

Reference implementation being ported:
`/n/lab_storage/huttenhower_lab/tools/biobakery_workflows/rocky8/v3.2/lib/python3.10/site-packages/biobakery_workflows/`
(`workflows/vis.py`, `workflows/stats.py`, `utilities.py`, `visualizations.py`).

---

## Decisions taken

1. **vis/stats take a folder, not channels.** Same as the AnADAMA workflows: point
   them at a bioBakery output folder and they discover the data files by name and
   by content. One code path for both standalone and chained runs.

2. **They run by default at the end of the read-based workflows**, toggled with
   `--run_vis false` / `--run_stats false`.

3. **anadama2 is being removed as a dependency**, not just as a workflow engine.
   This is the decision that shapes the rest:
   - The workflow engine was already gone by construction — every
     `workflow.add_task()` became a Nextflow process.
   - `anadama2.document.PweaveDocument` is *also* the plotting library (the
     templates and `visualizations.py` call ~15 methods on it), so it could not
     simply be dropped. It was **vendored** into `bin/lib/biobakery_document.py`.
   - The `.pmd` templates hard-code `from anadama2 import PweaveDocument` on line
     1, so all 18 were **vendored** into `assets/document_templates/` with that
     import swapped.
   - `biobakery_workflows` itself imports anadama2 at module scope
     (`utilities.py:35`, `files.py:30`), so a small stand-in package lives at
     `bin/lib/anadama2_fallback/anadama2/`. It is put *first* on `sys.path`, so
     the real framework is never imported even where it is installed.

4. **Report format defaults to `html`** (`--report_format pdf` still works).
   AnADAMA defaults to pdf; html avoids needing pdflatex at render time.

5. **Full stats parity** — feature tables, mantel, MaAsLin2 (+ figure tiles),
   HAllA, stratified pathway barplots, beta diversity / PERMANOVA, report, archive.

---

6. **The HUMAnN count tables are excluded from the stats feature set.**
   `get_input_files_for_study_type()` sweeps every `.tsv` in the input folder
   that is not the taxonomy or pathway file into `other_data_files`, so
   `humann_feature_counts.tsv` and `humann_read_and_species_count_table.tsv`
   become feature tables and get run through MaAsLin2, HAllA, the mantel test
   and beta diversity. They are per-sample QC summaries, not abundance tables,
   and they are written samples-as-rows while every consumer expects samples as
   columns — MaAsLin2 and `beta_diversity.R` detect the orientation and survive,
   the mantel test and HAllA do not and fail the whole run.

   Upstream already excludes the kneaddata read count table from stats for the
   same reason (`exclude_types`) and simply missed these two, so upstream 3.2
   cannot complete a wmgx stats run on a folder containing them either.
   `biobakery_identify_inputs.exclude_non_abundance_files()` extends the
   exclusion. Stats therefore analyses taxonomy, pathways and ECs.

   Nothing is lost: those tables are exactly what the *vis* report's HUMAnN
   read-count and feature-count sections are built from.

---

7. **Report artifacts are named for this pipeline's vocabulary.** The vis report
   is `mgx_report.html`, not upstream's `wmgx_report.html`, to match the
   pipeline's own workflow names (`mgx | mtx | mgx_mtx | 16s`). Upstream
   *identifiers* are untouched: the `files.ShotGun` type keys
   (`wmgx_taxonomy`, `wmgx_qc_readcounts`, …) and the packaged
   `wmgx_methods.txt` still carry the upstream spelling, because they are that
   API, not our naming.

---

## What exists

### Nextflow

| File | Contents |
|---|---|
| `workflows/vis.nf` | `VIS` — discovery → taxonomy feature table → alpha diversity → EC names → report → archive |
| `workflows/stats.nf` | `STATS` — discovery → feature tables → mantel / MaAsLin2 / HAllA / stratified barplots / beta diversity or PERMANOVA → report → archive |
| `modules/vis/{identify_inputs,alpha_diversity,add_ec_names,report}/main.nf` | |
| `modules/stats/{feature_table,maaslin2,halla,mantel,beta_diversity,covariate_equation,stratified_pathways,report}/main.nf` | |
| `assets/Rscripts/ggplot2_labs_shim.R` | in-memory patch for MaAsLin2's `ggplot2::labs("")`, sourced by the `maaslin2` module |
| `modules/utils/archive/main.nf` | shared zip archive, ports `workflow.add_archive()` |
| `assets/NO_FILE` | placeholder for optional `path` inputs |

### Python (`bin/scripts/`)

| File | Ports |
|---|---|
| `biobakery_bootstrap.py` | puts the vendored layer on `sys.path` **and** `PYTHONPATH` (pweave runs in a subprocess) |
| `biobakery_identify_inputs.py` | discovery from `vis.py` 76-131 / `stats.py` 85-103 → JSON manifest |
| `biobakery_vis_report.py` | document stage of `vis.py` 94-202 |
| `biobakery_stats_report.py` | document stage of `stats.py` 166-202 |
| `stats_stratified_metadata.py` | metadata prep + `create_merged_data_file()`, emits the barplot plan |
| `stats_humann_barplot.py` | `run_humann_barplot()` (utilities 900-936) |
| `stats_covariate_equation.py` | equation builder from `run_beta_diversity()` (utilities 406-421) |
| `maaslin_image_tiles.py` | `get_maaslin_image_files()` + `generate_tiles_of_maaslin_figures()` |
| `transpose_table.py` | the `transpose()` closure in `run_halla_on_input_file_set()` |

### Vendored (`bin/lib/`)

- `biobakery_document.py` — anadama2 `document.py` (1543 lines) with its only two
  internal imports replaced; `LICENSE-anadama2` alongside it (MIT).
- `biobakery_log.py` — replaces `LoggerReporter.read_log()` for `workflow_info.pmd`.
- `anadama2_fallback/anadama2/` — the stand-in package.

### Vendored (`assets/Rscripts/`)

R scripts are **not** vendored wholesale — most run unmodified and should keep
tracking upstream. Only a patched script lives here, and it shadows the package
copy of the same name via `biobakery_bootstrap.get_rscript()`, which otherwise
falls through to `utilities.get_package_file(name, "Rscript")`. Every module
that runs an R script now goes through that resolver, so patching a second one
later means dropping a file in, not editing a module.

- `alpha_diversity.R` — one line: dropped a no-op `labs("")`. ggplot2 4.x builds
  labels as an S7 object that rejects unnamed labels, so upstream aborts on the
  first continuous metadata variable.
- `beta_diversity.R` — the univariate branch's `adonis()` ported to `adonis2()`,
  which is not a rename: `adonis()` returned a list with the anova table under
  `$aov.tab`, `adonis2()` *is* the table. `adonis` is defunct in vegan 2.7.2.
- `mantel_test.R` — each input table's orientation is now checked against the
  metadata before transposing, the way upstream's own `beta_diversity.R`
  already does. Upstream transposes unconditionally, which empties the distance
  matrix for the samples-as-rows count tables.

---

## Test fixture

The netscratch data this port was previously driven by has been deleted. The
replacement lives outside the repo, under `~/biobakery_vis_stats_test/`:

| File | What it is |
|---|---|
| `make_fixture.py` | seeded generator; writes a bioBakery-standard folder — 16 samples, 23 species, 20 pathways (stratified), 15 ECs, QC and HUMAnN count tables, metadata with 2 categorical + 2 continuous variables, an AnADAMA-format log |
| `env.sh` | the full hand-run environment (see the recipe at the bottom of this file) |
| `input/` | the generated fixture |
| `nf_vis/`, `nf_stats/` | Nextflow run directories |

Half the samples carry a case effect on a few species and pathways, so MaAsLin2,
PERMANOVA and the mantel test have real signal rather than noise. Regenerate
with `python make_fixture.py --output input`.

`test/tutorial_output` in the repo is **not** usable for vis: it is a single
sample and 29 clades, so `document.read_table()` finds no sample set.

## Verified so far

Everything below was re-verified on a clean run (no `-resume`, work directories
deleted) against the current working tree.

- All Python helpers compile; all four `assets/Rscripts/*.R` parse.
- `biobakery_workflows` imports with **no real anadama2 present**; `import anadama2`
  resolves to the stand-in.
- `nextflow config` parses; both DAGs build under `-preview`.
- Input discovery finds all six file roles for vis and for stats.
- **The vis report renders end to end**, by hand and through Nextflow: all 45
  chunks, 0 chunk errors, every report section present (QC, taxonomy,
  ordination, heatmaps, barplots, pathways/ECs, features, software versions,
  tasks run).
- **The whole `VIS` Nextflow workflow completes**: `identify_inputs` →
  `feature_table` → `alpha_diversity` → `add_ec_names` → `vis_report` →
  `archive_output`, publishing `outdir/vis/{mgx_report.html,figures,data,
  alpha_diversity_plots}` and `outdir/vis.zip` — and nothing else, in
  particular no `outdir/stats/`.
- **The vis report is correct, not merely rendered**: 74 figures, every one of
  them a relative link that resolves on disk, no error or traceback text.
  Some figure names contain spaces and are percent-encoded in the HTML; that is
  correct and resolves in a browser.
- Published report links are relative, so the folder relocates as a unit.

- **The whole `STATS` Nextflow workflow completes**: `identify_inputs` →
  `feature_table` (3) → `mantel_test` → `maaslin2` (3) →
  `halla_transpose_metadata` → `halla` (3) → `stratified_metadata` →
  `stratified_barplot` (6) → `covariate_equation` → `beta_diversity` (9) →
  `stats_report` → `archive_output`, publishing
  `outdir/stats/{stats_report.html,features,beta_diversity,mantel_test,
  maaslin2_*,halla_*,stratified_pathways,data}` and `outdir/stats.zip`.
  Requires the HAllA environment override — see issue 2.
- **The stats report is correct, not merely rendered**: 16 figures, every one
  of them relative and resolving on disk, no error or traceback text, and every
  section present (mantel, MaAsLin2 heatmap + tiles, HAllA, stratified
  pathways, beta diversity across all three feature types and all three
  analyses).
- `halla_pathways` and `halla_ec` legitimately produce no `hallagram.png` —
  HAllA only draws one when there are associations to plot, which is why the
  module marks it optional output and the report shows only `halla_taxonomy`.

## Not yet verified

- No parity diff against `biobakery_workflows vis|stats` output yet. The 3.2
  reference output that existed on netscratch was deleted, so this needs an
  upstream run to be redone first. Note that an upstream 3.2 stats run cannot
  currently complete on a standard wmgx folder (see decision 6 and issues 1–2),
  so the comparison will have to be section by section rather than a diff.
- Chained mode (vis/stats at the end of mgx) is **not wired up**.
- `--report_format pdf` has not been exercised since the link changes.

---

## Open issues, in priority order

### ~~1. `beta_diversity.R` uses `adonis`, defunct in vegan 2.7.2~~ — fixed
Vendored and ported to `adonis2(..., by="terms")`, reading row 1 of the returned
anova table in place of `$aov.tab[1,]`. All 12 `STATS:beta_diversity` tasks
(4 feature types × univariate/multivariate/pairwise) now pass, and the
univariate table carries the fixture's seeded case effect (taxonomy R² 46%,
p=0.001).

### ~~2. HAllA cannot load libR through rpy2~~ — fixed, by giving HAllA its own module
The `LD_LIBRARY_PATH` line was the right diagnosis but the wrong scope. Once it
was set, HAllA failed earlier still, in argument parsing:

```
pkg_resources.ContextualVersionConflict: (numpy 2.2.6 (.../assembly_depends/...),
  Requirement.parse('numpy<2,>=1.18'), {'statsmodels'})
```

HAllA pins `numpy<2` through statsmodels and enforces it with
`pkg_resources.require('HAllA')` before it reads `sys.argv`. The
`biobakeryworkflows/3.2` module puts `assembly_depends` — numpy 2.2.6 — ahead of
`depends` on `PYTHONPATH`, deliberately, because the matplotlib in `depends` is
too old for numpy 2 and the report rendering needs the newer one. The two
requirements cannot both be met in one environment.

They do not have to be. There is a standalone HAllA install behind
`rocky8/halla/0.8.20` (`/n/lab_storage/huttenhower_lab/tools/halla/rocky8`) with
numpy 1.26.4, its own `R_LIBS` and R 4.3.3, and it resolves the `libRblas.so`
load too. **The `halla` process must load that module and not
`biobakeryworkflows/3.2`.** Wired up as a `withName: halla` `beforeScript` in
`conf/profiles/harvard_rc.config`; for `-profile local` hand runs the same
environment is set without lmod by
`~/biobakery_vis_stats_test/halla_env.config`, passed with `-c`.

Verified by hand: HAllA runs to completion on `taxonomy_features.txt` under that
environment and writes `hallagram.png` plus a significant association.

### ~~3. Duplicate feature-table names~~ — fixed
`STATS:mantel_test` failed with `input file name collision -- multiple input
files for each of the following file names: counts_features.txt`.

Root cause is upstream: `get_input_files_for_study_type()` types a file by the
last underscore-delimited part of its bioBakery type, so `wmgx_feature_counts`
and `wmgx_humann_counts` both become `counts`. `create_feature_table_inputs()`
then points both at `features/counts_features.txt` and keeps only the last in
`feature_tasks_info` — AnADAMA tolerates two tasks writing one target, Nextflow
does not.

`biobakery_identify_inputs.deduplicate_by_type()` makes that overwrite explicit
(last wins), so the manifest carries one file per type. Both of those files are
now excluded from the stats feature set outright (see the decision below), so
the collision cannot arise; `deduplicate_by_type()` stays as the general
safeguard, since any other pair of types colliding on a suffix would hit it.

### 4. Chained mode is not wired
`VIS`/`STATS` take a `ready` channel to order them after profiling, and `main.nf`
passes `Channel.of(true)` for standalone runs. The read-based workflows do not yet
call them. Two real problems have to be solved first:

- **`publishDir` is asynchronous.** Nextflow does not guarantee published files
  are on disk before a downstream process runs, so pointing vis at `params.outdir`
  mid-run is a race.
- **Layout mismatch.** `files.ShotGun` expects bioBakery names and locations —
  `metaphlan/merged/metaphlan_taxonomic_profiles.tsv`,
  `kneaddata/merged/kneaddata_read_count_table.tsv`,
  `humann/merged/{pathabundance_relab,ecs_relab}.tsv`,
  `humann/counts/humann_{feature_counts,read_and_species_count_table}.tsv`.
  This pipeline publishes different names, so the named lookups miss and only
  content-sniffed files (taxonomy, pathways, ECs) are found. The report loses its
  QC read-count, HUMAnN read-count and feature-count sections.

Both are fixed the same way: a `stage_report_input` process that builds a
bioBakery-standard folder from the actual output *channels*. That is the
recommended next piece of work.

### 5. Nextflow does not re-run on driver changes
`bin/scripts/*.py` and `bin/lib/*.py` are referenced by `${projectDir}` path,
not staged as task inputs, so editing them does not change a task hash and
`-resume` happily serves a stale cached result. Every driver change during this
work needed a full re-run. Either stage the scripts as inputs or remember to
drop `-resume` when a driver changes.

### 6. Not yet done
- `conf/base.config` and `conf/profiles/harvard_rc.config` have **no resource
  blocks or `beforeScript` module loads** for the new processes, with one
  exception: `withName: halla` is wired in `harvard_rc.config`, because that
  one is not a tuning knob but a correctness requirement (issue 2). Everything
  else still needs blocks.
- README / architecture docs not updated.
- The 0.0.4 hutlab modulefile is not written. Template it from
  `/n/lab_storage/huttenhower_lab/tools/hutlab/src/modules_rocky8/rocky8/biobakery-workflows-nextflow/0.0.3`.
  It must provide, beyond the 0.0.3 contents: **pweave** (a separate PyPI
  package that happens to live in the anadama2 module's site-packages here),
  **R 4.5.1 with the anadama2 `R_LIBS`** (that is where `vegan` is — not in the
  R module, not in biobakery's own `R_LIBS`), **`$R_HOME/lib64/R/lib` on
  `LD_LIBRARY_PATH`**, and **HUMAnN** on `PATH` for `add_ec_names`, which
  shells out to `humann_rename_table`. It must **not** try to accommodate
  HAllA: HAllA needs numpy < 2 and gets its own module on the `halla` process
  instead (issue 2).

---

## Known upstream defects worked around

Each of these is 3.2 code running against the current pinned R/Python stack, not
a porting mistake. They are recorded so the divergences stay deliberate.

| Where | Symptom | Handling |
|---|---|---|
| `document.show_pcoa()` continuous branch | `pyplot.colorbar(scalarmappaple)` with a detached ScalarMappable → `Unable to determine Axes to steal space for Colorbar` on matplotlib ≥ 3.6, killing the chunk and losing every continuous-metadata PCoA | fixed in the vendored `biobakery_document.py`: `ax=subplot` |
| `Rscripts/alpha_diversity.R` | `labs("")` → `<ggplot2::labels> object is invalid` on ggplot2 4.x, aborting on the first continuous variable | vendored in `assets/Rscripts/`, no-op call dropped |
| `Rscripts/beta_diversity.R` univariate branch | `adonis` defunct in vegan 2.7.2; the pairwise and multivariate branches of the same script already call `adonis2` | vendored in `assets/Rscripts/`, ported to `adonis2(..., by="terms")` and reading row 1 of the returned table instead of `$aov.tab[1,]` |
| `Rscripts/mantel_test.R` | transposes every feature table on the assumption that samples are its columns. `create_feature_table.py` never transposes, so the count tables (`humann_feature_counts.tsv`, `humann_read_and_species_count_table.tsv`) arrive samples-as-rows, transpose to rows named for the count columns, and intersect no other table — `vegdist()` returns an empty distance matrix and ade4 aborts in `bicenter.wt()` with "weights must be non-negative and not all zero" before any test runs | vendored in `assets/Rscripts/`, with the orientation check upstream's own `beta_diversity.R` already does (which is why beta diversity on the same tables succeeds) |
| `get_input_files_for_study_type()` | two files collapse to type `counts`, then `create_feature_table_inputs()` writes both to one target | `deduplicate_by_type()` in the discovery step |
| `get_input_files_for_study_type()` | sweeps the HUMAnN count tables into the stats feature set, though they are samples-as-rows QC summaries rather than abundance tables; `exclude_types` already drops the kneaddata read count table for the same reason and missed these | `exclude_non_abundance_files()` in the discovery step — see decision 6 |
| `Maaslin2:::maaslin2_association_plots` (MaAsLin2 1.22.0) | `ggplot2::labs("")` → `<ggplot2::labels> object is invalid` on ggplot2 4.x, on the first *continuous* metadata variable, after the results tables are already written | `assets/Rscripts/ggplot2_labs_shim.R`, sourced by the `maaslin2` module: rewrites that one function in memory, since the shared install cannot be patched |
| `hclust2` | `ValueError: Invalid vmin or vmax` on degenerate data (seen on the fixture's `ecs_zscore` heatmap only) | already tolerated: `show_hclust2()` catches it and prints "Unable to generate heatmap" |

The vendored `sh()` **raising** on a non-zero exit is *not* a divergence:
upstream `anadama2/util.py:sh` (line 313-321) raises `ShellException` the same
way. This was previously flagged as uncertain; it is settled.

---

## Port-specific fixes worth remembering

- **Report links are relativized.** The `.pmd` templates embed figures by
  absolute `document.figures_folder` path. Under AnADAMA that folder is the
  final destination; under Nextflow it is a work directory, so a published
  report pointed at paths that vanish on cleanup.
  `biobakery_bootstrap.relativize_report_links()` strips the report folder
  prefix from the rendered HTML (PDF needs nothing — pandoc embeds the images).
- **Everything the report references is staged into the report folder first**,
  so there is a prefix to strip: `stage_alpha_diversity_plots()` copies the
  plots in (which is also the layout upstream writes), and
  `stage_package_image()` copies the `wms_workflow` diagram out of the install
  prefix so the folder is self-contained off-cluster.
- **`alpha_diversity` has no `publishDir`.** It used to publish to
  `${params.outdir}/vis`, the same folder `vis_report` publishes wholesale, and
  the later publish clobbered the plots.
- **`feature_table` / `trim_taxonomy` take a `publish` flag** for the same
  reason. Upstream writes the feature tables into whichever output folder called
  `create_feature_table_inputs()`, so stats gets `stats/features/` and vis gets
  `vis/features/`. The shared module published to a fixed `stats/features/`,
  which put a `stats` folder in the output of every vis-only run. vis now passes
  `false`: the table is only an intermediate for the alpha diversity plots
  there, and publishing it into `vis/` would race `vis_report`'s wholesale copy
  exactly as `alpha_diversity` did.
- **The render scratch directory is cleaned up after the report.**
  `PweaveDocument.create()` renders inside a `tempfile.mkdtemp()` made in the
  report folder and clears it with `rmtree(ignore_errors=True)`; on NFS that
  routinely leaves the emptied directory behind, because deleting a still-open
  file makes the server silly-rename it, the `rmdir` then fails, and the error
  is swallowed. Upstream never notices — there the folder is just a location on
  disk. Here it is a published artifact, so the stray `tmp*` directory shipped
  in `outdir/<workflow>/` and inside the zip.
  `biobakery_bootstrap.remove_render_temp_dirs()` removes it, and only if empty.

---

## Running the report drivers by hand

Use `~/biobakery_vis_stats_test/env.sh`, which is the recipe below kept in one
place:

```bash
source ~/biobakery_vis_stats_test/env.sh
```

The module environment is fiddly. `PYTHONPATH` order matters: `matplotlib` in
`depends/` is too old for numpy 2.x, so `assembly_depends/` must come first.

```bash
T=/n/lab_storage/huttenhower_lab/tools/biobakery_workflows/rocky8
A=/n/lab_storage/huttenhower_lab/tools/anadama2/rocky8/v0.10.0
H=/n/lab_storage/huttenhower_lab/tools/humann4/rocky8
R_HOME_DIR=/n/sw/helmod-rocky8/apps/Core/R/4.5.1-fasrc01
REPO=/n/home04/smaharjan/biobakery-nextflow

export PATH=$R_HOME_DIR/bin:$T/v3.2/bin:$T/assembly_depends/bin:$T/depends/bin:$T/conda/pandoc_env/bin:$A/bin:$H/v4_alpha_1_final/bin:$H/depends:/n/sw/Mambaforge-22.11.1-4/bin:$PATH
export LD_LIBRARY_PATH=$T/depends/lib:$R_HOME_DIR/lib64/R/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$T/v3.2/lib/python3.10/site-packages:$T/assembly_depends/lib/python3.10/site-packages:$T/depends/lib/python3.10/site-packages:$T/depends/lib:$A/lib/python3.10/site-packages:$H/v4_alpha_1_final/lib/python3.10/site-packages
export R_LIBS=$A/../R_LIBS:$T/R_LIBS:$T/assembly_depends/R_LIBS
export MPLBACKEND=Agg
export PATH=$REPO/bin/scripts:$PATH
export PYTHONPATH=$REPO/bin/scripts:$PYTHONPATH
```

Three of these were missing from the earlier version of this recipe and each
cost a failed run:

- **R on `PATH`.** `document.compute_pcoa()` runs `subprocess.Popen(["R", ...])`;
  with no R it raises `FileNotFoundError` and kills the chunk. This was the
  long-standing "PCoA chunk failure".
- **`R_LIBS` including the *anadama2* one.** `vegan` is there, not in the R
  module (base packages only) and not in biobakery's `R_LIBS`.
- **`$R_HOME/lib64/R/lib` on `LD_LIBRARY_PATH`** for HAllA's rpy2.

`$A/lib/python3.10/site-packages` is needed **only for pweave**; the bootstrap
still forces `import anadama2` to the stand-in. HUMAnN's site-packages goes last
so it cannot disturb the numpy/matplotlib resolution order.

For Nextflow itself, JDK 11+ is required (the login shell default is Java 8):

```bash
export JAVA_HOME=/n/sw/helmod-rocky8/apps/Core/jdk/21.0.2-fasrc01
export PATH=$JAVA_HOME/bin:/n/lab_storage/huttenhower_lab/tools/nextflow/24.10.4/bin:$PATH

```

The two runs used to get here, verbatim (drop `-resume` whenever a driver in
`bin/` changed — see issue 5):

```bash
source ~/biobakery_vis_stats_test/env.sh
cd ~/biobakery_vis_stats_test/nf_vis

nextflow run $REPO/main.nf --workflow vis -profile local -w $PWD/work \
  --vis_input   ~/biobakery_vis_stats_test/input \
  --input_metadata ~/biobakery_vis_stats_test/input/metadata.tsv \
  --outdir $PWD/out --project_name "vis nf test"

cd ~/biobakery_vis_stats_test/nf_stats

nextflow run $REPO/main.nf --workflow stats -profile local \
  -c ~/biobakery_vis_stats_test/halla_env.config -w $PWD/work \
  --stats_input ~/biobakery_vis_stats_test/input \
  --input_metadata ~/biobakery_vis_stats_test/input/metadata.tsv \
  --outdir $PWD/out --project_name "stats nf test" \
  --stats_fixed_effects 'diagnosis,age'
```

The `-c ~/biobakery_vis_stats_test/halla_env.config` is **required** for stats:
env.sh cannot host HAllA, which pins numpy < 2, alongside the report rendering,
which needs numpy 2. That file gives the `halla` process its own environment
without lmod, and is the hand-run equivalent of the `withName: halla`
`beforeScript` in `conf/profiles/harvard_rc.config`. See issue 2.

## Next session, in order

1. Chained mode (issue 4) — the `stage_report_input` process. This is the last
   functional gap; both workflows are complete and verified standalone.
2. Configs and the modulefile (issue 6). Only `withName: halla` exists so far.
3. `--report_format pdf`, which has not been exercised since the link changes.
4. Section-by-section comparison against a fresh upstream 3.2 run.
