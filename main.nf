#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ── Workflow imports ───────────────────────────────────────────────────────
include { MGX }          from './workflows/mgx.nf'
include { MTX }          from './workflows/mtx.nf'
include { MGX_MTX }      from './workflows/mgx_mtx.nf'
include { SIXTEENS }     from './workflows/sixteens.nf'
include { VIS }          from './workflows/vis.nf'
include { STATS }        from './workflows/stats.nf'
include { ASSEMBLY } from './workflows/assembly.nf'
include { PHANTA }       from './workflows/phanta.nf'

// ── Router ─────────────────────────────────────────────────────────────────
workflow {

    // Validate required params for read-based workflows. mgx_mtx is the
    // exception: it takes two input folders instead of one, and validates them
    // itself, mirroring wmgx_wmtx.py replacing --input with --input-metagenome
    // and --input-metatranscriptome.
    if (!params.readsdir && params.workflow in ['mgx', 'mtx', '16s', 'assembly']) {
        error "ERROR: --readsdir is required. Example: --readsdir /path/to/fastqs"
    }

    // phanta is absent from that list on purpose: it is the one workflow that
    // takes a samplesheet rather than a reads folder, and it validates its own
    // input. See phanta_samplesheet_rows in workflows/phanta.nf.
    if (params.workflow == 'phanta') {
        phanta_samplesheet_guard(params.phanta_samplesheet)
    }

    switch (params.workflow) {

        case 'mgx':
            MGX()
            break

        case 'mtx':
            MTX()
            break

        case 'mgx_mtx':
            MGX_MTX()
            break

        case '16s':
            SIXTEENS()
            break

        case 'vis':
            VIS(Channel.value(report_input_dir(params.vis_input, 'vis')))
            break

        case 'stats':
            STATS(Channel.value(report_input_dir(params.stats_input, 'stats')))
            break

        case 'assembly':
            ASSEMBLY()
            break

        case 'phanta':
            PHANTA()
            break

        default:
            error "Unknown workflow '${params.workflow}'. Choose: mgx | mtx | mgx_mtx | 16s | vis | stats | assembly | phanta"
    }
}

// The folder a standalone vis or stats run reports on. Chained runs get the
// folder stage_report_input builds instead, so this is only for --workflow
// vis|stats.
def report_input_dir(setting, label) {
    def dir = file(setting ?: params.outdir)

    if (!dir.exists())
        error "ERROR: ${label} input folder not found: ${dir}\n" +
              "Set --${label}_input to a bioBakery output folder."

    return dir
}


// Fail early and clearly when --workflow phanta is run without a samplesheet,
// or with one that is not there. The contents are parsed and validated in
// workflows/phanta.nf; this only covers the case worth catching before the
// workflow is even entered.
//
// Checks existence and nothing more. `-preview` builds the DAG without running
// a task, and bin/make_diagrams.py previews every workflow in CI against a stub
// samplesheet pointing at empty files, so anything stricter here would turn
// into a CI failure rather than a useful error.
def phanta_samplesheet_guard(setting) {
    if (!setting) {
        error "ERROR: --phanta_samplesheet is required for --workflow phanta.\n" +
              "Example: --phanta_samplesheet /path/to/samples.csv"
    }

    def sheet = file(setting)
    if (!sheet.exists()) {
        error "ERROR: samplesheet not found: ${sheet}"
    }

    return sheet
}
