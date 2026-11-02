#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    biobakery/hvp — Human Virome Project workflows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/biobakery/HVP-nextflow

    Laid out as the nf-core pipeline template is: initialise → named pipeline
    workflow → complete. The pipeline itself is the HVP orchestrator
    (workflows/hvp.nf), which runs components chosen at launch:

      --components kneaddata,baqlava,assembly   the full orchestration (default)
      --only genomad                            one component on its own, plus the
                                                preprocessing and whatever it needs

    Component names, stages and status are in conf/components.json, and drawn
    in docs/diagrams/architecture.mmd.
----------------------------------------------------------------------------------------
*/

include { HVP                     } from './workflows/hvp.nf'
include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION; resolveComponents } from './subworkflows/local/utils_hvp_pipeline/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Named, so the run can be imported into a larger pipeline the way nf-core
// pipelines are.
workflow BIOBAKERY_HVP {

    take:
    reads       // channel: [ [id: sample, paired_end: bool], reads ]
    components  // List<String>: the components to run

    main:
    HVP(reads, components)

    emit:
    read_profiles = HVP.out.read_profiles
    contigs       = HVP.out.contigs
    viral_genomes = HVP.out.viral_genomes
    annotations   = HVP.out.annotations
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.readsdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    // --only runs a single component; otherwise --components decides.
    def components = resolveComponents(params.components, params.only)

    BIOBAKERY_HVP (
        PIPELINE_INITIALISATION.out.reads,
        components
    )

    PIPELINE_COMPLETION (
        params.monochrome_logs
    )
}
