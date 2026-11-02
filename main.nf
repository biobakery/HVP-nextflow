#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ═════════════════════════════════════════════════════════════════════════════
// biobakery/hvp — Human Virome Project workflows
//
// Laid out the way nf-core pipelines are: initialise → named pipeline workflow
// → complete. The pipeline itself is the HVP orchestrator (workflows/hvp.nf),
// which runs components chosen at launch:
//
//   --components kneaddata,baqlava,assembly   the full orchestration (default)
//   --plugin genomad                          one component on its own, plus the
//                                             preprocessing and whatever it needs
//
// Component names, stages and status are listed in conf/components.json, and
// drawn in docs/diagrams/architecture.mmd.
// ═════════════════════════════════════════════════════════════════════════════

include { HVP                     } from './workflows/hvp.nf'
include { PIPELINE_INITIALISATION; PIPELINE_COMPLETION; resolveComponents } from './subworkflows/local/utils_hvp_pipeline/main.nf'

// The named pipeline workflow, so the run can be imported into a larger
// pipeline the same way nf-core pipelines are.
workflow BIOBAKERY_HVP {

    take:
    reads       // Channel: [ [id: sample, paired_end: bool], reads ]
    components  // List<String>: the components to run

    main:
    HVP(reads, components)

    emit:
    read_profiles = HVP.out.read_profiles
    contigs       = HVP.out.contigs
    viral_genomes = HVP.out.viral_genomes
    annotations   = HVP.out.annotations
}

workflow {

    PIPELINE_INITIALISATION(params.validate_params, params.readsdir)

    // --plugin runs a single component; otherwise --components decides.
    def components = resolveComponents(params.components, params.plugin)

    BIOBAKERY_HVP(PIPELINE_INITIALISATION.out.reads, components)

    PIPELINE_COMPLETION()
}
