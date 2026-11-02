#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: genomad · stage: virus_identification · status: skeleton
//
// Viral contigs and proviruses (geNomad).
//
// Contract (every virus_identification component):
//   take  contigs  [ [id, assembler], fasta ]
//   emit  viral_contigs  [ [id, assembler, caller], fasta ]
//         scores         [ meta, tsv ]
//
// TODO: geNomad end-to-end: full contigs and provirus regions (fasta) plus
// virus scores (tsv). Contribute to the existing nf-core genomad module.
// Until then it emits empty channels.
workflow GENOMAD {

    take:
    contigs

    main:
    skeletonNotice('genomad')

    emit:
    viral_contigs = Channel.empty()
    scores        = Channel.empty()
}
