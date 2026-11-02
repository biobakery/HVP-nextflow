#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: checkv_completeness · stage: characterisation · status: skeleton
//
// Genome completeness estimates (CheckV).
//
// Contract (every characterisation component):
//   take  genomes  [ meta, fasta ]  curated viral genomes
//   emit  annotations  [ [id, tool], tsv ]
//
// TODO: CheckV completeness, then summarise it together with the
// virus_rename table into one completeness estimate per genome.
// Until then it emits empty channels.
workflow CHECKV_COMPLETENESS {

    take:
    genomes

    main:
    skeletonNotice('checkv_completeness')

    emit:
    annotations = Channel.empty()
}
