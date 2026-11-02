#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: sylph · stage: read_profiling · status: skeleton
//
// Containment-based profiling (Sylph).
//
// Contract (every read_profiling component):
//   take  reads  [ [id, paired_end], fastq ]  clean reads, one file per sample
//   emit  profiles  [ [id, tool], file ]
//
// TODO: Sylph profile against a virus-inclusive database. Not in nf-core yet.
// Until then it emits empty channels.
workflow SYLPH {

    take:
    reads

    main:
    skeletonNotice('sylph')

    emit:
    profiles = Channel.empty()
}
