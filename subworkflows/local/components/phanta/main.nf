#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: phanta · stage: read_profiling · status: skeleton
//
// Virus-inclusive k-mer profiling (Phanta).
//
// Contract (every read_profiling component):
//   take  reads  [ [id, paired_end], fastq ]  clean reads, one file per sample
//   emit  profiles  [ [id, tool], file ]
//
// TODO: Port the Phanta implementation on the jw_phanta branch (workflows/phanta.nf)
// into this component.
// Until then it emits empty channels.
workflow PHANTA {

    take:
    reads

    main:
    skeletonNotice('phanta')

    emit:
    profiles = Channel.empty()
}
