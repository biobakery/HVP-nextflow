#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { TAXONOMIC_PROFILING  } from '../../taxonomic_profiling/main.nf'
include { FUNCTIONAL_PROFILING } from '../../functional_profiling/main.nf'
include { baqlava              } from '../../../../modules/local/baqlava/main.nf'

// Component: baqlava · stage: read_profiling · status: ready
//
// Viral profiling with BAQLaVa. MetaPhlAn runs first because BAQLaVa uses its
// profile to deplete bacterial reads; HUMAnN runs alongside on the same reads
// (turn it off with --baqlava_humann false).
//
// Contract (every read_profiling component):
//   take  reads     [ [id, paired_end], fastq ]  clean reads, one file per sample
//   emit  profiles  [ [id, tool], file ]         one profile per sample
workflow BAQLAVA {

    take:
    reads    // Channel: [ [id: sample, paired_end: bool], fastq ]

    main:
    TAXONOMIC_PROFILING(reads, '')

    if (params.baqlava_humann) {
        FUNCTIONAL_PROFILING(reads, TAXONOMIC_PROFILING.out.profile, '')
    }

    // Join each sample's reads with its own MetaPhlAn profile
    joined = reads
        .map  { meta, r -> tuple(meta.id, r) }
        .join (TAXONOMIC_PROFILING.out.profile)

    baqlava_out = baqlava(
        joined.map { sample, r, profile -> tuple(sample, r) },
        joined.map { sample, r, profile -> tuple(sample, profile) }
    )

    emit:
    profiles  = baqlava_out.results.map { sample, dir -> [ [id: sample, tool: 'baqlava'], dir ] }
    metaphlan = TAXONOMIC_PROFILING.out.merged
}
