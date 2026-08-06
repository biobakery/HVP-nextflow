#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { single_end_kneaddata } from '../modules/kneaddata/main.nf'
include { paired_end_kneaddata } from '../modules/kneaddata/main.nf'

// Quality control subworkflow: host decontamination + trimming via KneadData
workflow QUALITY_CONTROL {

    take:
    reads  // Channel: [ [id: sample, paired_end: bool], reads ]

    main:
    // Split channel by library type
    paired_reads = reads
        .filter  { meta, r -> meta.paired_end }
        .map     { meta, r -> tuple(meta.id, r) }

    single_reads = reads
        .filter  { meta, r -> !meta.paired_end }
        .map     { meta, r -> tuple(meta.id, r) }

    paired_out = paired_end_kneaddata(paired_reads)
    single_out = single_end_kneaddata(single_reads)

    // Re-attach meta map so downstream subworkflows get consistent channel shape
    cleaned_paired = paired_out.kneads.map { sample, reads -> [ [id: sample, paired_end: true],  reads ] }
    cleaned_single = single_out.kneads.map { sample, reads -> [ [id: sample, paired_end: false], reads ] }

    all_cleaned = cleaned_paired.mix(cleaned_single)
    all_logs    = paired_out.log.mix(single_out.log)

    emit:
    reads = all_cleaned  // Channel: [ [id, paired_end], cleaned_reads ]
    // named "logs", not "log": Nextflow already binds "log" to its own logger,
    // so emitting that name fails at runtime with "No such variable: log"
    logs  = all_logs
}
