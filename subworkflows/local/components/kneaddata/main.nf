#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { single_end_kneaddata } from '../../../../modules/local/kneaddata/main.nf'
include { paired_end_kneaddata } from '../../../../modules/local/kneaddata/main.nf'
include { kneaddata_read_counts } from '../../../../modules/local/kneaddata/main.nf'

// Component: kneaddata · stage: preprocess · status: ready
//
// Trim reads and remove host reads with KneadData.
//
// Contract (every preprocess component):
//   take  reads  [ [id, paired_end], reads ]  one file, or [R1, R2]
//   emit  reads  [ [id, paired_end], fastq ]  one file per sample
//         parts  [ id, [paired_1, paired_2, unmatched_1, unmatched_2] ]  paired
//                samples only; what an assembler wants instead of one file
workflow KNEADDATA {

    take:
    reads    // Channel: [ [id: sample, paired_end: bool], reads ]

    main:
    def subdir   = ''
    def db_paths = [params.host_genome].findAll { it }
    if (!db_paths) {
        error "ERROR: no KneadData reference database configured. Set --host_genome, " +
              "or leave 'kneaddata' out of --components to skip read QC."
    }
    // Pre-formed argument string: see the note in modules/local/kneaddata/main.nf on
    // why a List must not be handed to a process val input.
    def db_args = db_paths.collect { "--reference-db ${it}" }.join(' ')

    // Split channel by library type
    paired_reads = reads
        .filter  { meta, r -> meta.paired_end }
        .map     { meta, r -> tuple(meta.id, r) }

    single_reads = reads
        .filter  { meta, r -> !meta.paired_end }
        .map     { meta, r -> tuple(meta.id, r) }

    paired_out = paired_end_kneaddata(paired_reads, db_args, subdir)
    single_out = single_end_kneaddata(single_reads, db_args, subdir)

    // Re-attach meta map so downstream subworkflows get consistent channel shape
    cleaned_paired = paired_out.kneads.map { sample, reads -> [ [id: sample, paired_end: true],  reads ] }
    cleaned_single = single_out.kneads.map { sample, reads -> [ [id: sample, paired_end: false], reads ] }

    all_cleaned = cleaned_paired.mix(cleaned_single)
    all_logs    = paired_out.log.mix(single_out.log)

    // One read count table for the run: reads in, trimmed, and left after host removal.
    read_counts_out = kneaddata_read_counts(all_logs.map { sample, l -> l }.collect(), subdir)

    // The four separate KneadData outputs, for consumers that need the pairing
    // rather than the concatenated file: the assembly component runs MEGAHIT in
    // paired mode (-1/-2 plus the unmatched reads as single-end).
    paired_parts = paired_out.paired1
        .join(paired_out.paired2)
        .join(paired_out.unpaired1)
        .join(paired_out.unpaired2)
        .map { sample, p1, p2, u1, u2 -> tuple(sample, [p1, p2, u1, u2]) }

    emit:
    reads       = all_cleaned          // Channel: [ [id, paired_end], cleaned_reads ]
    parts       = paired_parts         // Channel: [ sample, [p1, p2, u1, u2] ]
    read_counts = read_counts_out.counts
}
