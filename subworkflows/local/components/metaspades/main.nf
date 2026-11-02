#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: metaspades · stage: assembly · status: skeleton
//
// Alternative short-read assembler (metaSPAdes).
//
// Contract (every assembly component):
//   take  reads  [ [id, paired_end], fastq ]  clean reads, one file per sample
//         parts  [ id, [files] ]              paired samples: the separate mates
//   emit  contigs  [ [id, assembler], fasta ]
//
// TODO: metaSPAdes (seems a good default assembler); the nf-core spades module
// may be reusable. Emit the assembly graph too if variants are wanted later.
// Until then it emits empty channels.
workflow METASPADES {

    take:
    reads
    parts

    main:
    skeletonNotice('metaspades')

    emit:
    contigs = Channel.empty()
}
