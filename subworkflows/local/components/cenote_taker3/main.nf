#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Component: cenote_taker3 · stage: virus_identification · status: skeleton
//
// Hallmark-gene virus discovery (Cenote-Taker 3).
//
// Contract (every virus_identification component):
//   take  contigs  [ [id, assembler], fasta ]
//   emit  viral_contigs  [ [id, assembler, caller], fasta ]
//         scores         [ meta, tsv ]
//
// TODO: Cenote-Taker 3 discovery on the same contigs; tag meta.caller = 'ct'.
// Until then it emits empty channels.
workflow CENOTE_TAKER3 {

    take:
    contigs

    emit:
    viral_contigs = Channel.empty()
    scores        = Channel.empty()
}
