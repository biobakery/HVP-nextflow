#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: checkv_contamination · stage: virus_curation · status: skeleton
//
// Trim host regions from proviruses (CheckV).
//
// Contract (every virus_curation component):
//   take  in_genomes  [ meta, fasta ]  viral sequences from the step before
//   emit  genomes  [ meta, fasta ]
//         report   [ meta, tsv ]
//
// TODO: CheckV contamination: trim host regions flanking proviruses, and
// report what was removed (tsv).
// Until then it passes its input through unchanged, so the stage still works.
workflow CHECKV_CONTAMINATION {

    take:
    in_genomes

    main:
    skeletonNotice('checkv_contamination')

    emit:
    genomes = in_genomes
    report  = Channel.empty()
}
