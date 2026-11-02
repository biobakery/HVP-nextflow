#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: virus_rename · stage: virus_curation · status: skeleton
//
// Stable IDs and a per-genome summary table (to be designed).
//
// Contract (every virus_curation component):
//   take  in_genomes  [ meta, fasta ]  viral sequences from the step before
//   emit  genomes  [ meta, fasta ]
//         report   [ meta, tsv ]
//
// TODO: Rename sequences as <sample>_<caller>_<provirus?>_<circular?>_<n> and write
// a table mapping each new id to the original contig, coordinates, DTR length,
// caller and score. To be designed.
// Until then it passes its input through unchanged, so the stage still works.
workflow VIRUS_RENAME {

    take:
    in_genomes

    main:
    skeletonNotice('virus_rename')

    emit:
    genomes = in_genomes
    report  = Channel.empty()
}
