#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: virus_taxonomy · stage: characterisation · status: skeleton
//
// ICTV taxonomic assignment (MMseqs2 (ICTV)).
//
// Contract (every characterisation component):
//   take  genomes  [ meta, fasta ]  curated viral genomes
//   emit  annotations  [ [id, tool], tsv ]
//
// TODO: MMseqs2 against ICTV references (mmseqs2ictv); see the ICTV
// challenge submissions for alternatives.
// Until then it emits empty channels.
workflow VIRUS_TAXONOMY {

    take:
    genomes

    main:
    skeletonNotice('virus_taxonomy')

    emit:
    annotations = Channel.empty()
}
