#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { skeletonNotice } from '../../utils_hvp_pipeline/main.nf'

// Component: circular_contigs · stage: virus_curation · status: skeleton
//
// Trim DTRs and reorient circular genomes (tr-trimmer, reorient-circular-seq).
//
// Contract (every virus_curation component):
//   take  in_genomes  [ meta, fasta ]  viral sequences from the step before
//   emit  genomes  [ meta, fasta ]
//         report   [ meta, tsv ]
//
// TODO: Detect direct terminal repeats (min length = the longest assembly k-mer)
// on contigs with no host region, keep one copy (tr-trimmer), then reorient
// so no predicted CDS is split (reorient-circular-seq; Dnaapler needs a TerL).
// Until then it passes its input through unchanged, so the stage still works.
workflow CIRCULAR_CONTIGS {

    take:
    in_genomes

    main:
    skeletonNotice('circular_contigs')

    emit:
    genomes = in_genomes
    report  = Channel.empty()
}
