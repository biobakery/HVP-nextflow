#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Component: host_scrubber · stage: preprocess · status: skeleton
//
// Scrub residual human reads (sra-human-scrubber).
//
// Contract (every preprocess component):
//   take  in_reads  [ [id, paired_end], fastq ]  one file per sample
//         in_parts  [ id, [files] ]              paired samples: the separate mates
//   emit  reads  [ [id, paired_end], fastq ]
//         parts  [ id, [files] ]
//
// TODO: Run sra-human-scrubber on the KneadData output and emit the scrubbed reads,
// plus the removed reads as a chaff file. nf-core's detaxizer is an alternative.
// Until then it passes its input through unchanged, so the stage still works.
workflow HOST_SCRUBBER {

    take:
    in_reads
    in_parts

    emit:
    reads = in_reads
    parts = in_parts
}
