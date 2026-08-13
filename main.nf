#!/usr/bin/env nextflow

nextflow.enable.dsl=2

nextflow.enable.dsl = 2

include { single_end_kneaddata; paired_end_kneaddata } from "${projectDir}/processes/kneaddata.nf"
include { metaphlan; metaphlan_bzip } from "${projectDir}/processes/metaphlan.nf"
include { humann } from "${projectDir}/processes/humann.nf"

params.input_type    = params.input_type ?: 'local'   // local | sra (reads ingestion channel)
params.paired_lookup_pattern  = params.paired_lookup_pattern ?: null
params.single_lookup_pattern  = params.single_lookup_pattern ?: null
params.accessions    = params.accessions ?: null

workflow {

    sample_id_ch = Channel
        .fromPath(params.sample_list, checkIfExists: true)
        .splitText()
        .map { it.trim() }
        .filter { it && !it.startsWith('#') }

    // Optional materialization step where SRA accessions are prefetch'd and become ordinary local FASTQs.
    if (params.input_type == 'sra') {
        sra_prefetch(Channel.fromPath(params.accessions, checkIfExists: true))
    }

    /*
     * Universal lookup step.
     * This runs regardless of whether FASTQs already existed
     * or were just created by SRA_DOWNLOAD.
     */

    paired_ch = params.paired_reads
        ? Channel
            .fromFilePairs(params.paired_reads, checkIfExists: true)
            .map { sample_id, reads ->
                tuple([id: sample_id, layout: 'paired'], reads)
            }
        : Channel.empty()

    single_ch = params.single_reads
        ? Channel
            .fromPath(params.single_reads, checkIfExists: true)
            .map { fq ->
                tuple([id: fq.baseName, layout: 'single'], [fq])
            }
        : Channel.empty()

    reads_ch = paired_ch.mix(single_ch)

    reads_by_layout = reads_ch.branch {
        paired: meta.layout == 'paired'
        single: meta.layout == 'single'
    }

    paired_clean = PAIRED_KNEADDATA(reads_by_layout.paired)
    single_clean = SINGLE_KNEADDATA(reads_by_layout.single)

    cleaned_ch = paired_clean.out.reads.mix(single_clean.out.reads)

    BIOBAKERY_DOWNSTREAM(cleaned_ch)
}