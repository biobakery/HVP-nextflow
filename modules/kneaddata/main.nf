#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// KneadData QC — single-end reads
process single_end_kneaddata {
    tag "$sample"
    publishDir "${params.outdir}/kneaddata", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_kneaddata.fastq.gz"),      emit: kneads
    tuple val(sample), path("${sample}_kneaddata.log"),            emit: log
    tuple val(sample), path("${sample}_kneaddata*.fastq.gz"),      emit: others, optional: true

    when:
    params.run_qc

    script:

    def extra_args = params.kneaddata_options ?: ""
    def bypass      = params.kneaddata_bypass_trim               ? "--bypass-trim"                  : ""
    def remove_inter = params.kneaddata_remove_intermediate_files ? "--remove-intermediate-output"   : ""
    """
    kneaddata \\
        --unpaired $reads \\
        --reference-db ${params.host_genome} \\
        --output ./ \\
        --threads ${task.cpus} \\
        --output-prefix ${sample}_kneaddata \\
        $bypass \\
        $remove_inter \\
        $extra_args

    for f in *.fastq; do [ -f "\$f" ] && pigz -p ${task.cpus} "\$f"; done
    """
}

// KneadData QC — paired-end reads
process paired_end_kneaddata {
    tag "$sample"
    publishDir "${params.outdir}/kneaddata", mode: 'copy',
        saveAs: { fn -> fn.endsWith('_concatenated.fastq.gz') ? null : fn }

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_concatenated.fastq.gz"),       emit: kneads
    tuple val(sample), path("${sample}_kneaddata.log"),                emit: log
    tuple val(sample), path("${sample}_kneaddata_paired_1.fastq.gz"),  emit: paired1,   optional: true
    tuple val(sample), path("${sample}_kneaddata_paired_2.fastq.gz"),  emit: paired2,   optional: true
    tuple val(sample), path("${sample}_kneaddata_unmatched_1.fastq.gz"), emit: unpaired1, optional: true
    tuple val(sample), path("${sample}_kneaddata_unmatched_2.fastq.gz"), emit: unpaired2, optional: true

    when:
    params.run_qc

    script:

    def extra_args = params.kneaddata_options ?: ""
    def bypass       = params.kneaddata_bypass_trim               ? "--bypass-trim"                 : ""
    def remove_inter = params.kneaddata_remove_intermediate_files ? "--remove-intermediate-output"  : ""
    """
    kneaddata \\
        -i1 ${reads[0]} \\
        -i2 ${reads[1]} \\
        --reference-db ${params.host_genome} \\
        --output ./ \\
        --processes ${task.cpus} \\
        --output-prefix ${sample}_kneaddata \\
        $bypass \\
        $remove_inter \\
        $extra_args

    for f in *.fastq; do [ -f "\$f" ] && pigz -p ${task.cpus} "\$f"; done

    cat ${sample}_kneaddata_paired_1.fastq.gz \\
        ${sample}_kneaddata_paired_2.fastq.gz \\
        ${sample}_kneaddata_unmatched_1.fastq.gz \\
        ${sample}_kneaddata_unmatched_2.fastq.gz \\
        > ${sample}_concatenated.fastq.gz
    """
}
