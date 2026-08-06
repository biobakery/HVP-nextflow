#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MetaPhlAn — taxonomic profiling
process metaphlan {
    tag "$sample"
    publishDir "${params.outdir}/metaphlan/${params.metaphlan_index}", mode: 'copy'

    input:
    tuple val(sample), path(kneads)

    output:
    tuple val(sample), path("${sample}_profile_${params.metaphlan_index}.tsv"),   emit: profile
    tuple val(sample), path("${sample}_bowtie2_${params.metaphlan_index}.tsv"),   emit: bowtie2
    tuple val(sample), path("${sample}_${params.metaphlan_index}.sam"),            emit: sam

    when:
    params.run_taxonomic_profiling

    script:

    def extra_args = params.metaphlan_options ?: ""
    // Resolve CLI flag differences between MetaPhlAn versions
    def db_arg  = (params.metaphlan_version == 'metaphlan_v4')
                    ? "--db_dir ${params.metaphlan_db}"
                    : "--bowtie2db ${params.metaphlan_db}"
    def out_arg = (params.metaphlan_version == 'metaphlan_v4')
                    ? "--mapout ${sample}_bowtie2_${params.metaphlan_index}.tsv"
                    : "--bowtie2out ${sample}_bowtie2_${params.metaphlan_index}.tsv"
    """
    metaphlan $kneads \\
        --input_type fastq \\
        -t ${params.metaphlan_analysis_type} \\
        --index ${params.metaphlan_index} \\
        $db_arg \\
        $out_arg \\
        --samout ${sample}_${params.metaphlan_index}.sam \\
        --read_min_len ${params.metaphlan_read_min_len} \\
        --nproc ${task.cpus} \\
        -o ${sample}_profile_${params.metaphlan_index}.tsv \\
        $extra_args
    """
}

// Compress MetaPhlAn SAM file (saves significant disk space)
process metaphlan_bzip {
    tag "$sample"
    publishDir "${params.outdir}/metaphlan/bzip", mode: 'copy'

    input:
    tuple val(sample), path(sam)

    output:
    tuple val(sample), path("${sample}_${params.metaphlan_index}.sam.bz2"), emit: sam_bzip

    script:
    """
    bzip2 -c -- "${sam}" > "${sample}_${params.metaphlan_index}.sam.bz2"
    """
}

// Merge per-sample MetaPhlAn profiles into a single table
process metaphlan_merge {
    publishDir "${params.outdir}/metaphlan", mode: 'copy'

    input:
    path profiles

    output:
    path "merged_metaphlan_profiles.tsv", emit: merged

    script:
    """
    merge_metaphlan_tables.py $profiles > merged_metaphlan_profiles.tsv
    """
}
