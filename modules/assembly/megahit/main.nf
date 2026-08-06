#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// MEGAHIT — de novo metagenome assembly
process megahit {
    tag "$sample"
    publishDir "${params.outdir}/assembly/main/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(reads)   // reads may be a single file (single/concatenated) or [r1, r2, u1, u2]

    output:
    tuple val(sample), path("${sample}.final.contigs.fa"), emit: contigs

    script:
    def min_len = params.megahit_min_contig_length ?: 2500
    def extra   = params.megahit_options           ?: ""

    if (params.paired_end) {
        """
        megahit \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            -r ${reads[2]},${reads[3]} \\
            -o megahit_out \\
            --out-prefix ${sample} \\
            -m 0.99 \\
            -t ${task.cpus} \\
            --continue \\
            --min-contig-len ${min_len} \\
            ${extra}

        mv megahit_out/${sample}.contigs.fa ${sample}.final.contigs.fa
        """
    } else {
        """
        megahit \\
            -r ${reads} \\
            -o megahit_out \\
            --out-prefix ${sample} \\
            -m 0.99 \\
            -t ${task.cpus} \\
            --continue \\
            --min-contig-len ${min_len} \\
            ${extra}

        mv megahit_out/${sample}.contigs.fa ${sample}.final.contigs.fa
        """
    }
}
