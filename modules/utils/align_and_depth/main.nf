#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Align reads to assembled contigs and compute contig depth with jgi_summarize_bam_contig_depths
process align_and_depth {
    tag "$sample"
    publishDir "${params.outdir}/assembly/contig_depths", mode: 'copy'

    input:
    tuple val(sample), path(contigs)
    tuple val(sample2), path(reads)  // same order as megahit: [r1, r2, u1, u2] or single

    output:
    tuple val(sample), path("${sample}.contig_depths.txt"), emit: depths
    tuple val(sample), path("${sample}.sorted.bam"),        emit: bam,    optional: true

    script:
    def idx = "${sample}_btidx"
    if (params.paired_end) {
        """
        if [ ! -s ${contigs} ]; then
            touch ${sample}.sorted.bam ${sample}.contig_depths.txt
        else
            bowtie2-build ${contigs} ${idx}
            bowtie2 -x ${idx} \\
                -1 ${reads[0]} -2 ${reads[1]} \\
                -U ${reads[2]},${reads[3]} \\
                -S ${sample}.sam \\
                -p ${task.cpus} \\
                --very-sensitive-local --no-unal
            samtools view -bS -F 4 ${sample}.sam > ${sample}.unsorted.bam
            samtools sort ${sample}.unsorted.bam -o ${sample}.sorted.bam --threads ${task.cpus}
        fi
        jgi_summarize_bam_contig_depths \\
            --outputDepth ${sample}.contig_depths.txt \\
            ${sample}.sorted.bam
        """
    } else {
        """
        if [ ! -s ${contigs} ]; then
            touch ${sample}.sorted.bam ${sample}.contig_depths.txt
        else
            bowtie2-build ${contigs} ${idx}
            bowtie2 -x ${idx} \\
                -U ${reads} \\
                -S ${sample}.sam \\
                -p ${task.cpus} \\
                --very-sensitive-local --no-unal
            samtools view -bS -F 4 ${sample}.sam > ${sample}.unsorted.bam
            samtools sort ${sample}.unsorted.bam -o ${sample}.sorted.bam --threads ${task.cpus}
        fi
        jgi_summarize_bam_contig_depths \\
            --outputDepth ${sample}.contig_depths.txt \\
            ${sample}.sorted.bam
        """
    }
}
