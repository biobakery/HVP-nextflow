#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Normalise one sample's counts by its Bracken-classified read total
process phanta_normalize_counts {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy',
               enabled: params.save_intermediate

    input:
    tuple val(meta), path(counts_table)
    path total_reads
    path failed

    output:
    tuple val(meta), path("${meta.id}.krak.report_bracken_species.filtered.to_merge.norm_brack"), emit: normed

    script:
    def sample = meta.id
    // Downstream of the second cohort barrier: this needs the whole cohort's
    // total_reads.tsv to find its own row. Upstream does every sample in one
    // invocation by scanning a shared directory; one sample per task is what
    // makes the step isolatable.
    """
    phanta_prep_to_merge_normed.py \\
        --counts-table ${counts_table} \\
        --total-reads ${total_reads} \\
        --failed-file ${failed} \\
        --sample ${sample} \\
        --out ${sample}.krak.report_bracken_species.filtered.to_merge.norm_brack
    """
}
