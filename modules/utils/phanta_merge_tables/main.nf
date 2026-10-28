#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Merge the per-sample tables into the three cohort tables
process phanta_merge_tables {
    publishDir path: { "${params.outdir}/phanta/final_merged_outputs" }, mode: 'copy'

    input:
    path to_merge
    path normed
    path scaled
    val  ordered_samples
    path db

    output:
    path "counts.txt",                        emit: counts
    path "relative_read_abundance.txt",       emit: relative_read
    path "relative_taxonomic_abundance.txt",  emit: relative_taxonomic

    script:
    // Column order is the whole point of passing an explicit sample list.
    // Nextflow channel arrival order is nondeterministic, so the file lists
    // staged here are in no particular order; the workflow sorts them back into
    // samplesheet order before this runs, and the sample IDs below are in that
    // same order. The merge scripts pair the two positionally.
    def order     = ordered_samples.join(' ')
    def to_merge_files = ordered_samples.collect { "${it}.krak.report_bracken_species.filtered.to_merge" }.join(' ')
    def normed_files   = ordered_samples.collect { "${it}.krak.report_bracken_species.filtered.to_merge.norm_brack" }.join(' ')
    def scaled_files   = ordered_samples.collect { "${it}.krak.report.filtered.bracken.scaled" }.join(' ')
    """
    phanta_merge.py \\
        --input-file ${to_merge_files} \\
        --sample-id ${order} \\
        --output-file counts.txt

    phanta_merge.py \\
        --input-file ${normed_files} \\
        --sample-id ${order} \\
        --output-file relative_read_abundance.txt

    phanta_merging_scaled_bracken.py \\
        --db ${db} \\
        --input-file ${scaled_files} \\
        --sample-id ${order} \\
        --output-file relative_taxonomic_abundance.txt
    """
}
