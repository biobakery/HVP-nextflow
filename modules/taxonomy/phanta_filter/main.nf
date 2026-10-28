#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Filter false-positive species, and decide whether Bracken can run
process phanta_filter {
    tag "${meta.id}"
    // The filtering decisions are a published product: they are the record of
    // which species were dropped and why, which is the whole point of this step.
    // The filtered report, the .species working file and the zero-byte
    // substitutes are intermediates.
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy',
               saveAs: { fn -> fn.endsWith('.filtering_decisions.txt') ? fn
                             : (params.save_intermediate ? fn : null) }

    input:
    tuple val(meta), path(report)
    path db

    output:
    tuple val(meta), path("${meta.id}.krak.report.filtered"),                 emit: filtered
    tuple val(meta), path("${meta.id}.krak.report.filtering_decisions.txt"),  emit: decisions
    tuple val(meta), path("${meta.id}.bracken_status.txt"),                   emit: status
    tuple val(meta), path("${meta.id}.krak.report.filtered.bracken.temp"),
                     path("${meta.id}.krak.report_bracken_species.filtered.temp"),
                     path("${meta.id}.krak.report.filtered.bracken.scaled.temp"), emit: substitutes
    path "${meta.id}.krak.report.species",                                    emit: species

    script:
    def sample = meta.id
    // Two upstream rules folded into one process: filter_kraken and
    // process_filtered_kraken read the same input at the same point in the
    // dataflow with no barrier between them, so splitting them costs a
    // container start and a stage round-trip per sample and buys nothing.
    //
    // The threshold order below is NOT grouped by domain upstream, which is why
    // these are passed by name.
    """
    phanta_filter_kraken_reports.py \\
        --krak-report ${report} \\
        --db ${db} \\
        --cov-thresh-bacterial ${params.phanta_cov_thresh_bacterial} \\
        --cov-thresh-viral ${params.phanta_cov_thresh_viral} \\
        --minimizer-thresh-bacterial ${params.phanta_minimizer_thresh_bacterial} \\
        --minimizer-thresh-viral ${params.phanta_minimizer_thresh_viral} \\
        --cov-thresh-arc ${params.phanta_cov_thresh_arc} \\
        --cov-thresh-euk ${params.phanta_cov_thresh_euk} \\
        --minimizer-thresh-arc ${params.phanta_minimizer_thresh_arc} \\
        --minimizer-thresh-euk ${params.phanta_minimizer_thresh_euk} \\
        --out-species ${sample}.krak.report.species \\
        --out-filtered ${sample}.krak.report.filtered \\
        --out-decisions ${sample}.krak.report.filtering_decisions.txt

    phanta_process_filtered_kraken.py \\
        --filtered-report ${sample}.krak.report.filtered \\
        --filter-thresh ${params.phanta_filter_thresh} \\
        --sample ${sample} \\
        --out-status ${sample}.bracken_status.txt \\
        --out-bracken ${sample}.krak.report.filtered.bracken.temp \\
        --out-species-report ${sample}.krak.report_bracken_species.filtered.temp \\
        --out-scaled ${sample}.krak.report.filtered.bracken.scaled.temp

    # The substitutes are only written when the sample is ineligible, but the
    # output block declares them unconditionally, so make sure they exist.
    # They are meant to be empty: downstream copies them into place as the
    # sample's real products.
    touch ${sample}.krak.report.filtered.bracken.temp \\
          ${sample}.krak.report_bracken_species.filtered.temp \\
          ${sample}.krak.report.filtered.bracken.scaled.temp
    """
}
