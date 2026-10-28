#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Scale a Bracken table by genome length
process phanta_scale {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy'

    input:
    tuple val(meta), path(bracken), path(decisions), path(status), path(sub_scaled)
    path db

    output:
    tuple val(meta), path("${meta.id}.krak.report.filtered.bracken.scaled"), emit: scaled

    script:
    def sample = meta.id
    // Upstream passes the Python-literal spelling of the layout flag straight
    // through to the script, which compares it against the string 'True'.
    def paired = meta.paired_end ? 'True' : 'False'
    """
    if grep -q \$'\\tineligible\\t' ${status}; then
        cp ${sub_scaled} ${sample}.krak.report.filtered.bracken.scaled
    else
        phanta_scale_bracken.py \\
            --db ${db} \\
            --bracken-report ${bracken} \\
            --filtering-decisions ${decisions} \\
            --read-length ${params.phanta_read_length} \\
            --paired ${paired} \\
            --out ${sample}.krak.report.filtered.bracken.scaled
    fi
    """
}
