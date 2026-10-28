#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Kraken2 — classify reads against the Phanta database
process kraken2 {
    tag "${meta.id}"
    // Publishing contract. The raw per-read .krak is NEVER published: it is the
    // one artifact whose content depends on Kraken2's thread scheduling, so it
    // differs between identical runs, and it is the largest file here. The
    // minimizer report is a stable intermediate, so it follows --save_intermediate.
    // Nothing is deleted either way -- both stay in the task directory for -resume.
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy',
               saveAs: { fn -> (fn.endsWith('.krak') || !params.save_intermediate) ? null : fn }

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("${meta.id}.krak.report"), emit: report
    tuple val(meta), path("${meta.id}.krak"),        emit: assignments

    script:
    def sample = meta.id
    // --paired is per sample here, where upstream makes it a cohort-wide
    // setting. The samplesheet carries layout per row, so a cohort can be
    // mixed; the workflow still rejects that, but the module does not have to.
    def paired   = meta.paired_end ? '--paired' : ''
    def gzipped  = params.phanta_gzipped ? '--gzip-compressed' : ''
    def extra    = params.phanta_options ?: ''
    """
    kraken2 \\
        --db ${db} \\
        --threads ${task.cpus} \\
        --output ${sample}.krak \\
        --report ${sample}.krak.report \\
        --report-minimizer-data \\
        ${paired} \\
        ${gzipped} \\
        ${reads} \\
        --confidence ${params.phanta_confidence} \\
        ${extra}

    # The eight-column minimizer report is what the filtering step needs; a
    # Kraken2 older than 2.1.2 silently emits six and the run fails much later
    # with a confusing pandas error instead of here.
    columns=\$(awk -F'\\t' 'NR==1{print NF; exit}' ${sample}.krak.report)
    if [ "\$columns" != "8" ]; then
        echo "ERROR: ${sample}.krak.report has \$columns columns, expected 8." >&2
        echo "       --report-minimizer-data needs Kraken2 >= 2.1.2." >&2
        exit 1
    fi
    """
}
