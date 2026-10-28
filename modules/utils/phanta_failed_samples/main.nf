#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Cohort barrier — the list of samples Bracken could not process
process phanta_failed_samples {
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy'

    input:
    path statuses

    output:
    path "samples_that_failed_bracken.txt", emit: failed

    script:
    // The first of two cohort barriers: every sample's filtering has to have
    // finished before this list is complete, and prepare_counts downstream
    // needs the whole list. Sorted by name, matching upstream -- which gets the
    // order from a shell glob rather than from the samplesheet.
    //
    // Upstream builds this by globbing for a canary file; here it reads the
    // explicit status records, so an empty list means "every sample passed"
    // rather than "the glob matched nothing".
    """
    cat ${statuses} | awk -F'\\t' '\$2 == "ineligible" { print \$1 }' | sort \\
        > samples_that_failed_bracken.txt
    """
}
