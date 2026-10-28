#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Bracken — re-estimate species abundance from the filtered report
process bracken {
    tag "${meta.id}"
    // The Kraken-style species report is a published product; the .bracken table
    // it comes with is a standard Bracken intermediate and follows
    // --save_intermediate. For an ineligible sample both are genuinely empty
    // files, and publishing the empty one is the point -- that is the
    // observable contract the fixtures record.
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy',
               saveAs: { fn -> fn.endsWith('_bracken_species.filtered') ? fn
                             : (params.save_intermediate ? fn : null) }

    input:
    tuple val(meta), path(filtered), path(status), path(sub_bracken), path(sub_species)
    path db

    output:
    tuple val(meta), path("${meta.id}.krak.report.filtered.bracken"),          emit: table
    tuple val(meta), path("${meta.id}.krak.report_bracken_species.filtered"),  emit: species_report

    script:
    def sample = meta.id
    // Bracken declares two outputs but is only told about one. The second
    // appears because Bracken derives its name from the INPUT filename:
    //
    //   input    {samp}.krak.report.filtered
    //   -o       {samp}.krak.report.filtered.bracken
    //   implied  {samp}.krak.report_bracken_species.filtered
    //
    // So the filtered report must be staged under exactly that name. Any
    // wrapper that renames its input breaks the second output silently.
    """
    if grep -q \$'\\tineligible\\t' ${status}; then
        # Bracken would error out on this sample. Upstream's substitutes are
        # genuinely empty files, and the fixtures contain them, so copy rather
        # than skip: downstream steps expect the files to exist.
        cp ${sub_bracken} ${sample}.krak.report.filtered.bracken
        cp ${sub_species} ${sample}.krak.report_bracken_species.filtered
    else
        bracken \\
            -d ${db} \\
            -i ${filtered} \\
            -o ${sample}.krak.report.filtered.bracken \\
            -r ${params.phanta_read_length} \\
            -l 'S' \\
            -t ${params.phanta_filter_thresh}
    fi

    if [ ! -f ${sample}.krak.report_bracken_species.filtered ]; then
        echo "ERROR: Bracken did not produce ${sample}.krak.report_bracken_species.filtered." >&2
        echo "       That file is named after the INPUT, so the filtered report must be" >&2
        echo "       staged as exactly ${sample}.krak.report.filtered." >&2
        exit 1
    fi
    """
}
