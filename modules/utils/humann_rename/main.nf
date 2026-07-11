#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Rename HUMAnN output features to human-readable names
process humann_rename {
    tag "$sample"
    publishDir "${params.outdir}/humann/${params.humann_version}/rename", mode: 'copy'

    input:
    tuple val(sample), path(regrouped)

    output:
    tuple val(sample), path("${sample}_${params.humann_regroup_grouping}_named_${params.humann_version}.tsv"), emit: renamed

    when:
    params.run_humann_rename

    script:
    """
    humann_rename_table \\
        -i ${regrouped} \\
        -n ${params.humann_regroup_grouping} \\
        -o ${sample}_${params.humann_regroup_grouping}_named_${params.humann_version}.tsv
    """
}
