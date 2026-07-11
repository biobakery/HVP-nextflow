#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Regroup HUMAnN gene families to a different annotation scheme
process humann_regroup {
    tag "$sample"
    publishDir "${params.outdir}/humann/${params.humann_version}/regroup", mode: 'copy'

    input:
    tuple val(sample), path(genefamilies)

    output:
    tuple val(sample), path("${sample}_${params.humann_regroup_grouping}_${params.humann_version}.tsv"), emit: regrouped

    when:
    params.run_humann_regroup

    script:
    """
    humann_regroup_table \\
        -i ${genefamilies} \\
        -g ${params.humann_regroup_grouping} \\
        -o ${sample}_${params.humann_regroup_grouping}_${params.humann_version}.tsv
    """
}
