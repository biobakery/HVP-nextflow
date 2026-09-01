#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Join per-sample HUMAnN tables of one feature type into a single matrix.
//
// One invocation per feature type ("label"), rather than one invocation that
// greps for three basename patterns: the tables are staged into their own
// directory so humann_join_tables can simply take the whole folder, which is
// what it is built to do. Labels used by the workflow are
// genefamilies | ecs | pathabundance and their _relab variants.
process humann_join {
    tag "$label"
    publishDir path: { "${params.outdir}/${subdir}humann/${params.humann_version}/merged" }, mode: 'copy'

    input:
    tuple val(label), path(tables, stageAs: 'tables/*')
    val subdir

    output:
    tuple val(label), path("merged_${label}_${params.humann_version}.tsv"), emit: merged

    when:
    params.run_humann_merge

    script:
    """
    humann_join_tables \\
        --input tables \\
        --output merged_${label}_${params.humann_version}.tsv
    """
}

// Count how many features each sample has above zero, per feature type.
// --ignore-un-features drops UNMAPPED/UNINTEGRATED, --ignore-stratification
// counts community-level features only.
process humann_count_features {
    tag "$label"
    publishDir path: { "${params.outdir}/${subdir}humann/${params.humann_version}/counts" }, mode: 'copy'

    input:
    tuple val(label), path(merged_relab)
    val subdir

    output:
    tuple val(label), path("humann_${label}_relab_counts.tsv"), emit: counts

    when:
    params.run_humann_merge

    script:
    """
    count_features.py \\
        --input ${merged_relab} \\
        --output humann_${label}_relab_counts.tsv \\
        --reduce-sample-name \\
        --ignore-un-features \\
        --ignore-stratification
    """
}

// Merge the three per-feature-type count tables into one
process humann_feature_counts_merge {
    publishDir path: { "${params.outdir}/${subdir}humann/${params.humann_version}/counts" }, mode: 'copy'

    input:
    path counts, stageAs: 'counts/*'
    val subdir

    output:
    path "humann_feature_counts.tsv", emit: feature_counts

    when:
    params.run_humann_merge

    script:
    """
    humann_join_tables \\
        --input counts \\
        --output humann_feature_counts.tsv \\
        --file_name _relab_counts.tsv
    """
}
