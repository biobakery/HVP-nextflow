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
    // humann_regroup_table and humann_rename_table take different vocabularies:
    // regroup -g wants e.g. uniref90_rxn, rename -n wants one of
    // {kegg-orthology, ec, metacyc-rxn, metacyc-pwy, pfam, eggnog, go,
    // infogo1000, uniref90}. Passing the regroup value straight through made
    // humann_rename_table exit on "invalid choice", so map it here.
    def rename_names = [
        'uniref90_rxn'    : 'metacyc-rxn',
        'uniref90_ko'     : 'kegg-orthology',
        'uniref90_eggnog' : 'eggnog',
        'uniref90_pfam'   : 'pfam',
        'uniref90_go'     : 'go',
        'uniref50_rxn'    : 'metacyc-rxn',
        'uniref50_ko'     : 'kegg-orthology',
        'uniref50_eggnog' : 'eggnog',
        'uniref50_pfam'   : 'pfam',
        'uniref50_go'     : 'go',
    ]
    def names = params.humann_rename_names ?:
                rename_names[params.humann_regroup_grouping] ?:
                params.humann_regroup_grouping
    """
    humann_rename_table \\
        -i ${regrouped} \\
        -n ${names} \\
        -o ${sample}_${params.humann_regroup_grouping}_named_${params.humann_version}.tsv
    """
}
