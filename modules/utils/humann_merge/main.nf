#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Merge per-sample HUMAnN outputs into a single matrix
process humann_merge {
    publishDir "${params.outdir}/humann/${params.humann_version}/merged", mode: 'copy'

    input:
    path tables

    output:
    path "merged_genefamilies_${params.humann_version}.tsv",   emit: genefamilies,  optional: true
    path "merged_pathabundance_${params.humann_version}.tsv",  emit: pathabundance, optional: true
    path "merged_pathcoverage_${params.humann_version}.tsv",   emit: pathcoverage,  optional: true

    when:
    params.run_humann_merge

    script:
    """
    if ls *_genefamilies_*.tsv 1>/dev/null 2>&1; then
        humann_join_tables -i . --file-name genefamilies -o merged_genefamilies_${params.humann_version}.tsv
    fi
    if ls *_pathabundance_*.tsv 1>/dev/null 2>&1; then
        humann_join_tables -i . --file-name pathabundance -o merged_pathabundance_${params.humann_version}.tsv
    fi
    if ls *_pathcoverage_*.tsv 1>/dev/null 2>&1; then
        humann_join_tables -i . --file-name pathcoverage -o merged_pathcoverage_${params.humann_version}.tsv
    fi
    """
}
