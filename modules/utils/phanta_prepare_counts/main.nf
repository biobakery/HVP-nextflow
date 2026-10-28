#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Build one sample's mergeable counts table and its read-accounting row
process phanta_prepare_counts {
    tag "${meta.id}"
    publishDir path: { "${params.outdir}/phanta/classification" }, mode: 'copy',
               enabled: params.save_intermediate

    input:
    tuple val(meta), path(species_report), path(krak_report), path(bracken_table)
    path failed
    path db

    output:
    tuple val(meta), path("${meta.id}.krak.report_bracken_species.filtered.to_merge"), emit: to_merge
    path "${meta.id}_total_reads.txt",                                                 emit: total_reads

    script:
    def sample = meta.id
    // Two upstream rules folded together: prepare_to_merge_counts and
    // tot_reads_per_sample read the same sample's files at the same point, with
    // no barrier between them.
    //
    // calc_total_reads exits 64 when a sample has no reads at all. Letting that
    // propagate is deliberate -- it means the sample produced nothing, which is
    // a real error rather than the biological-empty case that the Bracken
    // eligibility path already handles.
    """
    phanta_prep_to_merge_counts.py \\
        --brack-report ${species_report} \\
        --db ${db} \\
        --failed-file ${failed} \\
        --sample ${sample} \\
        --out ${sample}.krak.report_bracken_species.filtered.to_merge

    phanta_calc_total_reads.sh \\
        --sample ${sample} \\
        --krak-report ${krak_report} \\
        --bracken-table ${bracken_table} \\
        --out ${sample}_total_reads.txt
    """
}

// Cohort barrier — assemble total_reads.tsv, rows sorted by sample name
process phanta_total_reads {
    publishDir path: { "${params.outdir}/phanta/final_merged_outputs" }, mode: 'copy'

    input:
    path rows

    output:
    path "total_reads.tsv", emit: total_reads

    script:
    // Note the ordering: rows here are sorted by NAME, while the three merged
    // tables put their sample columns in SAMPLESHEET order. Those two orders
    // are genuinely different and both match upstream, which gets this one from
    // a shell glob and the other from its sample file. Sorting is done at emit
    // time because Nextflow channel arrival order is nondeterministic.
    """
    printf 'Samp_Name\\tTot_Samp_Reads\\tUnassigned_Step_One\\tAssigned_Step_Three\\tUnassigned_Step_Three\\n' \\
        > total_reads.tsv
    for f in \$(ls ${rows} | sort); do cat "\$f"; done >> total_reads.tsv
    """
}
