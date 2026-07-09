#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { humann }         from '../modules/humann/main.nf'
include { humann_regroup } from '../modules/utils/humann_regroup/main.nf'
include { humann_rename }  from '../modules/utils/humann_rename/main.nf'
include { humann_merge }   from '../modules/utils/humann_merge/main.nf'

// Functional profiling subworkflow: HUMAnN + optional regroup / rename / merge
workflow FUNCTIONAL_PROFILING {

    take:
    reads    // Channel: [ [id: sample, paired_end: bool], reads ]
    profiles // Channel: [ sample, metaphlan_profile.tsv ]

    main:
    reads_flat = reads.map { meta, r -> tuple(meta.id, r) }

    // Each sample's reads must be paired with its MetaPhlAn profile
    humann_input = reads_flat.join(profiles)

    humann_reads   = humann_input.map { sample, reads, profile -> tuple(sample, reads) }
    humann_profile = humann_input.map { sample, reads, profile -> tuple(sample, profile) }

    humann_out = humann(humann_reads, humann_profile)

    // Regroup gene families to alternate annotation scheme
    regroup_out = humann_regroup(humann_out.genefamilies)

    // Add human-readable names to regrouped features
    rename_out = humann_rename(regroup_out.regrouped)

    // Merge all per-sample tables into a single matrix
    all_tables = humann_out.genefamilies.map { s, f -> f }
        .mix(humann_out.pathabundance.map { s, f -> f })
        .collect()
    merge_out  = humann_merge(all_tables)

    emit:
    genefamilies  = humann_out.genefamilies
    pathabundance = humann_out.pathabundance
    pathcoverage  = humann_out.pathcoverage
    regrouped     = regroup_out.regrouped
    renamed       = rename_out.renamed
}
