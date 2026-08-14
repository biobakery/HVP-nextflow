#!/usr/bin/env nextflow
nextflow.enable.dsl=2

import groovy.json.JsonSlurper

include { identify_inputs }              from '../modules/vis/identify_inputs/main.nf'
include { feature_table; trim_taxonomy } from '../modules/stats/feature_table/main.nf'
include { alpha_diversity }              from '../modules/vis/alpha_diversity/main.nf'
include { add_ec_names }                 from '../modules/vis/add_ec_names/main.nf'
include { vis_report }                   from '../modules/vis/report/main.nf'
include { archive_output }               from '../modules/utils/archive/main.nf'

// Visualization workflow — the Nextflow port of `biobakery_workflows vis`.
//
// Like the AnADAMA workflow this takes a folder of bioBakery output and
// discovers the data files in it by their known names, rather than being handed
// individual files. When run at the end of mgx/mtx/16s the folder is the run's
// own output directory, so there is a single code path either way.
//
// The `ready` channel exists only to order this after the profiling stages when
// chained; standalone runs pass a value channel that fires immediately.
workflow VIS {

    take:
    ready

    main:
    def no_file = file("${projectDir}/assets/NO_FILE")
    def input_dir = file(params.vis_input ?: params.outdir)

    if (!input_dir.exists())
        error "ERROR: vis input folder not found: ${input_dir}\nSet --vis_input to a bioBakery output folder."

    metadata_ch = params.input_metadata
        ? Channel.value(file(params.input_metadata))
        : Channel.value(no_file)

    // ── Identify the input files ───────────────────────────────────────────
    folder_ch = ready.first().map { input_dir }

    identify_inputs(folder_ch, metadata_ch, 'vis')

    manifest_ch = identify_inputs.out.manifest
    info_ch     = manifest_ch.map { json -> new JsonSlurper().parse(json.toFile()) }

    // ── Taxonomy feature table (input to the alpha diversity plots) ────────
    // generate_alpha_diversity_plots() calls create_feature_table_inputs first;
    // wmgx profiles go through create_feature_table.py and 16s profiles through
    // trim_taxonomy.py.
    taxonomy_ch = info_ch.map { m -> [ m.study_type, file("${input_dir}/${m.taxonomic_profile}") ] }

    // not published: vis_report publishes the whole `vis` folder wholesale, so
    // this stays an intermediate -- see the note in the feature_table module
    wmgx_taxonomy_ch = taxonomy_ch
        .filter { study_type, profile -> study_type == 'WGX' }
        .map { study_type, profile ->
            [ 'taxonomy', profile, "--sample-tag-column '_taxonomic_profile' --reduce-stratified-species-only", false ] }

    sixteens_taxonomy_ch = taxonomy_ch
        .filter { study_type, profile -> study_type == '16S' }
        .map { study_type, profile -> [ 'taxonomy', profile, false ] }

    features_ch = feature_table(wmgx_taxonomy_ch).features
        .mix(trim_taxonomy(sixteens_taxonomy_ch).features)

    // ── Alpha diversity plots ──────────────────────────────────────────────
    // Only generated when metadata are provided, matching the original.
    if (params.input_metadata) {
        alpha_diversity(features_ch, metadata_ch)
        alpha_ch = alpha_diversity.out.plots.ifEmpty(no_file)
    } else {
        alpha_ch = Channel.value(no_file)
    }

    // ── EC names ───────────────────────────────────────────────────────────
    ecs_input_ch = info_ch
        .filter { m -> m.study_type == 'WGX' && m.ecsabundance }
        .map { m -> file("${input_dir}/${m.ecsabundance}") }

    ecs_ch = add_ec_names(ecs_input_ch).ecs.ifEmpty(no_file)

    // ── Report and archive ─────────────────────────────────────────────────
    vis_report(
        Channel.value(input_dir),
        manifest_ch,
        metadata_ch,
        alpha_ch,
        ecs_ch)

    archive_output(vis_report.out.report_folder, metadata_ch, 'vis')

    emit:
    report  = vis_report.out.report
    archive = archive_output.out.archive
}
