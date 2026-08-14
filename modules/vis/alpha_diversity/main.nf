#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Alpha diversity plots from the taxonomy feature table.
// Ports utilities.generate_alpha_diversity_plots() (utilities.py 538-553).
//
// The feature table this consumes is built by the shared feature_table process,
// exactly as the original calls create_feature_table_inputs() first. Only runs
// when metadata are provided.
// No publishDir: the report driver copies these plots into the report folder,
// which vis_report publishes as ${params.outdir}/vis. Publishing them here too
// would target that same folder, and vis_report's later publish of the whole
// directory would clobber them.
process alpha_diversity {
    input:
    tuple val(feature_type), path(features)
    path metadata

    output:
    path "alpha_diversity_plots", emit: plots

    script:
    """
    Rscript \$(python -c "import biobakery_bootstrap; print(biobakery_bootstrap.get_rscript('alpha_diversity'))") \\
        ${features} \\
        ${metadata} \\
        alpha_diversity_plots \\
        --max_missing ${params.max_missing}

    # the report only checks the folder contents, so make sure it exists even
    # when every variable was filtered out for missing values
    mkdir -p alpha_diversity_plots
    """
}
