#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PhyloPhlAn metagenomic — phylogenetic placement of MAGs to nearest SGBs/GGBs/FGBs
process phylophlan_metagenomic {
    tag "$sample"
    publishDir "${params.outdir}/phylophlan/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bins_dir)

    output:
    tuple val(sample), path("phylophlan_out.tsv"), emit: placement

    script:
    def db_folder = params.phylophlan_path
    def db        = file(params.phylophlan_path).list()?.find { it.count('.') == 1 } ?: ""
    def extra     = params.phylophlan_metagenomic_options ?: ""
    """
    if ls ${bins_dir}/*.bin.[0-9]*.fa 2>/dev/null | grep -q .; then
        # Stage only the numbered bins, as the anadama2 workflow's copy_bins step
        # does. MetaBAT also emits empty lowDepth/tooShort/unbinned placeholders,
        # and mash cannot sketch a zero-byte file: phylophlan reports
        # "error while sketching ... expected str, bytes or os.PathLike object,
        # not NoneType" and exits.
        mkdir -p numbered_bins
        cp ${bins_dir}/*.bin.[0-9]*.fa numbered_bins/
        phylophlan_assign_sgbs_legacy \\
            -i numbered_bins \\
            -n 1 \\
            --add_ggb --add_fgb \\
            -d ${db} \\
            -o phylophlan_out \\
            --nproc ${task.cpus} \\
            --verbose \\
            -e fa \\
            --database_folder ${db_folder} \\
            ${extra}
        mv phylophlan_out.tsv phylophlan_out.tsv
    else
        echo -e "line1\nline2\nline3\n#mag\tsgb\tggb\tfgb\tref" > phylophlan_out.tsv
    fi
    """
}

// Merge per-sample PhyloPhlAn placements and add taxonomic labels
process phylophlan_merge {
    publishDir "${params.outdir}/phylophlan", mode: 'copy'

    input:
    path placements  // collected phylophlan_out.tsv files

    output:
    path "phylophlan_out.tsv",    emit: merged
    path "phylophlan_relab.tsv",  emit: relab

    script:
    """
    # Merge: keep first file header, append data rows from the rest
    head -1 \$(ls *.tsv | head -1) > phylophlan_out.tsv
    for f in *.tsv; do
        tail -n +5 "\$f" >> phylophlan_out.tsv
    done

    python ${projectDir}/bin/scripts/phylophlan_add_tax_assignment.py \\
        --table phylophlan_out.tsv \\
        --output phylophlan_relab.tsv
    """
}
