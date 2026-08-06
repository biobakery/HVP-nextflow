#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Mash sketch → paste → dist pipeline for pairwise genome distance estimation
process mash_sketch {
    publishDir "${params.outdir}/sgbs/mash", mode: 'copy'

    input:
    path mag_filepaths   // text file listing MAG FASTA paths (from mash_list_inputs.py)

    output:
    path "sketches.msh", emit: sketch

    script:
    def extra = params.mash_sketch_options ?: ""
    """
    if [ ! -s ${mag_filepaths} ]; then
        touch sketches.msh
    else
        mash sketch -p ${task.cpus} -l ${mag_filepaths} -o sketches ${extra}
    fi
    """
}

process mash_paste {
    publishDir "${params.outdir}/sgbs/mash", mode: 'copy'

    input:
    path sketch

    output:
    path "references.msh", emit: references

    script:
    """
    if [ ! -s ${sketch} ]; then
        touch references.msh
    else
        mash paste references ${sketch}
    fi
    """
}

process mash_dist {
    publishDir "${params.outdir}/sgbs/mash", mode: 'copy'

    input:
    path references
    path sketch

    output:
    path "mash_dist_out.tsv", emit: distances

    script:
    """
    if [ ! -s ${references} ]; then
        touch mash_dist_out.tsv
    else
        mash dist -p ${task.cpus} -t ${references} ${sketch} > mash_dist_out.tsv
    fi
    """
}

// List qualifying MAG FASTA paths to feed into Mash
process mash_list_inputs {
    publishDir "${params.outdir}/sgbs/mash", mode: 'copy'

    input:
    path checkm_qa_n50
    path phylophlan_relab
    path bins_root   // parent directory containing all per-sample bins/

    output:
    path "mags_filepaths.txt", emit: filepaths

    script:
    """
    python ${projectDir}/bin/scripts/mash_list_inputs.py \\
        --checkm ${checkm_qa_n50} \\
        --phylophlan ${phylophlan_relab} \\
        --bins ${bins_root} \\
        --mash . \\
        --threads ${task.cpus}
    """
}

// Cluster MAGs into SGBs using Mash distances + CheckM/PhyloPhlAn metadata
process sgb_cluster {
    publishDir "${params.outdir}/sgbs/sgbs", mode: 'copy'

    input:
    path mash_distances
    path phylophlan_relab
    path checkm_qa_n50
    path bins_root   // per-sample bins root; mash_list_inputs writes qc_bins/ inside it

    output:
    path "SGB_info.tsv",     emit: sgb_info
    path "SGB_list.txt",     emit: sgb_list, optional: true

    script:
    """
    if [ ! -s ${mash_distances} ]; then
        mkdir -p fastANI
        touch fastANI/SGB_list.txt
        echo -e "cluster\tgenome\tcluster_members\tn_genomes\tcompleteness\tcontamination\tstrain_heterogeneity\tn50\tquality\tkeep\tcluster_name\tsgb" > SGB_info.tsv
    else
        Rscript ${projectDir}/bin/Rscripts/mash_clusters.R \\
            --mash ${mash_distances} \\
            --checkm ${checkm_qa_n50} \\
            --phylo ${phylophlan_relab} \\
            --out_dir . \\
            --threads ${task.cpus} \\
            --mag_dir ${bins_root}/qc_bins

        # mash_clusters.R writes into <out_dir>/sgbs/ and <out_dir>/fastANI/, but
        # the declared outputs are at the top level (and that is where the branch
        # above writes them). Without this the task exits 0 and nextflow then
        # reports the outputs as missing.
        cp sgbs/SGB_info.tsv SGB_info.tsv
        if [ -f fastANI/SGB_list.txt ]; then cp fastANI/SGB_list.txt SGB_list.txt; fi
    fi
    """
}

// Merge abundance data with taxonomy and SGB assignments into the final profile
process merge_tax_abundance {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path abundance_tables   // collected *.abundance.tsv, staged into this directory
    path phylophlan_relab
    path checkm_qa_n50
    path sgb_info

    output:
    path "final_profile.tsv", emit: final_profile

    script:
    def abundance_type = params.sgb_abundance_type ?: "by_sample"
    """
    Rscript ${projectDir}/bin/Rscripts/merge_tax_and_abundance.R \\
        -i . \\
        --tax ${phylophlan_relab} \\
        --qa ${checkm_qa_n50} \\
        --sgbs ${sgb_info} \\
        -o final_profile.tsv
    """
}
