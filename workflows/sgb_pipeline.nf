#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUALITY_CONTROL }        from '../subworkflows/quality_control.nf'
include { megahit }                from '../modules/sgb_pipeline/megahit/main.nf'
include { align_and_depth }        from '../modules/utils/align_and_depth/main.nf'
include { metabat2 }               from '../modules/binning/metabat2/main.nf'
include { checkm2 }                from '../modules/qc/checkm2/main.nf'
include { checkm2_merge }          from '../modules/qc/checkm2/main.nf'
include { mag_n50 }                from '../modules/qc/checkm2/main.nf'
include { checkm2_wrangling }      from '../modules/qc/checkm2/main.nf'
include { phylophlan_metagenomic } from '../modules/phylogenomics/phylophlan_metagenomic/main.nf'
include { phylophlan_merge }       from '../modules/phylogenomics/phylophlan_metagenomic/main.nf'
include { mash_list_inputs }       from '../modules/utils/mash/main.nf'
include { mash_sketch }            from '../modules/utils/mash/main.nf'
include { mash_paste }             from '../modules/utils/mash/main.nf'
include { mash_dist }              from '../modules/utils/mash/main.nf'
include { sgb_cluster }            from '../modules/utils/mash/main.nf'
include { merge_tax_abundance }    from '../modules/utils/mash/main.nf'
include { version_log }            from '../modules/utils/version_log/main.nf'

// Full MAG assembly → binning → SGB clustering pipeline.
// Ported from anadama2 biobakery-workflows feature/sgb_pipeline branch.
workflow SGB_PIPELINE {

    main:

    // ── Build input channel ───────────────────────────────────────────────
    if (!params.readsdir) {
        error "ERROR: --readsdir is required for sgb_pipeline workflow"
    }
    if (params.paired_end) {
        reads = Channel
            .fromFilePairs("${params.readsdir}/${params.filepattern}", checkIfExists: true)
            .map { sample, files -> [ [id: sample, paired_end: true], files ] }
    } else {
        reads = Channel
            .fromPath("${params.readsdir}/${params.filepattern}", checkIfExists: true)
            .map { f ->
                def sample = f.baseName.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, '')
                [ [id: sample, paired_end: false], f ]
            }
    }

    // ── Step 1: Host decontamination (KneadData) ──────────────────────────
    QUALITY_CONTROL(reads)
    cleaned = QUALITY_CONTROL.out.reads

    cleaned_flat = cleaned.map { meta, r -> tuple(meta.id, r) }

    // ── Step 2: De novo assembly (MEGAHIT) ────────────────────────────────
    contigs_out = megahit(cleaned_flat)

    // ── Step 3: Align reads to contigs + contig depth ─────────────────────
    contigs_with_reads = contigs_out.contigs.join(cleaned_flat)
    depth_out = align_and_depth(
        contigs_with_reads.map { s, c, r -> tuple(s, c) },
        contigs_with_reads.map { s, c, r -> tuple(s, r) }
    )

    // ── Step 4: Bin contigs into MAGs (MetaBAT2) ──────────────────────────
    bins_input = contigs_out.contigs.join(depth_out.depths)
    bins_out = metabat2(
        bins_input.map { s, c, d -> tuple(s, c) },
        bins_input.map { s, c, d -> tuple(s, d) }
    )

    // ── Step 5: MAG quality assessment (CheckM2) ──────────────────────────
    checkm_out    = checkm2(bins_out.bins)

    all_n50_input = bins_out.bins.map { s, b -> b }.collect()
    n50_out       = mag_n50(all_n50_input)

    all_reports   = checkm_out.quality.map { s, r -> r }.collect()
    merged_checkm = checkm2_merge(all_reports)
    qa_n50        = checkm2_wrangling(merged_checkm.merged, n50_out.n50)

    // ── Step 6: Phylogenetic placement (PhyloPhlAn metagenomic) ───────────
    phylo_out      = phylophlan_metagenomic(bins_out.bins)
    all_placements = phylo_out.placement.map { s, p -> p }.collect()
    phylo_merged   = phylophlan_merge(all_placements)

    // ── Step 7: SGB clustering (Mash + R) ────────────────────────────────
    bins_root  = bins_out.bins.map { s, b -> b.parent }.first()
    mag_list   = mash_list_inputs(qa_n50.qa_n50, phylo_merged.relab, bins_root)
    sketch     = mash_sketch(mag_list.filepaths)
    references = mash_paste(sketch.sketch)
    distances  = mash_dist(references.references, sketch.sketch)

    sgbs = sgb_cluster(distances.distances, phylo_merged.relab, qa_n50.qa_n50, bins_root)

    // ── Step 8: Merge abundance + taxonomy → final profile ────────────────
    abundance_dir = depth_out.bam.map { s, b -> b.parent }.first()
    merge_tax_abundance(
        abundance_dir,
        phylo_merged.relab,
        qa_n50.qa_n50,
        sgbs.sgb_info
    )

    // ── Version / DB log ──────────────────────────────────────────────────
    if (params.log_versions) {
        version_log()
    }
}
