#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUALITY_CONTROL }      from '../subworkflows/quality_control.nf'
include { TAXONOMIC_PROFILING }  from '../subworkflows/taxonomic_profiling.nf'
include { FUNCTIONAL_PROFILING } from '../subworkflows/functional_profiling.nf'
include { VIRAL_PROFILING }      from '../subworkflows/viral_profiling.nf'
include { STRAIN_PROFILING }     from '../subworkflows/strain_profiling.nf'
include { version_log }          from '../modules/utils/version_log/main.nf'

// Whole Metagenome Shotgun (MGX) workflow
workflow MGX {

    main:
    // ── Build input channel ────────────────────────────────────────────────
    if (params.paired_end) {
        read_ch = Channel
            .fromFilePairs("${params.readsdir}/${params.filepattern}", checkIfExists: true)
            .map { sample, reads -> [ [id: sample, paired_end: true], reads ] }
    } else {
        read_ch = Channel
            .fromPath("${params.readsdir}/${params.filepattern}", checkIfExists: true)
            .map { f ->
                def sample = f.baseName.replaceFirst(/(\.fastq|\.fq)$/, '')
                [ [id: sample, paired_end: false], f ]
            }
    }

    // ── QC (KneadData) ────────────────────────────────────────────────────
    if (params.run_qc) {
        QUALITY_CONTROL(read_ch)
        cleaned = QUALITY_CONTROL.out.reads
    } else {
        cleaned = read_ch
    }

    // ── Taxonomic profiling (MetaPhlAn) ───────────────────────────────────
    if (params.run_taxonomic_profiling) {
        TAXONOMIC_PROFILING(cleaned)
    }

    // ── Functional profiling (HUMAnN) ─────────────────────────────────────
    if (params.run_functional_profiling && params.run_taxonomic_profiling) {
        FUNCTIONAL_PROFILING(cleaned, TAXONOMIC_PROFILING.out.profile)
    }

    // ── Viral profiling (BAQLaVa) ──────────────────────────────────────────
    if (params.run_viral_profiling && params.run_taxonomic_profiling) {
        VIRAL_PROFILING(cleaned, TAXONOMIC_PROFILING.out.profile)
    } else if (params.run_viral_profiling) {
        error "ERROR: run_viral_profiling requires run_taxonomic_profiling = true (MetaPhlAn profile needed by BAQLaVa)"
    }

    // ── Strain profiling (StrainPhlAn) ────────────────────────────────────
    if (params.run_strain_profiling && params.run_taxonomic_profiling) {
        STRAIN_PROFILING(TAXONOMIC_PROFILING.out.sam_bzip)
    }

    // ── Version logging ────────────────────────────────────────────────────
    if (params.log_versions) {
        version_log()
    }
}
