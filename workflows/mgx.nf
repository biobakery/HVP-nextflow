#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { read_input }           from '../subworkflows/read_input.nf'
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
    // Layout is detected from the filenames; see subworkflows/read_input.nf
    read_ch = read_input(params.readsdir, 'mgx')

    // ── QC (KneadData) ────────────────────────────────────────────────────
    // A metagenome run decontaminates against the host genome only.
    if (params.run_qc) {
        QUALITY_CONTROL(read_ch, [params.host_genome], '')
        cleaned = QUALITY_CONTROL.out.reads
    } else {
        cleaned = read_ch
    }

    // ── Taxonomic profiling (MetaPhlAn) ───────────────────────────────────
    if (params.run_taxonomic_profiling) {
        TAXONOMIC_PROFILING(cleaned, '')
    }

    // Functional, viral and strain profiling all consume MetaPhlAn output, so
    // fail fast and consistently rather than silently skipping the stage.
    if (!params.run_taxonomic_profiling) {
        def dependents = []
        if (params.run_functional_profiling) dependents << 'run_functional_profiling (HUMAnN needs the MetaPhlAn profile)'
        if (params.run_viral_profiling)      dependents << 'run_viral_profiling (BAQLaVa needs the MetaPhlAn profile)'
        if (params.run_strain_profiling)     dependents << 'run_strain_profiling (StrainPhlAn needs the MetaPhlAn SAM)'
        if (dependents) {
            error "ERROR: run_taxonomic_profiling = false, but these require it:\n  - " +
                  dependents.join("\n  - ") +
                  "\nEnable run_taxonomic_profiling, or turn the above off."
        }
    }

    // ── Functional profiling (HUMAnN) ─────────────────────────────────────
    if (params.run_functional_profiling) {
        FUNCTIONAL_PROFILING(cleaned, TAXONOMIC_PROFILING.out.profile, '')
    }

    // ── Viral profiling (BAQLaVa) ──────────────────────────────────────────
    if (params.run_viral_profiling) {
        VIRAL_PROFILING(cleaned, TAXONOMIC_PROFILING.out.profile)
    }

    // ── Strain profiling (StrainPhlAn) ────────────────────────────────────
    if (params.run_strain_profiling) {
        STRAIN_PROFILING(TAXONOMIC_PROFILING.out.sam_bzip)
    }

    // ── Version logging ────────────────────────────────────────────────────
    if (params.log_versions) {
        version_log()
    }
}
