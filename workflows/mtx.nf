#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { read_input }           from '../subworkflows/read_input.nf'
include { QUALITY_CONTROL }      from '../subworkflows/quality_control.nf'
include { TAXONOMIC_PROFILING }  from '../subworkflows/taxonomic_profiling.nf'
include { FUNCTIONAL_PROFILING } from '../subworkflows/functional_profiling.nf'
include { STRAIN_PROFILING }     from '../subworkflows/strain_profiling.nf'
include { version_log }          from '../modules/utils/version_log/main.nf'
include { mtx_kneaddata_dbs }    from '../subworkflows/mtx_common.nf'

// Whole Metatranscriptome (MTX) workflow.
//
// biobakery_workflows 3.2 has no standalone metatranscriptome script: the
// metatranscriptome half only exists inside wmgx_wmtx.py. This workflow is that
// half run on its own, for the case where there is no paired metagenome — the
// same taxonomic and functional profiling as mgx, but with the metatranscriptome
// KneadData database set (host genome + host mRNA + rRNA) rather than the host
// genome alone.
//
// Without a metagenome to borrow taxonomy from, the MetaPhlAn profile is
// computed from the RNA reads themselves and handed to HUMAnN, which is what
// wmgx_wmtx.py does when --input-mapping is not given.
workflow MTX {

    main:
    read_ch = read_input(params.readsdir, 'mtx')

    // ── QC (KneadData, metatranscriptome database set) ────────────────────
    if (params.run_qc) {
        QUALITY_CONTROL(read_ch, mtx_kneaddata_dbs(), '')
        cleaned = QUALITY_CONTROL.out.reads
    } else {
        cleaned = read_ch
    }

    // ── Taxonomic profiling (MetaPhlAn) ───────────────────────────────────
    if (params.run_taxonomic_profiling) {
        TAXONOMIC_PROFILING(cleaned, '')
    }

    if (!params.run_taxonomic_profiling) {
        def dependents = []
        if (params.run_functional_profiling) dependents << 'run_functional_profiling (HUMAnN needs the MetaPhlAn profile)'
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

    // ── Strain profiling (StrainPhlAn) ────────────────────────────────────
    if (params.run_strain_profiling) {
        STRAIN_PROFILING(TAXONOMIC_PROFILING.out.sam_bzip)
    }

    if (params.log_versions) {
        version_log()
    }
}
