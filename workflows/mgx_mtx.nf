#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MGX } from './mgx.nf'

// Paired MGX + MTX workflow — stub
// Planned: run QC → MetaPhlAn → HUMAnN on both MGX and MTX reads,
// then integrate the two profiles.
workflow MGX_MTX {
    main:
    log.warn "mgx_mtx is not yet fully implemented. Running MGX half only."
    MGX()
}
