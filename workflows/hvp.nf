#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ═════════════════════════════════════════════════════════════════════════════
// HVP orchestrator
//
// Six stages, in order. Each stage reads one kind of data and produces the next
// (conf/components.json names them), and runs whichever of its components were
// selected. The orchestrator only wires; every tool runs inside a component.
//
//   1 preprocess            raw reads      → clean reads     chain
//   2 read_profiling        clean reads    → read profiles   parallel
//   3 assembly              clean reads    → contigs         parallel
//   4 virus_identification  contigs        → viral contigs   parallel
//   5 virus_curation        viral contigs  → viral genomes   chain
//   6 characterisation      viral genomes  → annotations     parallel
//
// chain:    components run one after another; with none selected the stage
//           passes its input straight through.
// parallel: every selected component reads the same input; outputs are mixed.
//
// To add a component: register it in conf/components.json, write
// subworkflows/local/components/<name>/main.nf to its stage's contract, and add
// one `if (on('<name>'))` block in that stage below.
// ═════════════════════════════════════════════════════════════════════════════

include { merge_pairs          } from '../modules/local/merge_pairs/main.nf'
include { version_log          } from '../modules/local/version_log/main.nf'

// 1 · preprocess
include { KNEADDATA            } from '../subworkflows/local/components/kneaddata/main.nf'
include { HOST_SCRUBBER        } from '../subworkflows/local/components/host_scrubber/main.nf'
// 2 · read_profiling
include { BAQLAVA              } from '../subworkflows/local/components/baqlava/main.nf'
include { PHANTA               } from '../subworkflows/local/components/phanta/main.nf'
include { SYLPH                } from '../subworkflows/local/components/sylph/main.nf'
// 3 · assembly
include { ASSEMBLY             } from '../subworkflows/local/components/assembly/main.nf'
include { METASPADES           } from '../subworkflows/local/components/metaspades/main.nf'
// 4 · virus_identification
include { GENOMAD              } from '../subworkflows/local/components/genomad/main.nf'
include { CENOTE_TAKER3        } from '../subworkflows/local/components/cenote_taker3/main.nf'
// 5 · virus_curation
include { CHECKV_CONTAMINATION } from '../subworkflows/local/components/checkv_contamination/main.nf'
include { CIRCULAR_CONTIGS     } from '../subworkflows/local/components/circular_contigs/main.nf'
include { VIRUS_RENAME         } from '../subworkflows/local/components/virus_rename/main.nf'
// 6 · characterisation
include { CHECKV_COMPLETENESS  } from '../subworkflows/local/components/checkv_completeness/main.nf'
include { VIRUS_TAXONOMY       } from '../subworkflows/local/components/virus_taxonomy/main.nf'

workflow HVP {

    take:
    ch_reads    // Channel: [ [id: sample, paired_end: bool], reads ]
    components  // List<String>: resolved by resolveComponents(), in pipeline order

    main:
    def on = { String name -> name in components }

    // ── 1 · preprocess: raw reads → clean reads (chain) ─────────────────────
    // Two channels leave this stage: one read file per sample for the
    // profilers, and the separate mates of paired samples for the assemblers.
    if (on('kneaddata')) {
        KNEADDATA(ch_reads)
        ch_clean = KNEADDATA.out.reads
        ch_parts = KNEADDATA.out.parts
    } else {
        // No read QC: pairs are concatenated for the profilers (as KneadData
        // would have done), and the raw mates go to the assemblers.
        ch_parts = ch_reads.filter { meta, r -> meta.paired_end }.map { meta, r -> tuple(meta.id, r) }
        ch_clean = merge_pairs(ch_parts).merged
            .map { sample, f -> [ [id: sample, paired_end: true], f ] }
            .mix( ch_reads.filter { meta, r -> !meta.paired_end } )
    }
    if (on('host_scrubber')) {
        HOST_SCRUBBER(ch_clean, ch_parts)
        ch_clean = HOST_SCRUBBER.out.reads
        ch_parts = HOST_SCRUBBER.out.parts
    }

    // ── 2 · read_profiling: clean reads → read profiles (parallel) ──────────
    ch_profiles = Channel.empty()
    if (on('baqlava')) {
        BAQLAVA(ch_clean)
        ch_profiles = ch_profiles.mix(BAQLAVA.out.profiles)
    }
    if (on('phanta')) {
        PHANTA(ch_clean)
        ch_profiles = ch_profiles.mix(PHANTA.out.profiles)
    }
    if (on('sylph')) {
        SYLPH(ch_clean)
        ch_profiles = ch_profiles.mix(SYLPH.out.profiles)
    }

    // ── 3 · assembly: clean reads → contigs (parallel) ──────────────────────
    ch_contigs = Channel.empty()
    if (on('assembly')) {
        ASSEMBLY(ch_clean, ch_parts)
        ch_contigs = ch_contigs.mix(ASSEMBLY.out.contigs)
    }
    if (on('metaspades')) {
        METASPADES(ch_clean, ch_parts)
        ch_contigs = ch_contigs.mix(METASPADES.out.contigs)
    }

    // ── 4 · virus_identification: contigs → viral contigs (parallel) ────────
    ch_viral = Channel.empty()
    if (on('genomad')) {
        GENOMAD(ch_contigs)
        ch_viral = ch_viral.mix(GENOMAD.out.viral_contigs)
    }
    if (on('cenote_taker3')) {
        CENOTE_TAKER3(ch_contigs)
        ch_viral = ch_viral.mix(CENOTE_TAKER3.out.viral_contigs)
    }

    // ── 5 · virus_curation: viral contigs → viral genomes (chain) ───────────
    ch_genomes = ch_viral
    if (on('checkv_contamination')) {
        CHECKV_CONTAMINATION(ch_genomes)
        ch_genomes = CHECKV_CONTAMINATION.out.genomes
    }
    if (on('circular_contigs')) {
        CIRCULAR_CONTIGS(ch_genomes)
        ch_genomes = CIRCULAR_CONTIGS.out.genomes
    }
    if (on('virus_rename')) {
        VIRUS_RENAME(ch_genomes)
        ch_genomes = VIRUS_RENAME.out.genomes
    }

    // ── 6 · characterisation: viral genomes → annotations (parallel) ────────
    ch_annotations = Channel.empty()
    if (on('checkv_completeness')) {
        CHECKV_COMPLETENESS(ch_genomes)
        ch_annotations = ch_annotations.mix(CHECKV_COMPLETENESS.out.annotations)
    }
    if (on('virus_taxonomy')) {
        VIRUS_TAXONOMY(ch_genomes)
        ch_annotations = ch_annotations.mix(VIRUS_TAXONOMY.out.annotations)
    }

    // ── Run record: tool versions, databases, components ────────────────────
    if (params.log_versions) {
        version_log(components.join(','))
    }

    emit:
    clean_reads   = ch_clean
    read_profiles = ch_profiles
    contigs       = ch_contigs
    viral_contigs = ch_viral
    viral_genomes = ch_genomes
    annotations   = ch_annotations
}
