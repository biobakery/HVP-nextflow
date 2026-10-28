#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { kraken2 }                 from '../modules/kraken2/main.nf'
include { bracken }                 from '../modules/bracken/main.nf'
include { phanta_filter }           from '../modules/taxonomy/phanta_filter/main.nf'
include { phanta_scale }            from '../modules/taxonomy/phanta_scale/main.nf'
include { phanta_db_manifest }      from '../modules/utils/phanta_db_manifest/main.nf'
include { phanta_failed_samples }   from '../modules/utils/phanta_failed_samples/main.nf'
include { phanta_prepare_counts }   from '../modules/utils/phanta_prepare_counts/main.nf'
include { phanta_total_reads }      from '../modules/utils/phanta_prepare_counts/main.nf'
include { phanta_normalize_counts } from '../modules/utils/phanta_normalize_counts/main.nf'
include { phanta_merge_tables }     from '../modules/utils/phanta_merge_tables/main.nf'
include { version_log }             from '../modules/utils/version_log/main.nf'

// Phanta — read-level taxonomic profiling including the gut virome
workflow PHANTA {

    main:

    // ── Input ────────────────────────────────────────────────────────────
    // Phanta is the one workflow that takes a samplesheet instead of
    // --readsdir. Cohort aggregation needs an explicit, ordered sample list:
    // the merged tables' column order follows the samplesheet, and filename
    // globbing cannot express an order.
    def rows = phanta_samplesheet_rows(params.phanta_samplesheet)

    // The order is read once, eagerly, and carried as a value. Everything
    // downstream that has to be deterministic is ordered against this list
    // rather than against channel arrival order, which is nondeterministic.
    def ordered_samples = rows.collect { it.sample }

    ch_reads = Channel
        .fromList(rows)
        .map { row ->
            def files = row.fastq_2 ? [ file(row.fastq_1), file(row.fastq_2) ]
                                    : [ file(row.fastq_1) ]
            [ [id: row.sample, paired_end: (row.fastq_2 as boolean)], files ]
        }

    ch_db = Channel.value(file(params.phanta_db))

    // ── Provenance, before anything expensive ────────────────────────────
    // Validates every required database file and fails on a checksum mismatch.
    // Kraken2 takes ch_manifest as an input purely to order it after this, so
    // that a run against the wrong database stops here rather than an hour in.
    phanta_db_manifest(ch_db)
    ch_manifest = phanta_db_manifest.out.manifest

    // Everything downstream takes the database through THIS channel rather than
    // ch_db. It cannot emit until the manifest task has succeeded, which is what
    // actually enforces "a checksum mismatch fails the run before any Kraken2
    // task" -- a comment saying so is not a dependency.
    ch_db_checked = ch_manifest.map { file(params.phanta_db) }

    if (params.log_versions) {
        version_log()
    }

    // ── Fan-out: classify, filter, quantify, scale ───────────────────────
    kraken2(ch_reads, ch_db_checked)
    phanta_filter(kraken2.out.report, ch_db_checked)

    // Bracken needs the filtered report staged under its exact upstream name,
    // because it derives its second output's filename from its input.
    ch_bracken_in = phanta_filter.out.filtered
        .join(phanta_filter.out.status)
        .join(phanta_filter.out.substitutes.map { meta, b, s, _sc -> [ meta, b, s ] })
        .map { meta, filtered, status, sub_b, sub_s -> [ meta, filtered, status, sub_b, sub_s ] }

    bracken(ch_bracken_in, ch_db_checked)

    ch_scale_in = bracken.out.table
        .join(phanta_filter.out.decisions)
        .join(phanta_filter.out.status)
        .join(phanta_filter.out.substitutes.map { meta, _b, _s, sc -> [ meta, sc ] })

    phanta_scale(ch_scale_in, ch_db_checked)

    // ── Cohort barrier 1: which samples Bracken could not process ────────
    // collect() waits for every sample, then combine() puts the one resulting
    // file back onto each per-sample item. This is the explicit form of the
    // feedback upstream gets for free by having every rule read a shared
    // directory -- which is also what made upstream impossible to isolate.
    phanta_failed_samples(
        phanta_filter.out.status.map { _meta, status -> status }.collect()
    )
    ch_failed = phanta_failed_samples.out.failed

    ch_counts_in = bracken.out.species_report
        .join(kraken2.out.report)
        .join(bracken.out.table)

    // No .first() on ch_failed: a process whose inputs are all value channels
    // emits a value channel, and phanta_failed_samples takes a collect(). So
    // this is already reused for every sample. (Adding .first() here was
    // harmless but wrong, and Nextflow warns about it.)
    phanta_prepare_counts(ch_counts_in, ch_failed, ch_db_checked)

    // ── Cohort barrier 2: the cohort read totals ─────────────────────────
    phanta_total_reads(phanta_prepare_counts.out.total_reads.collect())
    ch_total_reads = phanta_total_reads.out.total_reads

    phanta_normalize_counts(phanta_prepare_counts.out.to_merge, ch_total_reads, ch_failed)

    // ── Aggregate ────────────────────────────────────────────────────────
    // Sorted by sample name here only to make the staged file list stable; the
    // merge step itself works from the explicit ordered_samples list, so the
    // published column order is the samplesheet's regardless of what order
    // these channels happened to emit in.
    phanta_merge_tables(
        phanta_prepare_counts.out.to_merge.map { _meta, f -> f }.collect(sort: true),
        phanta_normalize_counts.out.normed.map { _meta, f -> f }.collect(sort: true),
        phanta_scale.out.scaled.map { _meta, f -> f }.collect(sort: true),
        ordered_samples,
        ch_db_checked
    )

    emit:
    counts              = phanta_merge_tables.out.counts
    relative_read       = phanta_merge_tables.out.relative_read
    relative_taxonomic  = phanta_merge_tables.out.relative_taxonomic
    total_reads         = ch_total_reads
    failed              = ch_failed
    manifest            = ch_manifest
}

// Parse and validate the samplesheet.
//
// Eager Groovy rather than a channel operator, for the same reason
// read_input.nf is a function: the order has to be known before any channel is
// built, and the validation errors are far more useful raised once, up front,
// than as a failing task per bad row.
//
// Deliberately does NOT check that the FASTQ files exist. `-preview` builds the
// graph without running a task, and bin/make_diagrams.py previews every
// workflow in CI against stub inputs; a file-existence check here would make
// that fail. Existence is verified inside the Kraken2 task, where it belongs.
def phanta_samplesheet_rows(samplesheet) {

    if (!samplesheet) {
        error "ERROR: --phanta_samplesheet is required for --workflow phanta.\n" +
              "  A CSV with the header 'sample,fastq_1,fastq_2'; leave fastq_2 empty\n" +
              "  for single-end samples. Convert Phanta's own tab-delimited sample\n" +
              "  file with bin/phanta_convert_samplesheet.py."
    }

    def handle = file(samplesheet)
    if (!handle.exists()) {
        error "ERROR: samplesheet not found: ${samplesheet}"
    }

    def lines = handle.readLines().findAll { it.trim() }
    if (lines.size() < 2) {
        error "ERROR: samplesheet ${samplesheet} has no data rows."
    }

    def header = lines[0].split(',', -1).collect { it.trim() }
    if (header[0] != 'sample' || header[1] != 'fastq_1') {
        error "ERROR: samplesheet ${samplesheet} must start with the header\n" +
              "  sample,fastq_1,fastq_2\n" +
              "  Found: ${lines[0]}"
    }

    def rows = []
    def seen = [:]
    lines.drop(1).eachWithIndex { line, index ->
        def cells = line.split(',', -1).collect { it.trim() }
        def lineno = index + 2

        if (cells.size() < 2) {
            error "ERROR: ${samplesheet} line ${lineno}: expected at least 'sample,fastq_1'.\n" +
                  "  Found: ${line}"
        }

        def sample  = cells[0]
        def fastq_1 = cells[1]
        def fastq_2 = cells.size() > 2 ? cells[2] : ''

        if (!sample)  error "ERROR: ${samplesheet} line ${lineno}: empty sample name."
        if (!fastq_1) error "ERROR: ${samplesheet} line ${lineno}: sample '${sample}' has no fastq_1."

        // A duplicate ID would silently overwrite a column in the merged
        // tables, so it is fatal rather than a warning.
        if (seen.containsKey(sample)) {
            error "ERROR: ${samplesheet} line ${lineno}: duplicate sample '${sample}',\n" +
                  "  already seen on line ${seen[sample]}. Sample names become column\n" +
                  "  names in the merged tables and must be unique."
        }
        seen[sample] = lineno

        rows << [sample: sample, fastq_1: fastq_1, fastq_2: fastq_2]
    }

    // Kraken2 takes --paired for the whole run, and the read-length doubling in
    // the scaling step depends on the layout, so a mixed cohort would be
    // quietly wrong rather than merely unusual.
    def paired = rows.findAll { it.fastq_2 }
    if (paired && paired.size() != rows.size()) {
        def singles = rows.findAll { !it.fastq_2 }.collect { it.sample }
        error "ERROR: ${samplesheet} mixes paired-end and single-end samples.\n" +
              "  ${paired.size()} paired, ${singles.size()} single-end " +
              "(${singles.take(5).join(', ')}${singles.size() > 5 ? ', ...' : ''}).\n" +
              "  Phanta applies one library layout to the whole cohort. Split the run."
    }

    if (!params.phanta_db) {
        error "ERROR: --phanta_db is required for --workflow phanta.\n" +
              "  Set it on the command line, or use a -profile that provides it."
    }

    log.info "[phanta] ${rows.size()} samples, " +
             "${paired ? 'paired-end' : 'single-end'}, in samplesheet order: " +
             "${rows.collect { it.sample }.take(5).join(', ')}" +
             "${rows.size() > 5 ? ', ...' : ''}"

    return rows
}
