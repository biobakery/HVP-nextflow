#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Build the input read channel, detecting library layout from the filenames
// rather than making the user assert it.
//
// Behaviour, following the anadama2 workflow:
//   * look for the pair identifier in the read filenames
//   * if it is never found, warn once and run everything single-end
//   * if it is found, warn again listing any sample that turned up without a
//     mate, since that is almost always a file naming mistake rather than a
//     genuinely single-end sample
//   * --single_end forces single-end even when pairs are present
//
// Emits: [ [id: sample, paired_end: bool], reads ]
workflow READ_INPUT {

    main:
    if (!params.readsdir) {
        error "ERROR: --readsdir is required"
    }

    def readsdir = params.readsdir - ~/\/$/

    // Explicit filepattern wins; otherwise derive it from the layout so a
    // single-end run does not have to restate the pattern.
    def paired_glob = params.filepattern ?: params.filepattern_paired
    def single_glob = params.filepattern ?: params.filepattern_single

    // Decide the layout before building any channel, so the warnings are emitted
    // once for the run rather than once per file.
    def pair_files = params.single_end ? [] : file("${readsdir}/${paired_glob}")
    def use_paired = pair_files.size() > 0

    if (params.single_end) {
        log.info "Running single-end: --single_end was set."
    }
    else if (!use_paired) {
        log.warn "No files matched the paired-end pattern '${paired_glob}' in ${readsdir}. " +
                 "Running in single-end mode. Set --filepattern if your reads use a different " +
                 "naming convention, or --single_end true to silence this."
    }

    if (use_paired) {
        reads = Channel
            .fromFilePairs("${readsdir}/${paired_glob}", checkIfExists: true)
            .map { sample, files -> [ [id: sample, paired_end: true], files ] }

        // Any read file that did not pair up is reported: fromFilePairs silently
        // drops singletons, which makes a naming mixup very easy to miss.
        def all_reads = file("${readsdir}/${single_glob}").collect { it.name }
        def paired_names = pair_files.collect { it.name }
        def orphans = all_reads - paired_names
        if (orphans) {
            log.warn "${orphans.size()} read file(s) in ${readsdir} did not match a mate and " +
                     "will be skipped: ${orphans.sort().take(10).join(', ')}" +
                     (orphans.size() > 10 ? ", ..." : "") +
                     ". Check the pair identifier ('${params.pair_identifier}') and file naming."
        }
    }
    else {
        reads = Channel
            .fromPath("${readsdir}/${single_glob}", checkIfExists: true)
            .map { f ->
                def sample = f.baseName.replaceFirst(/(\.fastq|\.fq)(\.gz)?$/, '')
                [ [id: sample, paired_end: false], f ]
            }
    }

    emit:
    reads
}
