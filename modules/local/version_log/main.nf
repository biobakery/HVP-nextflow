#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Capture tool versions, database paths, components and parameters for reproducibility
process version_log {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    val components   // comma-separated: the components this run resolved to

    output:
    path "pipeline_log.txt", emit: pipeline_log

    script:
    def baqlava_db_str = params.baqlava_db     ?: 'default (bundled)'
    def phylophlan_str = params.phylophlan_path ?: 'not set'
    def selected       = components.tokenize(',')
    // Only ask for the versions of tools this run actually used.
    def probes = []
    if ('kneaddata' in selected) probes << 'kneaddata --version'
    if ('baqlava'   in selected) probes += ['metaphlan --version', 'humann --version', 'baqlava --version']
    if ('assembly'  in selected) probes += ['megahit --version', 'metabat2 --help', 'checkm2 --version', 'phylophlan --version', 'mash --version']
    def version_lines = probes.collect { cmd ->
        def tool = cmd.tokenize(' ')[0]
        "if which ${tool} >/dev/null 2>&1; then echo \"${tool}: \$(${cmd} 2>&1 | head -1)\"; else echo \"${tool}: not in PATH\"; fi"
    }.join('\n        ')
    """
    {
        echo "================================================================"
        echo "  biobakery/hvp — Pipeline Run Log"
        echo "================================================================"
        echo ""
        echo "Date:       \$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
        echo "Version:    ${workflow.manifest.version}"
        echo "Components: ${selected.join(' -> ')}"
        echo "Plugin:     ${params.plugin ?: '-'}"
        echo ""

        echo "----------------------------------------------------------------"
        echo "  Tool Versions"
        echo "----------------------------------------------------------------"
        ${version_lines}
        echo ""

        echo "----------------------------------------------------------------"
        echo "  Database Paths"
        echo "----------------------------------------------------------------"
        echo "host_genome:      ${params.host_genome ?: '-'}"
        echo "metaphlan_db:     ${params.metaphlan_db ?: '-'}"
        echo "metaphlan_index:  ${params.metaphlan_index}"
        echo "humann_db:        ${params.humann_db ?: '-'}"
        echo "baqlava_db:       ${baqlava_db_str}"
        echo "phylophlan_path:  ${phylophlan_str}"
        echo ""

        echo "----------------------------------------------------------------"
        echo "  Parameters"
        echo "----------------------------------------------------------------"
        echo "readsdir:                    ${params.readsdir}"
        echo "outdir:                      ${params.outdir}"
        echo "paired_end / single_end:     ${params.paired_end} / ${params.single_end}"
        echo "filepattern:                 ${params.filepattern ?: '-'}"
        echo ""
        echo "kneaddata_bypass_trim:       ${params.kneaddata_bypass_trim}"
        echo "metaphlan_version:           ${params.metaphlan_version}"
        echo "metaphlan_analysis_type:     ${params.metaphlan_analysis_type}"
        echo "baqlava_humann:              ${params.baqlava_humann}"
        echo "humann_version:              ${params.humann_version}"
        echo "baqlava_bypass_depletion:    ${params.baqlava_bypass_depletion}"
        echo "megahit_min_contig_length:   ${params.megahit_min_contig_length}"
        echo "sgb_completeness / contam.:  ${params.sgb_completeness} / ${params.sgb_contamination}"
        echo ""
        echo "max_memory / cpus / time:    ${params.max_memory} / ${params.max_cpus} / ${params.max_time}"
        echo "================================================================"

    } > pipeline_log.txt
    """
}
