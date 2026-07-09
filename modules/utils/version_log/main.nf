#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Capture tool and database versions for reproducibility
process version_log {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    output:
    path "versions.txt",  emit: versions
    path "db_paths.txt",  emit: db_paths

    script:
    def baqlava_db_str    = params.baqlava_db    ?: 'default (bundled)'
    def strainphlan_db_str = params.strainphlan_db ?: 'auto-detect from metaphlan_db'
    """
    {
        echo "=== Tool Versions ==="
        echo "Date: \$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
        echo ""
        kneaddata --version  2>&1 | head -1 || echo "kneaddata: not in PATH"
        metaphlan --version  2>&1 | head -1 || echo "metaphlan: not in PATH"
        humann --version     2>&1 | head -1 || echo "humann: not in PATH"
        which baqlava       >/dev/null 2>&1 && baqlava --version 2>&1 | head -1 || echo "baqlava: not installed"
        which strainphlan.py >/dev/null 2>&1 && strainphlan.py --version 2>&1 | head -1 || echo "strainphlan: not installed"
    } > versions.txt

    {
        echo "=== Database Paths ==="
        echo "human_genome:    ${params.human_genome}"
        echo "metaphlan_db:    ${params.metaphlan_db}"
        echo "metaphlan_index: ${params.metaphlan_index}"
        echo "humann_db:       ${params.humann_db}"
        echo "baqlava_db:      ${baqlava_db_str}"
        echo "strainphlan_db:  ${strainphlan_db_str}"
    } > db_paths.txt
    """
}
