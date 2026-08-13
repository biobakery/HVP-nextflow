process SRA_PREFETCH {

    tag "$sample_id"

    publishDir params.rawfastq_dir, mode: 'copy', overwrite: false

    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}*.fastq*"), emit: fastqs

    script:
    def fasterq_cmd = params.fasterq_dump_execline ?: 'fasterq-dump'
    def split_strat  = params.sra_split_strat ?: '--split-files'
    def n_threads    = task.cpus

    """
    mkdir -p ${params.rawfastq_dir}

    found_files=\$(find ${params.rawfastq_dir} -maxdepth 1 -type f -name "*${sample_id}*" | wc -l)

    if [ "\$found_files" -gt 0 ]; then
        echo "Found \$found_files FASTQ file(s) mapping to sample ID ${sample_id}; skipping fasterq-dump."

        cp ${params.rawfastq_dir}/*${sample_id}* ./
    else
        echo "No existing FASTQ files found for sample ID ${sample_id}; running fasterq-dump."

        ${fasterq_cmd} \\
            ${sample_id} \\
            -O . \\
            ${split_strat} \\
            -e ${n_threads}

        gzip -f *.fastq
    fi

    for f in *.fastq.gz; do
        case "\$f" in
            *_1.fastq.gz) mv "\$f" ${sample_id}_R1.fastq.gz ;;
            *_2.fastq.gz) mv "\$f" ${sample_id}_R2.fastq.gz ;;
            *.fastq.gz)   mv "\$f" ${sample_id}.fastq.gz ;;
        esac
    done
    """
}