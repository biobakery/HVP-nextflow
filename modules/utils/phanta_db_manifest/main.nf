#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Record which database this run used, and refuse to run against the wrong one
process phanta_db_manifest {
    publishDir path: { "${params.outdir}/pipeline_info" }, mode: 'copy'

    input:
    path db

    output:
    path "phanta_db_manifest.txt", emit: manifest

    script:
    // Upstream leaves no trace of which database produced a result: the
    // database arrives as a bare path from a config file, and nothing about it
    // is written anywhere. Given a finished Phanta output directory there is no
    // way to tell which build made it, so results cannot be audited or
    // reproduced without out-of-band notes.
    //
    // The manifest is a SIDECAR on purpose. Putting provenance in a header, a
    // comment or a column of the four tables would break byte-identity with the
    // baseline immediately, because provenance varies per run. It lives in
    // pipeline_info/, which the output comparison excludes.
    //
    // inspect.out is the fingerprint that matters: a text dump of every taxon
    // with its clade k-mer counts, so it is content-derived -- two different
    // builds cannot share its checksum -- and small enough to hash every run.
    // hash.k2d is 20+ GB, so a full checksum of it costs a complete read and is
    // gated behind phanta_db_checksum_full.
    // NOT `?: ''`. Nextflow coerces a CLI value that looks numeric into a
    // number, and Groovy's Elvis operator treats 0 as falsy -- so an all-digit
    // checksum became '' and the comparison below was silently skipped, letting
    // the run proceed against an unverified database. Test for null explicitly
    // and stringify.
    def expected = params.phanta_db_expected_inspect_sha256 == null
                     ? '' : params.phanta_db_expected_inspect_sha256.toString()
    """
    set -euo pipefail

    sha_of()  { sha256sum "\$1" | cut -d' ' -f1; }
    size_of() { wc -c < "\$1" | tr -d ' '; }

    # Every file the pipeline actually reads. Checked here, before anything
    # expensive runs, rather than failing inside a Kraken2 task an hour in.
    missing=0
    for f in hash.k2d taxo.k2d opts.k2d seqid2taxid.map inspect.out \\
             taxonomy/nodes.dmp taxonomy/names.dmp library/species_genome_size.txt; do
        if [ ! -s "${db}/\$f" ]; then
            echo "ERROR: required database file missing or empty: \$f" >&2
            missing=1
        fi
    done

    # Bracken's distribution is built per read length; the wrong one is simply
    # absent, so name the parameter that has to change.
    kmer_distrib="${db}/database${params.phanta_read_length}mers.kmer_distrib"
    if [ ! -s "\$kmer_distrib" ]; then
        echo "ERROR: no Bracken distribution for read length ${params.phanta_read_length}." >&2
        echo "       Expected: \$(basename \$kmer_distrib)" >&2
        echo "       Available: \$(ls ${db}/database*mers.kmer_distrib 2>/dev/null | xargs -n1 basename | tr '\\n' ' ')" >&2
        echo "       Set --phanta_read_length to one the database ships." >&2
        missing=1
    fi
    [ "\$missing" = "0" ] || exit 1

    inspect_sha=\$(sha_of "${db}/inspect.out")

    {
        echo -e "phanta_db\\t${params.phanta_db}"
        echo -e "phanta_db_realpath\\t\$(realpath ${db})"
        echo -e "generated_utc\\t\$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
        echo -e "hash.k2d.size\\t\$(size_of ${db}/hash.k2d)"
    """ + (params.phanta_db_checksum_full ? """
        echo -e "hash.k2d.sha256\\t\$(sha_of ${db}/hash.k2d)"
    """ : """
        echo -e "hash.k2d.sha256\\tnot computed (set --phanta_db_checksum_full)"
    """) + """
        for f in taxo.k2d opts.k2d seqid2taxid.map inspect.out library/species_genome_size.txt; do
            echo -e "\$f.size\\t\$(size_of ${db}/\$f)"
            echo -e "\$f.sha256\\t\$(sha_of ${db}/\$f)"
        done
        echo -e "kmer_distrib.file\\t\$(basename \$kmer_distrib)"
        echo -e "kmer_distrib.sha256\\t\$(sha_of \$kmer_distrib)"
        echo -e "kraken2.version\\t\$(kraken2 --version 2>&1 | head -1)"
        echo -e "bracken.version\\t\$(bracken -v 2>&1 | head -1)"
        echo -e "python.version\\t\$(python --version 2>&1)"
        echo -e "pandas.version\\t\$(python -c 'import pandas; print(pandas.__version__)')"
        echo -e "numpy.version\\t\$(python -c 'import numpy; print(numpy.__version__)')"
        echo -e "container\\t${task.container ?: 'not recorded (no container configured for this profile)'}"
        echo -e "param.phanta_confidence\\t${params.phanta_confidence}"
        echo -e "param.phanta_read_length\\t${params.phanta_read_length}"
        echo -e "param.phanta_filter_thresh\\t${params.phanta_filter_thresh}"
        echo -e "param.phanta_cov_thresh_bacterial\\t${params.phanta_cov_thresh_bacterial}"
        echo -e "param.phanta_cov_thresh_viral\\t${params.phanta_cov_thresh_viral}"
        echo -e "param.phanta_cov_thresh_arc\\t${params.phanta_cov_thresh_arc}"
        echo -e "param.phanta_cov_thresh_euk\\t${params.phanta_cov_thresh_euk}"
        echo -e "param.phanta_minimizer_thresh_bacterial\\t${params.phanta_minimizer_thresh_bacterial}"
        echo -e "param.phanta_minimizer_thresh_viral\\t${params.phanta_minimizer_thresh_viral}"
        echo -e "param.phanta_minimizer_thresh_arc\\t${params.phanta_minimizer_thresh_arc}"
        echo -e "param.phanta_minimizer_thresh_euk\\t${params.phanta_minimizer_thresh_euk}"
        echo -e "param.phanta_gzipped\\t${params.phanta_gzipped}"
    } > phanta_db_manifest.txt

    expected="${expected}"
    if [ -n "\$expected" ] && [ "\$expected" != "\$inspect_sha" ]; then
        echo "ERROR: this is not the database this run was configured for." >&2
        echo "  expected inspect.out sha256: \$expected" >&2
        echo "  actual:                      \$inspect_sha" >&2
        echo "  database:                    \$(realpath ${db})" >&2
        echo "Failing before any Kraken2 task rather than producing results" >&2
        echo "against an unknown database. Set --phanta_db_expected_inspect_sha256" >&2
        echo "to the actual value if the change is intended." >&2
        exit 1
    fi
    """
}
