#!/usr/bin/env python3
"""Scale a Bracken table by genome length.

Adapted from Phanta's pipeline_scripts/scale_bracken.py at commit ecbfcf9. The
only change is named arguments and an explicit output path; every calculation is
transcribed unchanged.

Note on dtypes, because it bites: the merge below joins a string column against
whatever pandas inferred for species_level_taxa in the database's
library/species_genome_size.txt. A real Phanta database holds genomes that
resolve no further than genus, for which calculate_genome_size.py writes the
literal 'error - taxid above species', so that column is object dtype and the
merge is string-to-string. Against a database where every row resolves cleanly
the column is int64 and pandas raises "You are trying to merge on object and
int64 columns". That is a property of the database, not of pandas -- 1.5.3 and
2.3.3 behave identically -- which is why the tiny test database deliberately
carries a genus-tagged genome.
"""

import argparse

import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--db', required=True)
    p.add_argument('--bracken-report', required=True)
    p.add_argument('--filtering-decisions', required=True)
    p.add_argument('--read-length', type=int, required=True)
    p.add_argument('--paired', required=True, choices=['True', 'False'],
                   help="Upstream passes the Python-literal spelling; kept for parity.")
    p.add_argument('--out', required=True)
    args = p.parse_args()

    read_length = args.read_length
    if args.paired == 'True':
        read_length = read_length * 2

    bracken_df = pd.read_csv(args.bracken_report, sep="\t")
    bracken_df['taxonomy_id'] = bracken_df['taxonomy_id'].astype('str')
    kraken_df = pd.read_csv(args.filtering_decisions, sep="\t")
    kraken_df['species_taxid'] = kraken_df['species_taxid'].astype('str')
    length_df = pd.read_csv(args.db + '/library/species_genome_size.txt', sep="\t")

    tmp = pd.merge(bracken_df,
                   length_df[['species_level_taxa', 'max', 'mean', 'median']],
                   left_on='taxonomy_id', right_on='species_level_taxa')

    per_million_scaling = (tmp['new_est_reads'].sum() / 1000000)
    tmp['fraction_total_reads'] = tmp['new_est_reads'] / tmp['new_est_reads'].sum()
    # Weight by genome length, so a short genome is not penalised for yielding
    # fewer reads at the same abundance.
    tmp['tmp_scaling'] = tmp['fraction_total_reads'] * (tmp['median'].sum() / tmp['median'])
    tmp['rel_taxon_abundance'] = tmp['tmp_scaling'] / tmp['tmp_scaling'].sum()
    tmp['genetic_abundance'] = tmp['fraction_total_reads']
    tmp['depth'] = (tmp['new_est_reads'] * read_length) / tmp['median']
    tmp['reads_per_million_bases'] = (tmp['new_est_reads'] / (tmp['median'] / 1000000)).round(4)
    tmp['reads_per_million_reads'] = (tmp['new_est_reads'] / per_million_scaling).round(4)
    tmp['RPMPM'] = tmp['new_est_reads'] / ((tmp['median'] / 1000000) * per_million_scaling)
    tmp['copies_per_million_reads'] = tmp['depth'] / (per_million_scaling).round(4)
    tmp = tmp.merge(kraken_df, how='left', left_on='taxonomy_id', right_on='species_taxid')

    tmp.rename(columns={'median': 'median_genome_length', 'max_cov': 'breadth'}, inplace=True)
    processed = tmp.drop(['tmp_scaling', 'species_level_taxa', 'max', 'mean',
                          'species_taxid', 'superkingdom', 'max_uniq_minimizers', 'kept'], axis=1)
    processed.to_csv(args.out, sep="\t", index=False)


if __name__ == '__main__':
    main()
