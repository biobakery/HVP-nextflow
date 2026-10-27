#!/usr/bin/env python3
"""Merge per-sample scaled Bracken tables into the relative-abundance table.

Adapted from Phanta's pipeline_scripts/merging_scaled_bracken.py at commit
ecbfcf9. Same two changes as phanta_merge.py: explicit ordered sample IDs
instead of filename parsing, and the ProcessPoolExecutor removed (it was used
only via .map(), which preserves input order, so the result is unchanged).

Upstream also ships merging_bracken_tables.R and merging_scaled_bracken.R.
Neither is wired into the Snakefile, and neither is ported.
"""

import argparse
import os

import pandas as pd


def read_scaled_table(table, sample_id):
    if not os.path.exists(table):
        print("%s doesn't exists" % table)
        return None
    try:
        df = pd.read_csv(table, sep="\t", dtype={'taxonomy_id': 'string'}) \
            .loc[:, ["taxonomy_id", "rel_taxon_abundance"]] \
            .rename(columns={"taxonomy_id": "TaxID", "rel_taxon_abundance": sample_id}) \
            .set_index("TaxID")
    except pd.errors.EmptyDataError:
        print("%s is empty, pleae check" % table)
        return None
    return df if not df.empty else None


def taxid_to_lineages(taxid_to_rank, taxid_to_name, child_to_parent, taxid):
    if taxid == '0':
        return 'unclassified'

    rank = taxid_to_rank[taxid]
    lineage_taxids = rank + '_' + taxid
    lineage_names = rank + '_' + taxid_to_name[taxid]
    child, parent = taxid, None

    while not parent == '1':
        parent = child_to_parent[child]
        parent_rank = taxid_to_rank[parent]
        parent_name = taxid_to_name[parent]
        lineage_taxids = parent_rank + '_' + parent + '|' + lineage_taxids
        lineage_names = parent_rank + '_' + parent_name + '|' + lineage_names
        child = parent
    return lineage_names, lineage_taxids


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-file", dest="input_file", nargs="+", required=True,
                   help="Per-sample scaled tables, in samplesheet order.")
    p.add_argument("--sample-id", dest="sample_id", nargs="+", required=True,
                   help="Sample IDs, in the same order as --input-file.")
    p.add_argument("--output-file", dest="output_file", required=True)
    p.add_argument("--db", dest="db", required=True)
    args = p.parse_args()

    if len(args.input_file) != len(args.sample_id):
        raise SystemExit(
            "ERROR: %d input files but %d sample IDs. They are positional and must "
            "correspond one to one." % (len(args.input_file), len(args.sample_id)))

    child_to_parent = {}
    taxid_to_rank = {}
    with open(os.path.join(args.db, 'taxonomy/nodes.dmp')) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            child, parent, rank = line[0], line[2], line[4]
            child_to_parent[str(child)] = parent
            taxid_to_rank[str(child)] = rank

    taxid_to_name = {}
    with open(os.path.join(args.db, 'inspect.out')) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            taxid, name = line[4], line[5].strip()
            taxid_to_name[str(taxid)] = name

    dfs = [df for df in (read_scaled_table(t, s)
                         for t, s in zip(args.input_file, args.sample_id))
           if df is not None]

    if not dfs:
        raise SystemExit(
            "ERROR: every scaled table was empty or missing, so there is nothing to merge.")

    prof_df = pd.concat(dfs, axis=1).fillna(0).reset_index()

    samples_list = list(prof_df.columns)[1:]
    headers = ["Taxon_Lineage_with_Names", "Taxon_Lineage_with_IDs"]
    prof_df[headers] = prof_df.apply(
        lambda x: taxid_to_lineages(taxid_to_rank, taxid_to_name, child_to_parent, x["TaxID"]),
        axis=1, result_type="expand")
    prof_df = prof_df.loc[:, headers + samples_list]
    prof_df.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
