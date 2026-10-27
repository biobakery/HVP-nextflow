#!/usr/bin/env python3
"""Turn a Kraken-style Bracken species report into a mergeable counts table.

Adapted from Phanta's pipeline_scripts/prep_to_merge_counts.py at commit
ecbfcf9. Named arguments and an explicit output path; the lineage construction
is transcribed unchanged.

Note which file this reads: the Kraken-STYLE report Bracken emits as a side
effect ({samp}.krak.report_bracken_species.filtered), not the .bracken table.
It takes the taxid from column 5 and the read count from column 2.
"""

import argparse


def read_nodes(nodes_file):
    child_to_parent = {}
    taxid_to_rank = {}
    with open(nodes_file) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            child, parent, rank = line[0], line[2], line[4]
            child_to_parent[child] = parent
            taxid_to_rank[child] = rank
    return child_to_parent, taxid_to_rank


def read_names(inspect_file):
    """taxid -> name, from inspect.out.

    inspect.out is the only place Phanta reads taxon names from, so an
    inspect.out whose name column is blank -- which is what Kraken2 produces from
    a names.dmp with the wrong line terminator -- yields empty lineage names
    everywhere, with no error.
    """
    taxid_to_name = {}
    with open(inspect_file) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            taxid, name = line[4], line[5].strip()
            taxid_to_name[taxid] = name
    return taxid_to_name


def taxid_to_lineages(taxid, child_to_parent, taxid_to_rank, taxid_to_name):
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
    return lineage_taxids, lineage_names


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--brack-report', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--failed-file', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    child_to_parent, taxid_to_rank = read_nodes(args.db + '/taxonomy/nodes.dmp')
    taxid_to_name = read_names(args.db + '/inspect.out')

    with open(args.failed_file) as infile:
        failed_samples = {line.rstrip('\n') for line in infile}

    with open(args.brack_report) as infile, open(args.out, 'w') as outfile:
        if args.sample in failed_samples:
            pass  # a failed sample yields an empty table, not an absent one
        else:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                taxid = line[4]
                lin_taxids, lin_names = taxid_to_lineages(
                    taxid, child_to_parent, taxid_to_rank, taxid_to_name)
                num_reads = line[1]
                outfile.write('\t'.join([lin_names, lin_taxids, num_reads]) + '\n')


if __name__ == '__main__':
    main()
