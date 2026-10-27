#!/usr/bin/env python3
"""Normalise one sample's counts table by its Bracken-classified read count.

Adapted from Phanta's pipeline_scripts/prep_to_merge_normed.py at commit
ecbfcf9. This is the one script whose contract genuinely changes.

Upstream scans a shared classification directory: it reads every sample name out
of total_reads.tsv and writes a .norm_brack beside each sample's counts table,
in one invocation. That cannot be isolated to a task, so this processes ONE
sample per invocation, taking that sample's counts table explicitly and looking
up only its own row in the cohort total_reads.tsv.

Two edge cases carried over deliberately:
  * an empty Bracken-classified read count is treated as 0
  * a failed sample produces an empty output file, not an absent one

The arithmetic is reproduced exactly: round(reads / classified, 10) formatted
with plain str(). It is NOT done through pandas, and must not be -- the string
formatting is what lands in the published table.
"""

import argparse


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--counts-table', required=True)
    p.add_argument('--total-reads', required=True, help='The cohort total_reads.tsv.')
    p.add_argument('--failed-file', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--out', required=True)
    args = p.parse_args()

    # Column 4 (Assigned_Step_Three) is the Bracken-classified read count.
    brack_reads = None
    with open(args.total_reads) as infile:
        infile.readline()  # header
        for line in infile:
            line = line.rstrip('\n').split('\t')
            if line[0] == args.sample:
                brack_reads = line[3]
                break

    if brack_reads is None:
        raise SystemExit(
            "ERROR: sample '%s' has no row in %s. The cohort total_reads.tsv must "
            "cover every sample being normalised." % (args.sample, args.total_reads))

    with open(args.failed_file) as infile:
        failed_samples = {line.rstrip('\n') for line in infile}

    with open(args.counts_table) as infile, open(args.out, 'w') as outfile:
        if args.sample in failed_samples:
            pass
        else:
            brack_classified_reads = 0 if brack_reads == '' else int(brack_reads)
            for line in infile:
                line = line.rstrip('\n').split('\t')
                lin_names, lin_taxids, num_reads = line[0], line[1], float(line[2])
                frac = str(round(num_reads / brack_classified_reads, 10))
                outfile.write('\t'.join([lin_names, lin_taxids, frac]) + '\n')


if __name__ == '__main__':
    main()
