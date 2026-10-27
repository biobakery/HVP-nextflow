#!/usr/bin/env python3
"""Decide whether Bracken can run on a filtered Kraken report.

Adapted from Phanta's pipeline_scripts/process_filtered_kraken.py at commit
ecbfcf9.

Bracken errors out when no species in the report clears the read threshold.
Upstream protects against that by writing three zero-byte `.temp` files whose
names the next two rules test for with `[ -f ... ]`, copying them into place
instead of running the tool. Those files are substitute OUTPUTS, not status
flags -- an ineligible sample's Bracken table, Kraken-style species report and
scaled table are all genuinely empty files, and the fixtures contain them.

Two changes:

  1. Explicit output paths, so nothing is derived from the input path and the
     task can be isolated.
  2. An explicit status file, so the pipeline can branch on a recorded decision
     and a reason rather than on the existence of a canary file. The zero-byte
     products are still written, byte-identically, because they are part of the
     observable contract.

Upstream also loads nodes.dmp here and never uses it; that dead read is dropped,
which is why this takes no --db.
"""

import argparse


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--filtered-report', required=True)
    p.add_argument('--filter-thresh', type=int, required=True,
                   help='Same value as Bracken -t. One parameter drives both.')
    p.add_argument('--sample', required=True)
    p.add_argument('--out-status', required=True)
    p.add_argument('--out-bracken', required=True,
                   help='Zero-byte substitute for the Bracken table; written only when ineligible.')
    p.add_argument('--out-species-report', required=True,
                   help='Zero-byte substitute for the Kraken-style species report.')
    p.add_argument('--out-scaled', required=True,
                   help='Zero-byte substitute for the scaled Bracken table.')
    args = p.parse_args()

    # Eligible as soon as ONE species line clears the threshold.
    eligible = False
    best = 0
    with open(args.filtered_report) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            rank, assigned = line[5], int(line[1])
            if rank == 'S':
                best = max(best, assigned)
                if assigned >= args.filter_thresh:
                    eligible = True
                    break

    with open(args.out_status, 'w') as status:
        if eligible:
            status.write('%s\teligible\t-\n' % args.sample)
        else:
            status.write('%s\tineligible\tno species line reaches filter_thresh=%d '
                         '(highest was %d assigned reads)\n' % (args.sample, args.filter_thresh, best))

    if not eligible:
        # Genuinely empty, not absent: downstream copies these into place.
        for path in (args.out_bracken, args.out_species_report, args.out_scaled):
            open(path, 'w').close()


if __name__ == '__main__':
    main()
