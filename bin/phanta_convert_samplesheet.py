#!/usr/bin/env python3
"""Convert Phanta's tab-delimited sample file into this pipeline's CSV.

Upstream Phanta takes a headerless, tab-delimited file:

    sample_name<TAB>read1_file[<TAB>read2_file]

This pipeline takes a CSV with a header:

    sample,fastq_1,fastq_2

Two formats exist because cohort aggregation needs an ordered sample list and
the rest of this repository works from --readsdir plus filename globbing; the
CSV is the new input model, and this exists so an existing Phanta sample file
does not have to be retyped. Row ORDER is preserved, because it decides the
column order of the merged tables.

Usage:  phanta_convert_samplesheet.py --input samp_file.txt --output samples.csv
"""

import argparse
import csv
import os
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help="Phanta's tab-delimited sample file")
    p.add_argument('--output', required=True, help='CSV to write')
    p.add_argument('--check-exists', action='store_true',
                   help='Also verify every read file is present (off by default, '
                        'because upstream sample files often carry placeholder paths)')
    args = p.parse_args()

    rows = []
    seen = {}
    with open(args.input) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue

            cells = [c.strip() for c in line.split('\t')]
            if len(cells) < 2:
                raise SystemExit(
                    "ERROR: %s line %d: expected at least 'sample<TAB>read1'.\n"
                    "  Found: %r" % (args.input, lineno, line))

            sample, fastq_1 = cells[0], cells[1]
            fastq_2 = cells[2] if len(cells) > 2 else ''

            if sample in seen:
                raise SystemExit(
                    "ERROR: %s line %d: duplicate sample '%s', already seen on line %d.\n"
                    "  Sample names become column names in the merged tables and must "
                    "be unique." % (args.input, lineno, sample, seen[sample]))
            seen[sample] = lineno

            if '/full/path/to' in fastq_1:
                sys.stderr.write(
                    "warning: %s line %d still has an upstream placeholder path (%s).\n"
                    % (args.input, lineno, fastq_1))

            if args.check_exists:
                for path in filter(None, (fastq_1, fastq_2)):
                    if not os.path.exists(path):
                        raise SystemExit("ERROR: %s line %d: no such file: %s"
                                         % (args.input, lineno, path))

            rows.append((sample, fastq_1, fastq_2))

    if not rows:
        raise SystemExit("ERROR: %s has no sample rows." % args.input)

    paired = [r for r in rows if r[2]]
    if paired and len(paired) != len(rows):
        raise SystemExit(
            "ERROR: %s mixes paired-end and single-end samples (%d paired, %d single).\n"
            "  Phanta applies one library layout to the whole cohort."
            % (args.input, len(paired), len(rows) - len(paired)))

    with open(args.output, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['sample', 'fastq_1', 'fastq_2'])
        writer.writerows(rows)

    print("%s: %d samples, %s -> %s"
          % (args.input, len(rows), 'paired-end' if paired else 'single-end', args.output))


if __name__ == '__main__':
    main()
