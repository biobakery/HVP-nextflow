#!/usr/bin/env python3
"""Merge per-sample counts tables into one cohort table.

Adapted from Phanta's pipeline_scripts/merge.py at commit ecbfcf9.

Two changes:

  1. Sample IDs are passed explicitly and positionally alongside the input
     files, instead of being parsed out of each filename with
     basename.split(".krak.report")[0]. Recovering a sample ID from a filename
     is exactly what the port is meant to stop doing, and passing an ordered
     list is also what fixes the column order: Nextflow channel arrival order is
     nondeterministic, so the caller supplies samplesheet order and this script
     preserves it.

  2. The ProcessPoolExecutor is gone. It was only ever used with .map(), which
     yields in input order, so removing it cannot change the result -- and a
     pool of worker processes for reading a handful of small TSVs costs more
     than it saves. Determinism becomes structural rather than incidental.

Behaviour preserved exactly, including the part that is easy to mistake for a
bug: a sample whose table is empty -- which is what a Bracken-ineligible sample
produces -- is dropped entirely and gets NO column in the output. It is not a
column of zeros.
"""

import argparse
import os
import sys

import pandas as pd


def read_table(table, sample_id):
    if not os.path.exists(table):
        print("%s doesn't exists" % table)
        return None
    try:
        df = pd.read_csv(
            table, sep="\t",
            names=["Taxon_Lineage_with_Names", "Taxon_Lineage_with_IDs", sample_id]
        ).set_index(["Taxon_Lineage_with_Names", "Taxon_Lineage_with_IDs"])
    except pd.errors.EmptyDataError:
        print("%s is empty, pleae check" % table)
        return None
    return df if not df.empty else None


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-file", dest="input_file", nargs="+", required=True,
                   help="Per-sample tables, in samplesheet order.")
    p.add_argument("--sample-id", dest="sample_id", nargs="+", required=True,
                   help="Sample IDs, in the same order as --input-file.")
    p.add_argument("--output-file", dest="output_file", required=True)
    args = p.parse_args()

    if len(args.input_file) != len(args.sample_id):
        raise SystemExit(
            "ERROR: %d input files but %d sample IDs. They are positional and must "
            "correspond one to one." % (len(args.input_file), len(args.sample_id)))

    dfs = [df for df in (read_table(t, s) for t, s in zip(args.input_file, args.sample_id))
           if df is not None]

    if not dfs:
        raise SystemExit(
            "ERROR: every input table was empty or missing, so there is nothing to "
            "merge. If this is a real cohort, check whether every sample failed Bracken.")

    df_all = pd.concat(dfs, axis=1).fillna(0).reset_index()
    df_all.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
