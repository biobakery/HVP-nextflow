#!/usr/bin/env python3
"""Filter false-positive species out of a Kraken2 minimizer report.

Adapted from Phanta's pipeline_scripts/filter_kraken_reports.py at commit
ecbfcf9. Two things changed, and nothing else:

  1. Arguments are named, and every output path is explicit. Upstream takes ten
     bare positional arguments and derives its three output paths from the input
     path, which makes task isolation impossible. The positional order is also a
     trap: upstream's thresholds are not grouped by domain but as bacterial and
     viral coverage, then bacterial and viral minimizers, then archaeal and
     eukaryotic coverage, then archaeal and eukaryotic minimizers. Transcribing
     them in the intuitive order silently swaps thresholds, which is exactly the
     kind of error named arguments remove.

  2. The domain-rank vocabulary is detected from the database instead of being
     hard-coded to the string 'superkingdom'. See below.

The filtering logic itself is transcribed unchanged, including the round-trip
through the .species file (upstream writes it, then re-reads and re-parses it as
text, so the float formatting of the 'cov' column is load-bearing) and the
dictionary iteration order that decides row order in the decisions file.

## The rank vocabulary

Upstream asks the taxonomy for the ancestor at rank 'superkingdom'. NCBI has
since renamed that rank to 'domain', and moved Viruses to 'acellular root'. A
database built from a current NCBI taxonomy therefore answers "no such rank" for
every species -- and upstream's response to that is the catch-all `else` branch
below, which KEEPS the species. Filtering silently does nothing, the run exits
0, and the output tables look entirely normal.

So the vocabulary is detected once, from the ranks actually present in
nodes.dmp, and a database that has no domain-level vocabulary at all is a hard
error rather than a silent pass.

Per-taxon lookup failures are deliberately NOT fatal: a real database contains
taxa that legitimately sit outside the four domains ('other sequences', for
one), upstream keeps them, and making that fatal would break real runs. They are
counted and reported on stderr instead.

On a database whose nodes.dmp says 'superkingdom' -- which is what every
released Phanta database ships -- the detection selects 'superkingdom' and this
script is byte-for-byte equivalent to upstream.
"""

import argparse
import sys

import numpy as np
import pandas as pd

# Taxids of the four domains. Upstream compares against these literals, and they
# are stable across the superkingdom -> domain rename.
BACTERIA, ARCHAEA, EUKARYOTA, VIRUSES = '2', '2157', '2759', '10239'

# Rank spellings that denote a domain, newest NCBI first within each vocabulary.
# 'acellular root' and 'realm' are how current NCBI ranks Viruses, which has no
# 'domain' node -- without them a modern-taxonomy database would apply the wrong
# thresholds to every virus.
LEGACY_VOCABULARY = {'superkingdom'}
MODERN_VOCABULARY = {'domain', 'acellular root', 'realm'}

NOT_FOUND = 'error - taxid above desired rank, or not annotated at desired rank'


def make_dicts(nodes_file):
    """child -> parent and taxid -> rank, straight from nodes.dmp."""
    child_to_parent = {}
    taxid_to_rank = {}
    with open(nodes_file) as infile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            child, parent, rank = line[0], line[2], line[4]
            child_to_parent[child] = parent
            taxid_to_rank[child] = rank
    return child_to_parent, taxid_to_rank


def detect_domain_ranks(taxid_to_rank):
    """Which spelling of the domain rank does this database use?

    Returns the set of rank strings to treat as domain-level. Raises if the
    database has none of them, which means the taxonomy cannot answer the only
    question the thresholds depend on.
    """
    present = set(taxid_to_rank.values())

    if present & LEGACY_VOCABULARY:
        return LEGACY_VOCABULARY, 'superkingdom'

    modern = present & MODERN_VOCABULARY
    if modern:
        sys.stderr.write(
            "phanta_filter_kraken_reports: this database uses the modern NCBI rank "
            "vocabulary (%s) rather than 'superkingdom'. Falling back to it. Note that "
            "upstream Phanta would silently keep every species against this database, so "
            "results here will NOT match an upstream run on the same input.\n"
            % ', '.join(sorted(modern))
        )
        return MODERN_VOCABULARY, ', '.join(sorted(modern))

    raise SystemExit(
        "ERROR: no domain-level rank found in the database taxonomy.\n"
        "  Looked for: %s\n"
        "  The thresholds in this step are chosen per domain, so without one of these\n"
        "  ranks every species would fall through to the catch-all branch and be kept,\n"
        "  which disables filtering entirely while still exiting 0.\n"
        "  Check that the database's taxonomy/nodes.dmp is the one it shipped with."
        % ', '.join(sorted(LEGACY_VOCABULARY | MODERN_VOCABULARY))
    )


def taxid_to_desired_rank(taxid, desired_ranks, child_parent, taxid_rank):
    """Walk up the lineage to the first ancestor at one of `desired_ranks`.

    `desired_ranks` is a set rather than upstream's single string, which is the
    whole of the vocabulary change. Behaviour is otherwise identical, including
    returning the NOT_FOUND sentinel string rather than raising.
    """
    if taxid not in taxid_rank:
        return NOT_FOUND
    if taxid_rank[taxid] in desired_ranks:
        return taxid
    if taxid == '0':
        return 'unclassified'

    child, parent = taxid, None
    while not parent == '1':
        parent = child_parent.get(child)
        if parent is None:
            return NOT_FOUND
        if taxid_rank.get(parent) in desired_ranks:
            return parent
        child = parent
    return NOT_FOUND


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--krak-report', required=True)
    p.add_argument('--db', required=True)
    p.add_argument('--cov-thresh-bacterial', type=float, required=True)
    p.add_argument('--cov-thresh-viral', type=float, required=True)
    p.add_argument('--minimizer-thresh-bacterial', type=int, required=True)
    p.add_argument('--minimizer-thresh-viral', type=int, required=True)
    p.add_argument('--cov-thresh-arc', type=float, required=True)
    p.add_argument('--cov-thresh-euk', type=float, required=True)
    p.add_argument('--minimizer-thresh-arc', type=int, required=True)
    p.add_argument('--minimizer-thresh-euk', type=int, required=True)
    p.add_argument('--out-species', required=True)
    p.add_argument('--out-filtered', required=True)
    p.add_argument('--out-decisions', required=True)
    args = p.parse_args()

    # ── STEP THREE: species-level view of the report, with coverage ──────────
    kraken_file = pd.read_csv(
        args.krak_report, sep="\t",
        names=['fraction', 'fragments', 'assigned', 'minimizers', 'uniqminimizers',
               'classification_rank', 'ncbi_taxa', 'scientific name'])
    inspect_file = pd.read_csv(
        args.db + '/inspect.out', sep="\t",
        names=['frac', 'minimizers_clade', 'minimizers_taxa', 'rank', 'ncbi_taxonomy', 'sci_name'])

    kraken_file = kraken_file.merge(inspect_file, left_on='ncbi_taxa', right_on='ncbi_taxonomy')
    kraken_file['cov'] = kraken_file['uniqminimizers'] / kraken_file['minimizers_taxa']
    # Fewer than 5 minimizers on a taxon makes its coverage meaningless.
    kraken_file.loc[kraken_file.minimizers_taxa < 5, 'cov'] = np.nan

    species_kraken = kraken_file.copy()[kraken_file['classification_rank'].str.startswith('S', na=False)]

    child_parent, taxid_rank = make_dicts(args.db + '/taxonomy/nodes.dmp')
    domain_ranks, vocabulary = detect_domain_ranks(taxid_rank)

    species_kraken['species_level_taxa'] = species_kraken.apply(
        lambda x: taxid_to_desired_rank(str(x['ncbi_taxonomy']), {'species'}, child_parent, taxid_rank),
        axis=1)
    species_kraken.to_csv(args.out_species, sep="\t", index=False)

    # ── STEP FOUR: read that file back ───────────────────────────────────────
    # Upstream re-reads the file it just wrote rather than reusing the frame, so
    # the values here are the *round-tripped* ones. Keep it: changing it would
    # change float formatting and with it the decisions file.
    species_to_superkingdom = {}
    species_genome_info = {}
    unresolved_domains = 0

    with open(args.out_species) as infile:
        header = infile.readline().rstrip('\n').split('\t')
        cols = {name: i for i, name in enumerate(header)}
        for needed in ('uniqminimizers', 'cov', 'species_level_taxa', 'rank', 'ncbi_taxa'):
            if needed not in cols:
                raise SystemExit("ERROR: column '%s' missing from %s" % (needed, args.out_species))

        for line in infile:
            line = line.rstrip('\n').split('\t')
            unique_minimizers = int(line[cols['uniqminimizers']])
            cov = line[cols['cov']]
            species_taxid = line[cols['species_level_taxa']]
            rank = line[cols['rank']]
            taxid = line[cols['ncbi_taxa']]

            cov = 0 if cov == '' else float(cov)

            if rank == 'S':
                superkingdom = taxid_to_desired_rank(species_taxid, domain_ranks, child_parent, taxid_rank)
                if superkingdom == NOT_FOUND:
                    unresolved_domains += 1
                species_to_superkingdom[species_taxid] = superkingdom

            species_genome_info.setdefault(species_taxid, []).append(
                [taxid, cov, rank, unique_minimizers])

    if unresolved_domains:
        # Not fatal: upstream keeps these, and a real database legitimately has
        # some. Loud enough to notice if the count is suspiciously large.
        sys.stderr.write(
            "phanta_filter_kraken_reports: %d species could not be resolved to a domain "
            "(vocabulary: %s) and were kept unfiltered, as upstream does.\n"
            % (unresolved_domains, vocabulary))

    # ── STEP FIVE: decide what to keep ───────────────────────────────────────
    species_to_keep = set()

    with open(args.out_decisions, 'w') as outfile:
        outfile.write('\t'.join(
            ['species_taxid', 'superkingdom', 'max_cov', 'max_uniq_minimizers', 'kept']) + '\n')

        for species in species_genome_info:
            lines_to_consider = species_genome_info[species]
            relevant_coverages, relevant_unique_minimizer_vals = [], []

            # Only the LOWEST-rank entry under each species carries a usable
            # coverage: an S line followed by S1 lines is a summary of them.
            for i in range(len(lines_to_consider)):
                if i == (len(lines_to_consider) - 1):
                    relevant_coverages.append(lines_to_consider[i][1])
                    relevant_unique_minimizer_vals.append(lines_to_consider[i][3])
                else:
                    this_rank = lines_to_consider[i][2]
                    next_rank = lines_to_consider[i + 1][2]
                    if this_rank == 'S':
                        continue
                    if int(this_rank[1:]) < int(next_rank[1:]):
                        continue
                    relevant_coverages.append(lines_to_consider[i][1])
                    relevant_unique_minimizer_vals.append(lines_to_consider[i][3])

            max_cov = max(relevant_coverages)
            max_minimizers = max(relevant_unique_minimizer_vals)
            superkingdom = species_to_superkingdom[species]

            if superkingdom == BACTERIA:
                if max_cov >= args.cov_thresh_bacterial and max_minimizers >= args.minimizer_thresh_bacterial:
                    species_to_keep.add(species)
            elif superkingdom == ARCHAEA:
                if max_cov >= args.cov_thresh_arc and max_minimizers >= args.minimizer_thresh_arc:
                    species_to_keep.add(species)
            elif superkingdom == EUKARYOTA:
                if max_cov >= args.cov_thresh_euk and max_minimizers >= args.minimizer_thresh_euk:
                    species_to_keep.add(species)
            elif superkingdom == VIRUSES:
                if max_cov >= args.cov_thresh_viral and max_minimizers >= args.minimizer_thresh_viral:
                    species_to_keep.add(species)
            else:
                # Not one of the four domains. Upstream keeps it; so do we.
                species_to_keep.add(species)

            keep = 'True' if species in species_to_keep else 'False'
            outfile.write('\t'.join(
                [species, superkingdom, str(max_cov), str(max_minimizers), keep]) + '\n')

    # ── STEP SIX: rewrite the report without the rejected species ────────────
    with open(args.krak_report) as infile, open(args.out_filtered, 'w') as outfile:
        for line in infile:
            line = line.rstrip('\n').split('\t')
            if line[5].startswith('S'):
                species_taxid = taxid_to_desired_rank(line[6], {'species'}, child_parent, taxid_rank)
                if species_taxid in species_to_keep:
                    outfile.write('\t'.join(line) + '\n')
            else:
                outfile.write('\t'.join(line) + '\n')


if __name__ == '__main__':
    main()
