#!/usr/bin/env python3
"""Generate the synthetic read sets for the Phanta tiny-database tests.

Three samples, all paired-end, all drawn from the same four genomes the tiny
database is built from:

  sample_a, sample_b  ordinary cohort members. Enough reads per genome to clear
                      every cov_thresh_*/minimizer_thresh_* in Appendix B and
                      the filter_thresh=10 read floor, so Bracken runs normally.
  sample_c            deliberately Bracken-INELIGIBLE. Four read pairs, which is
                      below filter_thresh, so process_filtered_kraken.py writes
                      its zero-byte substitutes and Bracken never runs. This is
                      the only way to exercise the empty-product path and the
                      failed-sample list; Phanta's own test reads never trigger
                      it.

Determinism is the point of this script, not convenience. The outputs are
committed fixtures that published baselines are compared against byte for byte,
so every source of variation is nailed down: a fixed PRNG seed per sample, a
fixed iteration order, and gzip headers written with mtime=0 and no embedded
filename (Python's gzip module otherwise stamps the current time into every
file, which alone would make two runs differ).

Usage:  make_reads.py --sources DIR --outdir DIR
"""

import argparse
import gzip
import pathlib
import random

READ_LEN = 100          # must match the Bracken database: bracken-build -l 100
INSERT = 300            # outer distance between the two mates
QUAL = "I" * READ_LEN   # flat Q40; quality is never read by this pipeline

# sample -> genome basename -> number of read PAIRS. Hand-tuned so a/b clear the
# thresholds comfortably and c clears nothing. Changing any number here
# invalidates the committed baselines.
COMPOSITION = {
    "sample_a": {"NC_000908.2": 2000, "NC_000909.1": 600, "NC_001133.9": 400,
                 "NC_001422.1": 250, "NC_000912.1": 200},
    "sample_b": {"NC_000908.2": 1200, "NC_000909.1": 900, "NC_001133.9": 250,
                 "NC_001422.1": 60, "NC_000912.1": 150},
    "sample_c": {"NC_000908.2": 4},
}
SEEDS = {"sample_a": 20251014, "sample_b": 20251015, "sample_c": 20251016}

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(path):
    """Return the concatenated sequence of a single-record FASTA, uppercased."""
    parts = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts).upper()


def open_deterministic_gzip(path):
    """gzip.open, minus the two fields that would differ between runs."""
    raw = open(path, "wb")
    return gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sources", required=True, help="directory of fetched FASTAs")
    ap.add_argument("--outdir", required=True, help="where to write the .fq.gz files")
    args = ap.parse_args()

    sources = pathlib.Path(args.sources)
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genomes = {}
    for accession in sorted({a for comp in COMPOSITION.values() for a in comp}):
        genomes[accession] = read_fasta(sources / f"{accession}.fna")

    for sample in sorted(COMPOSITION):
        rng = random.Random(SEEDS[sample])
        r1_path = outdir / f"{sample}_R1.fq.gz"
        r2_path = outdir / f"{sample}_R2.fq.gz"

        with open_deterministic_gzip(r1_path) as r1, open_deterministic_gzip(r2_path) as r2:
            index = 0
            # sorted() so the interleaving of genomes is fixed, not dict-order dependent
            for accession in sorted(COMPOSITION[sample]):
                seq = genomes[accession]
                for _ in range(COMPOSITION[sample][accession]):
                    index += 1
                    start = rng.randrange(0, len(seq) - INSERT)
                    fwd = seq[start:start + READ_LEN]
                    rev = revcomp(seq[start + INSERT - READ_LEN:start + INSERT])
                    name = f"{sample}_{index:06d}"
                    r1.write(f"@{name}/1\n{fwd}\n+\n{QUAL}\n".encode())
                    r2.write(f"@{name}/2\n{rev}\n+\n{QUAL}\n".encode())

        print(f"{sample}: {index} pairs -> {r1_path.name}, {r2_path.name}")


if __name__ == "__main__":
    main()
