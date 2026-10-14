#!/usr/bin/env bash
#
# Build the Phanta tiny test database, and the synthetic reads that go with it.
#
# Run this INSIDE the pinned Phanta image, so that kraken2, bracken, dustmasker,
# pandas and biopython are the exact versions the baselines were generated with:
#
#   docker run --rm --platform linux/amd64 \
#     -v "$PWD":/w -w /w phanta:0.0.5-dev \
#     bash tests/fixtures/phanta/tiny_db/build.sh -o /w/build/tiny
#
# Why this database exists at all: the real Phanta database is 20-25 GB, which
# is too much to stage for every test run at either site, and it contains no
# sample that Bracken refuses to process — so the zero-byte-substitute path,
# which is load-bearing in the port, would never be exercised. See
# docs/phanta_migration_decisions.md section 8.
#
# Two things make this reproducible rather than an opaque blob:
#
#  1. The four genomes are pinned by VERSIONED accession and verified by SHA256
#     against accessions.tsv.
#  2. The taxonomy is NOT downloaded. taxonomy/nodes.dmp and names.dmp are
#     committed, hand-authored files. Appendix D of the spec says to run
#     `kraken2-build --download-taxonomy`; do not. NCBI's dump is a moving
#     target, so a download would make this database differ run to run, and —
#     worse — NCBI has renamed the rank `superkingdom` to `domain`, while
#     filter_kraken_reports.py looks up the literal string `superkingdom`. With
#     a freshly downloaded taxonomy every lookup falls through to the catch-all
#     `else` on line 248 that KEEPS the species, so filtering silently does
#     nothing and the output still looks entirely normal. The committed taxonomy
#     uses `superkingdom`, matching what the real Phanta database ships.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FIXTURES="$(dirname "$HERE")"
OUT=""
RECORD=0

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o|--out) OUT="$2"; shift 2 ;;
    --record-checksums) RECORD=1; shift ;;   # first-run bootstrap only
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$OUT" ]] || { echo "error: -o OUTDIR is required" >&2; exit 2; }

for tool in kraken2-build kraken2-inspect bracken-build python curl sha256sum; do
  command -v "$tool" >/dev/null || { echo "error: $tool not on PATH — run this inside the Phanta image" >&2; exit 3; }
done

DB="$OUT/tiny_db"
SRC="$OUT/sources"
READS="$OUT/reads"
mkdir -p "$DB/taxonomy" "$SRC" "$READS"

EUTILS="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# ── 1. fetch and verify the four genomes ─────────────────────────────────────
# Tab-separated, comments stripped: class, accession, seq_taxid, species_taxid,
# organism, sha256.
while IFS=$'\t' read -r class accession seq_taxid species_taxid organism sha; do
  [[ -z "${class:-}" || "$class" == \#* ]] && continue

  raw="$SRC/$accession.raw.fna"
  if [[ ! -s "$raw" ]]; then
    echo "fetching $accession ($organism)"
    curl -sS --retry 3 --max-time 120 \
      "$EUTILS?db=nuccore&id=$accession&rettype=fasta&retmode=text" -o "$raw"
    [[ -s "$raw" ]] || { echo "error: empty download for $accession" >&2; exit 4; }
  fi

  got="$(sha256sum "$raw" | cut -d' ' -f1)"
  if [[ "$RECORD" == 1 ]]; then
    echo "RECORD	$accession	$got"
  elif [[ "$sha" == "PENDING" ]]; then
    echo "error: accessions.tsv still says PENDING for $accession." >&2
    echo "       Run once with --record-checksums and paste the values in." >&2
    exit 5
  elif [[ "$got" != "$sha" ]]; then
    echo "error: checksum mismatch for $accession" >&2
    echo "       expected $sha" >&2
    echo "       got      $got" >&2
    echo "       NCBI changed the record, or the download truncated. Do not" >&2
    echo "       proceed: every committed baseline assumes these exact bytes." >&2
    exit 6
  fi

  # Retag the header for Kraken2. The taxid must sit before any whitespace,
  # because calculate_genome_size.py parses it out of Biopython's rec.id, which
  # stops at the first space.
  awk -v acc="$accession" -v taxid="$seq_taxid" \
    'NR==1 {print ">" acc "|kraken:taxid|" taxid; next} {print}' \
    "$raw" > "$SRC/$accession.fna"
done < <(grep -v '^#' "$HERE/accessions.tsv" | grep -v '^[[:space:]]*$')

[[ "$RECORD" == 1 ]] && { echo "checksums recorded above; nothing else was built"; exit 0; }

# ── 2. taxonomy: copy the committed dumps in, never download ─────────────────
cp "$HERE/taxonomy/nodes.dmp" "$HERE/taxonomy/names.dmp" "$DB/taxonomy/"

# ── 3. Kraken2 index ─────────────────────────────────────────────────────────
for fna in "$SRC"/*.fna; do
  [[ "$fna" == *.raw.fna ]] && continue
  kraken2-build --add-to-library "$fna" --db "$DB"
done
kraken2-build --build --db "$DB"

# inspect.out is what filter_kraken_reports.py and prep_to_merge_counts.py read,
# and it is the database's identity fingerprint in phanta_db_manifest.txt.
kraken2-inspect --db "$DB" > "$DB/inspect.out"

# ── 4. Bracken distribution ──────────────────────────────────────────────────
# -l 100 must match READ_LEN in make_reads.py and params.phanta_read_length.
bracken-build -d "$DB" -t 4 -k 35 -l 100

# ── 5. per-species genome sizes, read by scale_bracken.py ────────────────────
python "$HERE/calculate_genome_size.py" "$DB"

# ── 6. synthetic reads ───────────────────────────────────────────────────────
python "$FIXTURES/tiny_inputs/make_reads.py" --sources "$SRC" --outdir "$READS"

# ── 7. check every file the spec's "Required database files" table lists ─────
missing=0
for f in hash.k2d taxo.k2d opts.k2d seqid2taxid.map database100mers.kmer_distrib \
         inspect.out taxonomy/nodes.dmp taxonomy/names.dmp library/species_genome_size.txt; do
  if [[ ! -s "$DB/$f" ]]; then echo "MISSING: $f" >&2; missing=1; fi
done
[[ "$missing" == 0 ]] || { echo "error: database is incomplete" >&2; exit 7; }

echo
echo "tiny database built: $DB"
du -sh "$DB"
echo "inspect.out sha256: $(sha256sum "$DB/inspect.out" | cut -d' ' -f1)"
