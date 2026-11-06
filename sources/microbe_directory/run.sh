#!/bin/sh
# Refresh The Microbe Directory source: download the pinned table and NCBI names, then build.
# Usage: sh run.sh <scratch_dir>
# The scratch dir holds the downloads (about 170 MB) and is not committed.
set -e
raw="$1"; here="$(cd "$(dirname "$0")" && pwd)"
[ -n "$raw" ] || { echo "usage: sh run.sh <scratch_dir>" >&2; exit 1; }
mkdir -p "$raw"
# The Microbe Directory 2018 release (MIT licence); unchanged since 2018-01
MD_COMMIT=02fbf2633dfab0c418701197b6f869fa038d6776
# NCBI new_taxdump of the same date as the pinned ncbi-taxon-db release, for synonyms only
TAXDUMP=new_taxdump_2026-11-01.zip
[ -s "$raw/microbe-directory.csv" ] || curl -sSf -o "$raw/microbe-directory.csv" \
    "https://raw.githubusercontent.com/microbe-directory/microbe-directory/$MD_COMMIT/data/microbe-directory.csv"
[ -s "$raw/names.dmp" ] || {
    curl -sSf -o "$raw/$TAXDUMP" "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/$TAXDUMP"
    unzip -o -q "$raw/$TAXDUMP" names.dmp -d "$raw"
}
sh "$here/python.sh" "$here/build.py" "$raw/microbe-directory.csv" "$raw/names.dmp" "$here"
