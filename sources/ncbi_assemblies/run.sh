#!/bin/sh
# Refresh the NCBI RefSeq assembly source: download the assembly summary, then build the tables.
# Usage: sh run.sh <scratch_dir>
# The scratch dir holds assembly_summary_refseq.txt (about 250 MB) and is not committed.
set -e
raw="$1"; here="$(cd "$(dirname "$0")" && pwd)"
[ -n "$raw" ] || { echo "usage: sh run.sh <scratch_dir>" >&2; exit 1; }
mkdir -p "$raw"
summary="$raw/assembly_summary_refseq.txt"
# -R keeps the server's modification time, which is the date to record in the README
curl -sSf -R -o "$summary" https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
echo "assembly summary dated $(date -r "$summary" +%Y-%m-%d)"
sh "$here/python.sh" "$here/build.py" "$summary" "$here"
