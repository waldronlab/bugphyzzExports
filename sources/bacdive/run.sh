#!/bin/sh
# Refresh the BacDive source: fetch every strain through the public API, then build the tables.
# Usage: sh run.sh <scratch_dir> <ncbi_taxdump_dir>
# The scratch dir holds the raw JSON (about 1 GB gzipped) and is not committed.
set -e
raw="$1"; taxdump="$2"; here="$(cd "$(dirname "$0")" && pwd)"
for range in "1 45000" "45001 90000" "90001 135000" "135001 180000"; do
    set -- $range
    python3 -I "$here/fetch.py" "$raw" "$2" 50 "$1" > "$raw/fetch_$1.log" 2>&1 &
done
wait
python3 -I "$here/build.py" "$raw" "$taxdump" "$here"
