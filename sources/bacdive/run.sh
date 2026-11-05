#!/bin/sh
# Refresh the BacDive source: fetch every strain through the public API, then build the tables.
# Usage: sh run.sh <scratch_dir> <ncbi_taxdump_dir>
# The scratch dir holds the raw JSON (about 65 MB gzipped) and is not committed.
# Each fetch worker exits non-zero on any failed block; build.py also refuses
# to run on an incomplete scratch dir, so a partial table is never written.
set -e
raw="$1"; taxdump="$2"; here="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$raw"
pids=""
for range in "1 45000" "45001 90000" "90001 135000" "135001 180000"; do
    set -- $range
    python3 -I "$here/fetch.py" "$raw" "$2" 50 "$1" > "$raw/fetch_$1.log" 2>&1 &
    pids="$pids $!"
done
status=0
for pid in $pids; do wait "$pid" || status=1; done
if [ "$status" -ne 0 ]; then echo "a fetch worker failed; see $raw/fetch_*.log; re-run to resume" >&2; exit 1; fi
python3 -I "$here/build.py" "$raw" "$taxdump" "$here"
