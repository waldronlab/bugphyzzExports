#!/usr/bin/env python3
"""Fetch every BacDive strain record through the public API (no key needed).

Usage: python3 -I fetch.py <raw_dir> [max_id] [chunk] [start_id]
Run several copies with different start_id/max_id ranges to parallelise.
Writes <raw_dir>/<start>-<end>.json.gz, one file per chunk of ids, and skips
chunks already present, so it can be re-run to resume. Exits non-zero if any
chunk could not be fetched. Raw files are not committed; build.py turns them
into the tables.
"""
import sys, os, json, gzip, time, urllib.request, urllib.error
raw_dir = sys.argv[1]; max_id = int(sys.argv[2]) if len(sys.argv) > 2 else 180000
chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 50
start_id = int(sys.argv[4]) if len(sys.argv) > 4 else 1
os.makedirs(raw_dir, exist_ok=True)
BASE = "https://api.bacdive.dsmz.de/fetch/"
hdr = {"User-Agent": "bugphyzzExports-sources (https://github.com/waldronlab/bugphyzzExports)", "Accept": "application/json"}
total = 0; failed = 0
for start in range(start_id, max_id + 1, chunk):
    end = min(start + chunk - 1, max_id)
    out = os.path.join(raw_dir, f"{start:06d}-{end:06d}.json.gz")
    if os.path.exists(out):
        continue
    ids = ";".join(str(i) for i in range(start, end + 1))
    for attempt in range(6):
        try:
            req = urllib.request.Request(BASE + ids, headers=hdr)
            data = json.load(urllib.request.urlopen(req, timeout=180))
            break
        except urllib.error.HTTPError as e:
            if e.code == 404:  # the API answers 404 when none of the ids exist
                data = {"results": {}}; break
            time.sleep(15 * (attempt + 1)); data = None
        except Exception as e:  # truncated chunked responses, timeouts, bad JSON
            time.sleep(15 * (attempt + 1)); data = None
    if data is None:
        print(f"{start}-{end}: FAILED", flush=True); failed += 1; continue
    res = data.get("results") or {}
    n = len(res); total += n
    with gzip.open(out, "wt") as fh:
        json.dump(res, fh)
    ## BacDive ids have long gaps, so empty chunks are normal; stop only at max_id.
    if start % 5000 == 1:
        print(f"{start}-{end}: {n} records (total {total})", flush=True)
    time.sleep(0.2)
print(f"done; total records {total}; failed chunks {failed}", flush=True)
sys.exit(1 if failed else 0)
