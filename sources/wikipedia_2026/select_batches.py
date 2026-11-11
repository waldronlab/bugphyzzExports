#!/usr/bin/env python3 -I
"""Batches for a second extraction pass: only the taxa whose text matches a pattern.

Usage: select_batches.py <records.jsonl> <batches_dir> <regex>
Writes <batches_dir>/batch_NNNN.json (50 taxa each) in the format collect.py writes.
"""
import sys, os, json, re
records_file, out_dir, pattern = sys.argv[1:4]
rx = re.compile(pattern, re.I)
keep = []
for line in open(records_file):
    r = json.loads(line)
    if rx.search(r["text"]):
        keep.append({"taxid": r["taxid"], "name": r["name"], "rank": r["rank"], "text": r["text"]})
os.makedirs(out_dir, exist_ok=True)
for i in range(0, len(keep), 50):
    json.dump(keep[i:i + 50], open(os.path.join(out_dir, f"batch_{i // 50:04d}.json"), "w"), ensure_ascii=False)
print(f"{len(keep)} taxa in {(len(keep) + 49) // 50} batches")
