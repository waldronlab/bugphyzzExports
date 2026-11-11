#!/usr/bin/env python3 -I
"""Collect the English Wikipedia introductions that taxoniq carries for Bacteria and Archaea.

Usage: collect.py <ncbi_taxon_db_package_dir> <scratch_dir>
<ncbi_taxon_db_package_dir> is the directory that contains the `ncbi_taxon_db`
package of a taxoniq build with the Wikipedia stage (description, en_wiki_title
and wikidata indexes); the taxonomy-only release wheels do not have them.
Writes <scratch_dir>/records.jsonl (one taxon per line) and
<scratch_dir>/batches/batch_NNNN.json (50 taxa each, the input of extract.py).
Neither is committed: the text is CC BY-SA.
"""
import sys, os, json, re, html, datetime
db_pkg_dir, scratch = sys.argv[1:3]
sys.path.insert(0, db_pkg_dir)
import taxoniq, ncbi_taxon_db
assert os.path.dirname(ncbi_taxon_db.__file__).startswith(os.path.abspath(db_pkg_dir)), ncbi_taxon_db.__file__
BATCH = 50
PROKARYOTES = {2, 2157}   # Bacteria, Archaea
MIN_CHARS = 40            # shorter than this is a redirect stub or an empty extract

def plain(extract):
    """The extract is HTML (<p>, <i>, <b>, empty <span> left from pronunciation templates)."""
    t = re.sub(r"</p>|<br\s*/?>", "\n", extract)
    t = re.sub(r"<[^>]+>", "", t)
    t = html.unescape(t).replace("\xa0", " ")
    t = re.sub(r"\(\s*[;,]?\s*\)", "", t)          # "( )" left where a pronunciation was
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r"\s*\n\s*", "\n", t).strip()
    return t

probe = taxoniq.Taxon(2)
taxids = sorted(int(k) for k in probe._get_db("description_pos").keys())
records, skipped = [], {"not Bacteria or Archaea": 0, "taxid not in taxonomy": 0, "empty extract": 0}
for taxid in taxids:
    try:
        t = taxoniq.Taxon(taxid)
        lineage = {x.tax_id for x in t.lineage}
    except KeyError:
        skipped["taxid not in taxonomy"] += 1; continue
    if not lineage & PROKARYOTES:
        skipped["not Bacteria or Archaea"] += 1; continue
    text = plain(t.description)
    if len(text) < MIN_CHARS:
        skipped["empty extract"] += 1; continue
    try: wikidata = t.wikidata_id
    except KeyError: wikidata = ""
    try: title = t._get_str_attr("en_wiki_title")
    except (KeyError, taxoniq.NoValue): title = ""
    records.append({"taxid": str(taxid), "name": t.scientific_name, "rank": t.rank.name,
                    "domain": "Archaea" if 2157 in lineage else "Bacteria",
                    "wikidata_id": wikidata, "en_wiki_title": title, "text": text})

os.makedirs(os.path.join(scratch, "batches"), exist_ok=True)
with open(os.path.join(scratch, "records.jsonl"), "w") as fh:
    for r in records:
        fh.write(json.dumps(r, ensure_ascii=False) + "\n")
for i in range(0, len(records), BATCH):
    batch = [{k: r[k] for k in ("taxid", "name", "rank", "text")} for r in records[i:i + BATCH]]
    json.dump(batch, open(os.path.join(scratch, "batches", f"batch_{i // BATCH:04d}.json"), "w"), ensure_ascii=False)
ts = ncbi_taxon_db.db_timestamp
built = datetime.datetime.fromtimestamp(int(ts), datetime.timezone.utc).date().isoformat() if str(ts).isdigit() else ts
ranks = {}
for r in records: ranks[r["rank"]] = ranks.get(r["rank"], 0) + 1
print(f"taxoniq {taxoniq.__version__}, ncbi_taxon_db taxdump of {built}")
print(f"{len(taxids)} taxa with a description; kept {len(records)} Bacteria and Archaea in {(len(records) + BATCH - 1) // BATCH} batches")
print("ranks:", dict(sorted(ranks.items(), key=lambda kv: -kv[1])))
print("titles shared by more than one taxid:", len(records) - len({r["en_wiki_title"] for r in records}))
print("skipped:", skipped)
