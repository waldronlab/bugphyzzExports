#!/usr/bin/env python3 -I
"""Turn the Gemini extraction (one JSON list per batch) into bugphyzz curation-sheet rows.

Usage: build_traits.py <extraction_dir> <genus_taxids.json> <out_dir> [<arrangement_dir>]
<arrangement_dir> holds the second pass (system_prompt_arrangement.txt): cell
arrangement and spore shape.
Writes one TSV per physiology in the Google Sheet column layout
(NCBI_ID, Genome_ID, Accession_ID, Taxon_name, Attribute, Attribute_value,
Attribute_source, Frequency, Note) plus overrides applied from overrides.tsv.
"""
import json, sys, os, glob, csv, collections
ext_dir, taxid_file, out_dir = sys.argv[1:4]
arr_dir = sys.argv[4] if len(sys.argv) > 4 else None
os.makedirs(out_dir, exist_ok=True)
tax = json.load(open(taxid_file))
dois = {d["doi"]: d for d in json.load(open(os.path.join(os.path.dirname(__file__), "dois.json")))}
SOURCE = "Bergeys_2026"
# extraction vocabulary -> bugphyzz attributes.tsv vocabulary (None = drop)
SHAPE = {"bacillus": "bacillus", "coccus": "coccus", "coccobacillus": "coccobacillus", "spirillum": "spirillum",
         "filamentous": "filamentous", "elliptic": "elliptic", "pleomorphic": "pleiomorphic cell", "vibrio": "vibrio", "other": None}
AERO = {"aerobic": "aerobic", "anaerobic": "anaerobic", "facultatively anaerobic": "facultatively anaerobic", "microaerophilic": "microaerophilic"}
FREQ = {"always", "usually", "sometimes"}
ARRANGEMENT = {"single", "paired cells", "cell chain", "cell cluster", "tetrads", "palisade", "rosettes", "v-formation", "branched", "filamentous"}
SPORE_SHAPE = {"elliptic", "coccus", "bacillus"}
UNIT = {"growth temperature": "Celcius", "gc content": "mol%", "width": "um", "length": "um"}
rows = collections.defaultdict(list); dropped = collections.Counter(); chapters = {}; arr = {}
for f in sorted(glob.glob(os.path.join(ext_dir, "batch_*.json"))):
    for g in json.load(open(f)):
        chapters[g["doi"]] = g
if arr_dir:
    for f in sorted(glob.glob(os.path.join(arr_dir, "batch_*.json"))):
        for g in json.load(open(f)):
            arr[g["doi"]] = g
# newest chapter per genus taxid (a revision .pub2, or a chapter on a genus NCBI
# has since merged, would otherwise give the taxid two sets of rows)
by_genus = collections.defaultdict(list)
def taxon(doi, g):
    # the map is keyed by the chapter title in dois.json; the model may echo a tidied name
    title = dois.get(doi, {}).get("genus")
    return tax.get(title) or tax.get(g["genus"]) or {}
for doi, g in chapters.items():
    key = taxon(doi, g).get("taxid") or ("title", g["genus"])
    by_genus[key].append((dois.get(doi, {}).get("year") or 0, doi))
# Per trait, the newest chapter wins where it states the trait; where it is
# silent the next-older chapter fills in (#102). `src` records which chapter
# each trait came from, for the Note column.
LIST_FIELDS = ("aerophilicity", "gram_stain", "shape", "motility", "spore_formation")
RANGE_FIELDS = ("growth_temperature_c", "gc_content_mol_percent", "cell_length_um", "cell_width_um")
def stated(rec, field):
    v = (rec or {}).get(field)
    if field in LIST_FIELDS:
        return bool(v)
    return isinstance(v, dict) and any(x is not None for x in v.values())
filled_from_older = collections.Counter()
for key, lst in sorted(by_genus.items(), key=lambda kv: str(kv[0])):
    order = [d for _, d in sorted(lst, reverse=True)]   # newest first
    doi = order[0]; g0 = chapters[doi]; t = taxon(doi, g0)
    if not t.get("taxid"):
        dropped["genus without NCBI taxid"] += 1; continue
    g, arr_rec, src = {"genus": g0["genus"]}, {}, {}
    for field in LIST_FIELDS + RANGE_FIELDS:
        for d in order:
            if stated(chapters[d], field):
                g[field] = chapters[d][field]; src[field] = d
                if d != doi: filled_from_older[field] += 1
                break
    for field in ("arrangement", "spore_shape"):
        for d in order:
            if (arr.get(d) or {}).get(field):
                arr_rec[field] = arr[d][field]; src[field] = d
                if d != doi: filled_from_older[field] += 1
                break
    base = {"NCBI_ID": t["taxid"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": t["clean"], "Attribute_source": SOURCE, "Note": doi}
    newest = doi
    def emit(phys, attr, val, freq, note=""):
        if freq not in FREQ: dropped[f"frequency {freq}"] += 1; return
        rows[phys].append({**base, "Attribute": attr, "Attribute_value": val, "Frequency": freq, "Note": (doi + ("; " + note if note else ""))})
    doi = src.get("aerophilicity", newest)
    for e in g.get("aerophilicity") or []:
        a = AERO.get(e["value"]); emit("aerophilicity", a, "TRUE", e["frequency"]) if a else dropped.update([f"aerophilicity {e['value']}"])
    doi = src.get("gram_stain", newest)
    for e in g.get("gram_stain") or []:
        emit("gram stain", e["value"], "TRUE", e["frequency"])
    doi = src.get("shape", newest)
    for e in g.get("shape") or []:
        s = SHAPE.get(e["value"])
        if s: emit("shape", s, "TRUE", e["frequency"])
        else: dropped[f"shape {e['value']}"] += 1
    doi = src.get("motility", newest)
    for e in g.get("motility") or []:
        emit("motility", "motility", "TRUE" if e["value"] else "FALSE", e["frequency"])
    doi = src.get("spore_formation", newest)
    for e in g.get("spore_formation") or []:
        emit("spore formation", "spore formation", "TRUE" if e["value"] else "FALSE", e["frequency"])
    doi = src.get("growth_temperature_c", newest)
    t = g.get("growth_temperature_c") or {}
    if t.get("optimum") is not None:
        rows["growth temperature"].append({**base, "Attribute": "growth temperature", "Attribute_value": t["optimum"], "Frequency": "always", "Unit": UNIT["growth temperature"], "Note": doi + "; optimum"})
    elif t.get("min") is not None or t.get("max") is not None:
        dropped["growth temperature range without optimum"] += 1
    doi = src.get("gc_content_mol_percent", newest)
    gc = g.get("gc_content_mol_percent") or {}
    lo, hi = gc.get("min"), gc.get("max")
    if lo is None: lo = hi
    if hi is None: hi = lo
    if lo is not None and hi is not None and lo > hi:
        lo, hi = hi, lo
    if lo is not None and 20 <= lo <= 80 and 20 <= hi <= 80:
        lo, hi = round(lo, 1), round(hi, 1)
        rows["gc content"].append({**base, "Attribute": "gc content", "Attribute_value": (f"{lo}-{hi}" if lo != hi else f"{lo}"), "Frequency": "always", "Unit": "mol%", "Note": doi})
    else:
        if lo is not None: dropped["gc content outside 20-80"] += 1
    for phys, key in (("length", "cell_length_um"), ("width", "cell_width_um")):
        doi = src.get(key, newest)
        v = g.get(key) or {}
        lo, hi = v.get("min"), v.get("max")
        if lo is None: lo = hi
        if hi is None: hi = lo
        if lo is None: continue
        if lo > hi: lo, hi = hi, lo
        if lo <= 0:
            dropped[f"{phys} not positive"] += 1; continue
        lo, hi = round(lo, 2), round(hi, 2)
        rows[phys].append({**base, "Attribute": phys, "Attribute_value": (f"{lo}-{hi}" if lo != hi else f"{lo}"), "Frequency": "always", "Unit": "um", "Note": doi})
    doi = src.get("arrangement", newest)
    for e in arr_rec.get("arrangement") or []:
        if e["value"] in ARRANGEMENT: emit("arrangement", e["value"], "TRUE", e["frequency"])
        else: dropped[f"arrangement {e['value']}"] += 1
    doi = src.get("spore_shape", newest)
    for e in arr_rec.get("spore_shape") or []:
        if e["value"] in SPORE_SHAPE: emit("spore shape", e["value"], "TRUE", e["frequency"])
        else: dropped[f"spore shape {e['value']}"] += 1
# overrides: hand corrections that survive re-extraction
ov_path = os.path.join(os.path.dirname(__file__), "overrides.tsv")
if os.path.exists(ov_path):
    for o in csv.DictReader(open(ov_path), delimiter="\t"):
        if o["action"] == "drop":
            rows[o["physiology"]] = [r for r in rows[o["physiology"]] if not (r["NCBI_ID"] == o["NCBI_ID"] and r["Attribute"] == o["Attribute"])]
        elif o["action"] == "add":
            rows[o["physiology"]].append({"NCBI_ID": o["NCBI_ID"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": o["Taxon_name"], "Attribute": o["Attribute"], "Attribute_value": o["Attribute_value"], "Attribute_source": SOURCE, "Frequency": o["Frequency"], "Note": "override: " + o["reason"]})
for phys, rr in rows.items():
    cols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency"] + (["Unit"] if phys in UNIT else []) + ["Note"]
    rr.sort(key=lambda r: (int(r["NCBI_ID"]), r["Attribute"]))
    with open(os.path.join(out_dir, phys.replace(" ", "_") + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rr)
    print(f"{phys}: {len(rr)} rows, {len({r['NCBI_ID'] for r in rr})} genera")
print("dropped:", dict(dropped))
print("filled from an older chapter:", dict(filled_from_older))
