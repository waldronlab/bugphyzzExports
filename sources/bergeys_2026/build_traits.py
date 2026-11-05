#!/usr/bin/env python3 -I
"""Turn the Gemini extraction (one JSON list per batch) into bugphyzz curation-sheet rows.

Usage: build_traits.py <extraction_dir> <genus_taxids.json> <out_dir>
Writes one TSV per physiology in the Google Sheet column layout
(NCBI_ID, Genome_ID, Accession_ID, Taxon_name, Attribute, Attribute_value,
Attribute_source, Frequency, Note) plus overrides applied from overrides.tsv.
"""
import json, sys, os, glob, csv, collections
ext_dir, taxid_file, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)
tax = json.load(open(taxid_file))
dois = {d["doi"]: d for d in json.load(open(os.path.join(os.path.dirname(__file__), "dois.json")))}
SOURCE = "Bergeys_2026"
# extraction vocabulary -> bugphyzz attributes.tsv vocabulary (None = drop)
SHAPE = {"bacillus": "bacillus", "coccus": "coccus", "coccobacillus": "coccobacillus", "spirillum": "spirillum",
         "filamentous": "filamentous", "elliptic": "elliptic", "pleomorphic": "pleiomorphic cell", "vibrio": "bacillus", "other": None}
AERO = {"aerobic": "aerobic", "anaerobic": "anaerobic", "facultatively anaerobic": "facultatively anaerobic", "microaerophilic": "microaerophilic"}
FREQ = {"always", "usually", "sometimes"}
rows = collections.defaultdict(list); dropped = collections.Counter(); chapters = {}
for f in sorted(glob.glob(os.path.join(ext_dir, "batch_*.json"))):
    for g in json.load(open(f)):
        chapters[g["doi"]] = g
# newest chapter per genus
by_genus = collections.defaultdict(list)
for doi, g in chapters.items():
    by_genus[g["genus"]].append((dois.get(doi, {}).get("year") or 0, doi))
keep = {max(v)[1] for v in by_genus.values()}
for doi in sorted(keep):
    g = chapters[doi]; t = tax.get(g["genus"]) or {}
    if not t.get("taxid"):
        dropped["genus without NCBI taxid"] += 1; continue
    base = {"NCBI_ID": t["taxid"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": t["clean"], "Attribute_source": SOURCE, "Note": doi}
    def emit(phys, attr, val, freq, note=""):
        if freq not in FREQ: dropped[f"frequency {freq}"] += 1; return
        rows[phys].append({**base, "Attribute": attr, "Attribute_value": val, "Frequency": freq, "Note": (doi + ("; " + note if note else ""))})
    for e in g.get("aerophilicity") or []:
        a = AERO.get(e["value"]); emit("aerophilicity", a, "TRUE", e["frequency"]) if a else dropped.update([f"aerophilicity {e['value']}"])
    for e in g.get("gram_stain") or []:
        emit("gram stain", e["value"], "TRUE", e["frequency"])
    for e in g.get("shape") or []:
        s = SHAPE.get(e["value"])
        if s: emit("shape", s, "TRUE", e["frequency"], "vibrioid" if e["value"] == "vibrio" else "")
        else: dropped[f"shape {e['value']}"] += 1
    for e in g.get("motility") or []:
        emit("motility", "motility", "TRUE" if e["value"] else "FALSE", e["frequency"])
    for e in g.get("spore_formation") or []:
        emit("spore formation", "spore formation", "TRUE" if e["value"] else "FALSE", e["frequency"])
    t = g.get("growth_temperature_c") or {}
    if t.get("optimum") is not None:
        rows["growth temperature"].append({**base, "Attribute": "growth temperature", "Attribute_value": t["optimum"], "Frequency": "always", "Unit": "Celcius", "Note": doi + "; optimum"})
    elif t.get("min") is not None and t.get("max") is not None:
        rows["growth temperature"].append({**base, "Attribute": "growth temperature", "Attribute_value": (t["min"] + t["max"]) / 2, "Frequency": "sometimes", "Unit": "Celcius", "Note": f"{doi}; midpoint of range {t['min']}-{t['max']}"})
# overrides: hand corrections that survive re-extraction
ov_path = os.path.join(os.path.dirname(__file__), "overrides.tsv")
if os.path.exists(ov_path):
    for o in csv.DictReader(open(ov_path), delimiter="\t"):
        if o["action"] == "drop":
            rows[o["physiology"]] = [r for r in rows[o["physiology"]] if not (r["NCBI_ID"] == o["NCBI_ID"] and r["Attribute"] == o["Attribute"])]
        elif o["action"] == "add":
            rows[o["physiology"]].append({"NCBI_ID": o["NCBI_ID"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": o["Taxon_name"], "Attribute": o["Attribute"], "Attribute_value": o["Attribute_value"], "Attribute_source": SOURCE, "Frequency": o["Frequency"], "Note": "override: " + o["reason"]})
for phys, rr in rows.items():
    cols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency"] + (["Unit"] if phys == "growth temperature" else []) + ["Note"]
    rr.sort(key=lambda r: (int(r["NCBI_ID"]), r["Attribute"]))
    with open(os.path.join(out_dir, phys.replace(" ", "_") + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rr)
    print(f"{phys}: {len(rr)} rows, {len({r['NCBI_ID'] for r in rr})} genera")
print("dropped:", dict(dropped))
