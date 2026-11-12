#!/usr/bin/env python3 -I
"""Turn the Gemini extraction (one JSON list per batch) into bugphyzz curation-sheet rows.

Usage: build_traits.py <extraction_dir> <records.jsonl> <out_dir> [<extremophily_dir> [<revisions.tsv>]]
<records.jsonl> is collect.py's output: the taxa, their rank, article title and text.
<extremophily_dir> is the second pass (system_prompt_extremophily.txt): halophily,
temperature class and ph class.
<revisions.tsv> is revisions.py's output; with it the Note links the revision the
text came from, without it the live article.
Writes one TSV per physiology in the Google Sheet column layout
(NCBI_ID, Genome_ID, Accession_ID, Taxon_name, Attribute, Attribute_value,
Attribute_source, Frequency, [Unit,] Note), with overrides.tsv applied.
"""
import json, sys, os, glob, csv, collections, re, urllib.parse
ext_dir, records_file, out_dir = sys.argv[1:4]
extremo_dir = sys.argv[4] if len(sys.argv) > 4 else None
revisions = {}
if len(sys.argv) > 5:
    for r in csv.DictReader(open(sys.argv[5]), delimiter="\t"):
        if r["revid"]: revisions[r["taxid"]] = r["revid"]
os.makedirs(out_dir, exist_ok=True)
SOURCE = "Wikipedia_2026"
records = {}
for line in open(records_file):
    r = json.loads(line); records[r["taxid"]] = r

def article_note(rec):
    """Link to the article (the revision the text came from, if known) and the Wikidata item."""
    title = urllib.parse.quote(rec["en_wiki_title"].replace(" ", "_"), safe="()_,'")
    if rec["taxid"] in revisions:
        url = f"https://en.wikipedia.org/w/index.php?title={title}&oldid={revisions[rec['taxid']]}"
    else:
        url = "https://en.wikipedia.org/wiki/" + title
    return url + ("; wikidata " + rec["wikidata_id"] if rec["wikidata_id"] else "")
# bugphyzz exports genus, species and strain (ADR 0004); NCBI ranks below species count as strain
RANKS = {"genus", "species", "subspecies", "strain", "biotype", "serotype", "serogroup", "varietas", "forma_specialis", "isolate", "pathogroup"}
# extraction vocabulary -> bugphyzz attributes.tsv vocabulary (None = drop)
SHAPE = {"bacillus": "bacillus", "coccus": "coccus", "coccobacillus": "coccobacillus", "spirillum": "spirillum",
         "filamentous": "filamentous", "elliptic": "elliptic", "pleomorphic": "pleiomorphic cell", "vibrio": "vibrio", "other": None}
AERO = {"aerobic", "anaerobic", "facultatively anaerobic", "microaerophilic"}
GRAM = {"gram stain positive", "gram stain negative", "gram stain variable"}
FREQ = {"always", "usually", "sometimes"}
ARRANGEMENT = {"single", "paired cells", "cell chain", "cell cluster", "tetrads", "palisade", "rosettes", "v-formation", "branched", "filamentous"}
SPORE_SHAPE = {"elliptic", "coccus", "bacillus"}
HABITAT = {"soil", "rhizosphere", "root", "root nodule", "phylloplane", "plant matter", "fresh water", "sea water", "sediment", "mud",
           "groundwater", "hot spring", "hydrothermal vent", "saline evaporation pond", "waste water", "sludge", "food product", "dairy",
           "rumen", "gastrointestinal tract", "feces", "oral", "respiratory system", "nasopharyngeal cavity", "vagina", "human urine",
           "blood", "insect", "fish", "mammal", "air", "mine", "wetland", "ocean trench", "oil field production water",
           "skin", "intertidal zone"}
# habitats the first run returned as "other: ..." that are terms of the vocabulary under another name
OTHER_HABITAT = {"skin": "skin", "skin flora": "skin", "human skin": "skin",
                 "tidal flat": "intertidal zone", "tidal flats": "intertidal zone", "intertidal zone": "intertidal zone"}
HALOPHILY = {"halophilic", "halotolerant", "non-halophilic"}
TEMPERATURE_CLASS = {"psychrophile", "psychrotolerant", "mesophile", "thermotolerant", "thermophile", "hyperthermophile"}
PH_CLASS = {"acidophile", "acidotolerant", "neutrophile", "alkalitolerant", "alkaliphile"}
UNIT = {"growth temperature": "Celcius", "gc content": "mol%", "width": "um", "length": "um", "optimal ph": "pH"}
# A row is kept only if the introduction has a word for the trait: the model is told to use the
# text alone, and this catches values it supplied from what it knows about the taxon instead.
SUPPORT = {
    "aerophilicity": r"aerob|aeroph|aerotoler|oxygen|anoxi",
    "gram stain": r"gram",
    "shape": r"rod|bacill|cocc|spher|curved|vibrio|comma|spiral|helic|spiroch|filament|pleomorph|ovoid|oval|ellips|shape|hypha|mycel|trichome|corkscrew",
    "motility": r"motil|flagell|swim|glid|swarm|twitch",
    "spore formation": r"spor",
    "arrangement": r"pair|chain|cluster|tetrad|single|singly|diplo|strepto|staphylo|packet|sarcin|aggregat|clump|palisade|rosette|branch|hypha|mycel|filament|trichome|grape",
    "spore shape": r"spor",
    "animal pathogen": r"pathogen|diseas|infect|caus|virulen|toxi|harmless|sepsis|septic|itis\b|osis\b|abscess|zoono|parasit|plague|fever",
    "plant pathogenicity": r"pathogen|diseas|infect|caus|blight|rot\b|wilt|canker|gall|spot|scab|phytoplasma|virulen",
    "halophily class": r"halo|salt",
    "temperature class": r"thermo|psychro|cryo|mesoph|cold",
    "ph class": r"acid|alkal|neutroph",
    "host-associated": r"host|commensal|symbio|parasit|pathogen|flora|microbio|gut|intestin|skin|nodul|endophyt|free-living|infect|associat|colon|inhabit|oral|mouth|rumen|vagin|tract",
}
rows = collections.defaultdict(list); dropped = collections.Counter(); other_habitat = collections.Counter(); taxa = {}
for f in sorted(glob.glob(os.path.join(ext_dir, "batch_*.json"))):
    for g in json.load(open(f)):
        taxa[str(g["taxid"])] = g
missing = set(records) - set(taxa)
if missing: print(f"WARNING: {len(missing)} taxa have no extraction", file=sys.stderr)
for taxid in sorted(taxa, key=int):
    g = taxa[taxid]; rec = records.get(taxid)
    if rec is None:
        dropped["taxid not in records"] += 1; continue
    if rec["rank"] not in RANKS:
        dropped["rank above genus"] += 1; continue
    if not g.get("about_taxon"):
        dropped["text not about the taxon"] += 1; continue
    low = rec["text"].lower()
    note = article_note(rec)
    base = {"NCBI_ID": taxid, "Genome_ID": "", "Accession_ID": "", "Taxon_name": rec["name"], "Attribute_source": SOURCE, "Note": note}
    def emit(phys, attr, val, freq):
        if freq not in FREQ: dropped[f"frequency {freq}"] += 1; return
        if phys in SUPPORT and not re.search(SUPPORT[phys], low):
            dropped[f"{phys}: no word for it in the text"] += 1; return
        rows[phys].append({**base, "Attribute": attr, "Attribute_value": val, "Frequency": freq})
    def entries(key):
        return [e for e in (g.get(key) or []) if isinstance(e, dict) and "value" in e and "frequency" in e]
    def tf(v): return "TRUE" if v is True else "FALSE" if v is False else None
    for e in entries("aerophilicity"):
        emit("aerophilicity", e["value"], "TRUE", e["frequency"]) if e["value"] in AERO else dropped.update([f"aerophilicity {e['value']}"])
    for e in entries("gram_stain"):
        emit("gram stain", e["value"], "TRUE", e["frequency"]) if e["value"] in GRAM else dropped.update([f"gram stain {e['value']}"])
    for e in entries("shape"):
        s = SHAPE.get(e["value"])
        if s: emit("shape", s, "TRUE", e["frequency"])
        else: dropped[f"shape {e['value']}"] += 1
    for key, phys in (("motility", "motility"), ("spore_formation", "spore formation"), ("animal_pathogen", "animal pathogen"),
                      ("plant_pathogen", "plant pathogenicity"), ("host_associated", "host-associated")):
        for e in entries(key):
            v = tf(e["value"])
            if v: emit(phys, phys, v, e["frequency"])
            else: dropped[f"{phys} not true or false"] += 1
    for e in entries("arrangement"):
        if e["value"] in ARRANGEMENT: emit("arrangement", e["value"], "TRUE", e["frequency"])
        else: dropped[f"arrangement {e['value']}"] += 1
    for e in entries("spore_shape"):
        if e["value"] in SPORE_SHAPE: emit("spore shape", e["value"], "TRUE", e["frequency"])
        else: dropped[f"spore shape {e['value']}"] += 1
    for e in entries("habitat"):
        v = str(e["value"])
        if v in HABITAT: emit("habitat", v, "TRUE", e["frequency"])
        elif v.startswith("other:"):
            o = v[6:].strip().lower()
            if o in OTHER_HABITAT: emit("habitat", OTHER_HABITAT[o], "TRUE", e["frequency"])
            else: other_habitat[o] += 1; dropped["habitat outside the vocabulary"] += 1
        else: dropped[f"habitat {v}"] += 1
    t = g.get("growth_temperature_c") or {}
    if isinstance(t.get("optimum"), (int, float)) and -5 <= t["optimum"] <= 122:
        rows["growth temperature"].append({**base, "Attribute": "growth temperature", "Attribute_value": t["optimum"], "Frequency": "always", "Unit": UNIT["growth temperature"], "Note": note + "; optimum"})
    elif t.get("min") is not None or t.get("max") is not None or t.get("optimum") is not None:
        dropped["growth temperature without a usable optimum"] += 1
    ph = g.get("optimal_ph")
    if isinstance(ph, (int, float)) and 0 < ph < 14:
        rows["optimal ph"].append({**base, "Attribute": "optimal ph", "Attribute_value": ph, "Frequency": "always", "Unit": UNIT["optimal ph"]})
    elif ph is not None: dropped["optimal ph not a number in 0-14"] += 1
    def rng(v):
        lo, hi = (v or {}).get("min"), (v or {}).get("max")
        if not isinstance(lo, (int, float)): lo = None
        if not isinstance(hi, (int, float)): hi = None
        if lo is None: lo = hi
        if hi is None: hi = lo
        if lo is not None and lo > hi: lo, hi = hi, lo
        return lo, hi
    lo, hi = rng(g.get("gc_content_mol_percent"))
    if lo is not None and 20 <= lo <= 80 and 20 <= hi <= 80:
        lo, hi = round(lo, 1), round(hi, 1)
        rows["gc content"].append({**base, "Attribute": "gc content", "Attribute_value": (f"{lo}-{hi}" if lo != hi else f"{lo}"), "Frequency": "always", "Unit": "mol%"})
    elif lo is not None: dropped["gc content outside 20-80"] += 1
    for phys, key in (("length", "cell_length_um"), ("width", "cell_width_um")):
        lo, hi = rng(g.get(key))
        if lo is None: continue
        if lo <= 0:
            dropped[f"{phys} not positive"] += 1; continue
        lo, hi = round(lo, 2), round(hi, 2)
        rows[phys].append({**base, "Attribute": phys, "Attribute_value": (f"{lo}-{hi}" if lo != hi else f"{lo}"), "Frequency": "always", "Unit": "um"})
# second pass: halophily, temperature class, ph class
if extremo_dir:
    for f in sorted(glob.glob(os.path.join(extremo_dir, "batch_*.json"))):
        for g in json.load(open(f)):
            taxid = str(g["taxid"]); rec = records.get(taxid)
            if rec is None or rec["rank"] not in RANKS or not g.get("about_taxon"): continue
            low = rec["text"].lower()
            note = article_note(rec)
            for key, phys, vocab in (("halophily", "halophily class", HALOPHILY), ("temperature_class", "temperature class", TEMPERATURE_CLASS), ("ph_class", "ph class", PH_CLASS)):
                for e in g.get(key) or []:
                    if not isinstance(e, dict) or e.get("value") not in vocab: dropped[f"{phys} {e.get('value') if isinstance(e, dict) else e}"] += 1; continue
                    if e.get("frequency") not in FREQ: dropped[f"frequency {e.get('frequency')}"] += 1; continue
                    if not re.search(SUPPORT[phys], low): dropped[f"{phys}: no word for it in the text"] += 1; continue
                    n = note + (f"; {e['degree']}" if e.get("degree") in ("slight", "moderate", "extreme") else "")
                    rows[phys].append({"NCBI_ID": taxid, "Genome_ID": "", "Accession_ID": "", "Taxon_name": rec["name"], "Attribute": e["value"],
                                       "Attribute_value": "TRUE", "Attribute_source": SOURCE, "Frequency": e["frequency"], "Note": n})
# a taxon given both TRUE and FALSE, or the same value twice, for one binary trait
for phys in ("motility", "spore formation", "animal pathogen", "plant pathogenicity", "host-associated"):
    seen = set(); out = []
    for r in rows[phys]:
        k = (r["NCBI_ID"], r["Attribute_value"])
        if k in seen: dropped[f"{phys} duplicate"] += 1; continue
        seen.add(k); out.append(r)
    rows[phys] = out
for phys in list(rows):
    seen = set(); out = []
    for r in rows[phys]:
        k = (r["NCBI_ID"], r["Attribute"], str(r["Attribute_value"]))
        if k in seen: dropped[f"{phys} duplicate"] += 1; continue
        seen.add(k); out.append(r)
    rows[phys] = out
# overrides: hand corrections that survive re-extraction
ov_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "overrides.tsv")
if os.path.exists(ov_path):
    for o in csv.DictReader(open(ov_path), delimiter="\t"):
        if o["action"] == "drop":
            rows[o["physiology"]] = [r for r in rows[o["physiology"]] if not (r["NCBI_ID"] == o["NCBI_ID"] and r["Attribute"] == o["Attribute"])]
        elif o["action"] == "add":
            rows[o["physiology"]].append({"NCBI_ID": o["NCBI_ID"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": o["Taxon_name"], "Attribute": o["Attribute"], "Attribute_value": o["Attribute_value"], "Attribute_source": SOURCE, "Frequency": o["Frequency"], "Note": "override: " + o["reason"]})
for phys, rr in sorted(rows.items()):
    if not rr: continue
    cols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency"] + (["Unit"] if phys in UNIT else []) + ["Note"]
    rr.sort(key=lambda r: (int(r["NCBI_ID"]), r["Attribute"], str(r["Attribute_value"])))
    with open(os.path.join(out_dir, phys.replace(" ", "_") + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rr)
    rk = collections.Counter("genus" if records[r["NCBI_ID"]]["rank"] == "genus" else "species" if records[r["NCBI_ID"]]["rank"] == "species" else "strain" for r in rr if r["NCBI_ID"] in records)
    print(f"{phys}: {len(rr)} rows, {len({r['NCBI_ID'] for r in rr})} taxa ({rk['genus']} genus rows, {rk['species']} species, {rk['strain']} below species)")
print("dropped:", dict(sorted(dropped.items(), key=lambda kv: -kv[1])))
print("habitats outside the vocabulary (20 most frequent):", other_habitat.most_common(20))
