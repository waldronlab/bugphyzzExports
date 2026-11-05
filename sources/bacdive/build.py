#!/usr/bin/env python3 -I
"""Build the BacDive source files from raw API records.

Usage: build_bacdive.py <raw_dir> <taxdump_dir> <out_dir>

Writes:
  bacdive.csv.gz   one row per strain in the layout of the former BacDive
                   Google Sheet, read by bugphyzz:::.getBacDive()
  gc_content.tsv, length.tsv, width.tsv   curation-sheet rows for the three
                   attributes the sheet never carried (read via manifest.tsv)
"""
import sys, os, json, gzip, glob, csv, re, collections
raw_dir, taxdump, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)

# NCBI taxdump: rank per taxid, scientific name per taxid, species name -> taxid
rank, sci = {}, {}
with open(os.path.join(taxdump, "nodes.dmp")) as f:
    for line in f:
        p = line.split("|"); rank[p[0].strip()] = p[2].strip()
name2id = {}
with open(os.path.join(taxdump, "names.dmp")) as f:
    for line in f:
        p = [x.strip() for x in line.split("|")]
        if p[3] == "scientific name":
            sci[p[0]] = p[1]
            if rank.get(p[0]) == "species": name2id.setdefault(p[1], p[0])
merged = {}
with open(os.path.join(taxdump, "merged.dmp")) as f:
    for line in f:
        p = [x.strip() for x in line.split("|")]; merged[p[0]] = p[1]

def aslist(x):
    return x if isinstance(x, list) else ([] if x in (None, {}) else [x])
def first(x, key):
    for e in aslist(x):
        if isinstance(e, dict) and e.get(key) not in (None, ""): return str(e[key])
    return ""
def joined(x, key, sep=";"):
    vals = []
    for e in aslist(x):
        if isinstance(e, dict) and e.get(key) not in (None, "") and str(e[key]) not in vals: vals.append(str(e[key]))
    return sep.join(vals)

COLS = ["BacDive_ID", "Taxon_name", "NCBI_ID", "Rank", "Parent_name", "Parent_NCBI_ID", "Parent_rank", "Seq16S_NCBI_ID", "Genome_ID", "Type_strain",
        "Gram_stain", "Cell_shape", "Motility", "Hemolysis", "Colony_color", "Incubation_period", "Cultivation_medium_used",
        "Culture_temperature_growth", "Culture_temperature", "Culture_temperature_range", "Oxygen_tolerance", "Spore_formation", "Halophily",
        "Metabolite_utiilization", "Metabolite_production", "Sample_type", "Geographic_location", "Country",
        "Biosafety_level", "Biosafety_level_comment", "Pathogenicity_human", "Pathogenicity_animal"]
rows, gc_rows, len_rows, wid_rows = [], [], [], []
stats = collections.Counter()
for f in sorted(glob.glob(os.path.join(raw_dir, "*.json.gz"))):
    with gzip.open(f, "rt") as fh:
        recs = json.load(fh)
    for sid, d in recs.items():
        gen = d.get("General", {}); tax = d.get("Name and taxonomic classification", {})
        morph = d.get("Morphology", {}); phys = d.get("Physiology and metabolism", {})
        cult = d.get("Culture and growth conditions", {}); iso = d.get("Isolation, sampling and environmental information", {})
        safe = d.get("Safety information", {}); seq = d.get("Sequence information", {})
        species = tax.get("species") or ""; desig = tax.get("strain designation") or ""
        taxon_name = (species + " " + desig).strip() if desig else species
        # NCBI ids: BacDive gives the best match and its level; strains below species get their own id
        ncbi = aslist(gen.get("NCBI tax id"))
        strain_id = species_id = ""
        for e in ncbi:
            if not isinstance(e, dict): continue
            tid = str(e.get("NCBI tax id", "")); tid = merged.get(tid, tid)
            lvl = e.get("Matching level", "")
            if lvl == "strain" and not strain_id: strain_id = tid
            elif lvl == "species" and not species_id: species_id = tid
        if not species_id and species: species_id = name2id.get(species, "")
        if not species_id and strain_id and rank.get(strain_id) == "species": species_id, strain_id = strain_id, ""
        # the former sheet also took a strain-level taxid from the genome assembly record
        if not strain_id:
            for e in aslist(d.get("Sequence information", {}).get("Genome sequences")):
                tid = str(e.get("NCBI tax ID", "")) if isinstance(e, dict) else ""
                tid = merged.get(tid, tid)
                if tid and tid != species_id and rank.get(tid) not in (None, "species", "genus", "family"):
                    strain_id = tid; stats["strain taxid from genome record"] += 1; break
        stats["strains"] += 1; stats["with strain taxid"] += bool(strain_id); stats["with species taxid"] += bool(species_id)
        cm = aslist(morph.get("cell morphology"))
        temps = [e for e in aslist(cult.get("culture temp")) if isinstance(e, dict)]
        opt = [e for e in temps if e.get("type") == "optimum" and e.get("growth") in (None, "positive")]
        grow = [e for e in temps if e.get("type") == "growth" and e.get("growth") == "positive"]
        t_entry = (opt or grow or temps or [{}])[0]
        halo = []
        for e in aslist(phys.get("halophily")):
            if isinstance(e, dict) and e.get("salt"):
                halo.append(f"{e['salt']} {e.get('growth','')} {e.get('tested relation','')} {e.get('concentration','')}".replace("  ", " ").strip())
        mu = []
        for e in aslist(phys.get("metabolite utilization")):
            if isinstance(e, dict) and e.get("metabolite") and e.get("utilization activity"):
                mu.append(f"{e['metabolite']} {e['utilization activity']} {e.get('kind of utilization tested','')}".strip())
        mp = []
        for e in aslist(phys.get("metabolite production")):
            if isinstance(e, dict) and e.get("metabolite") and e.get("production"):
                mp.append(f"{e['metabolite']} {e['production']}")
        risk = aslist(safe.get("risk assessment"))
        isol = aslist(iso.get("isolation"))
        genomes = aslist(seq.get("Genome sequences"))
        row = {
            "BacDive_ID": sid, "Taxon_name": taxon_name, "NCBI_ID": strain_id, "Rank": "strain",
            "Parent_name": species, "Parent_NCBI_ID": species_id, "Parent_rank": "species" if species_id else "",
            "Seq16S_NCBI_ID": first(seq.get("16S sequences"), "accession"), "Genome_ID": first(genomes, "accession"),
            "Type_strain": tax.get("type strain", ""),
            "Gram_stain": joined(cm, "gram stain"), "Cell_shape": joined(cm, "cell shape"), "Motility": joined(cm, "motility"),
            "Hemolysis": joined(morph.get("colony morphology"), "type of hemolysis"),
            "Colony_color": joined(morph.get("colony morphology"), "colony color"), "Incubation_period": joined(morph.get("colony morphology"), "incubation period"),
            "Cultivation_medium_used": ";".join(e["name"] for e in aslist(cult.get("culture medium")) if isinstance(e, dict) and e.get("name") and e.get("growth") == "yes"),
            "Culture_temperature_growth": t_entry.get("growth", ""), "Culture_temperature": t_entry.get("temperature", ""),
            "Culture_temperature_range": ";".join(f"{e.get('type','')} {e.get('growth','')} {e.get('temperature','')}".strip() for e in temps),
            "Oxygen_tolerance": joined(phys.get("oxygen tolerance"), "oxygen tolerance"),
            "Spore_formation": joined(phys.get("spore formation"), "spore formation"),
            "Halophily": ";".join(halo), "Metabolite_utiilization": ";".join(mu), "Metabolite_production": ";".join(mp),
            "Sample_type": joined(isol, "sample type"), "Geographic_location": joined(isol, "geographic location"), "Country": joined(isol, "country"),
            "Biosafety_level": first(risk, "biosafety level"), "Biosafety_level_comment": first(risk, "biosafety level comment"),
            "Pathogenicity_human": first(risk, "pathogenicity human"), "Pathogenicity_animal": first(risk, "pathogenicity animal"),
        }
        rows.append(row)
        # narrow files for attributes the sheet never had; rows need an NCBI id to be useful
        tid, tname = (strain_id, taxon_name) if strain_id else (species_id, species)
        if not tid: continue
        base = {"NCBI_ID": tid, "Genome_ID": "", "Accession_ID": "", "Taxon_name": tname, "Attribute_source": "BacDive", "Frequency": "always", "Note": f"BacDive-ID {sid}"}
        gc = first(seq.get("GC content"), "GC-content")
        m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?\s*", gc or "")
        if m and 20 <= float(m.group(1)) <= 80:
            gc_rows.append({**base, "Attribute": "gc content", "Attribute_value": m.group(1) + (f"-{m.group(2)}" if m.group(2) else ""), "Unit": "mol%"})
        for key, store, attr in (("cell length", len_rows, "length"), ("cell width", wid_rows, "width")):
            v = joined(cm, key, sep="|").split("|")[0]
            m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?\s*(?:µm|um|μm)?\s*", v or "")
            if m: store.append({**base, "Attribute": attr, "Attribute_value": m.group(1) + (f"-{m.group(2)}" if m.group(2) else ""), "Unit": "um"})
rows.sort(key=lambda r: int(r["BacDive_ID"]))
with gzip.open(os.path.join(out_dir, "bacdive.csv.gz"), "wt", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=COLS, lineterminator="\n"); w.writeheader(); w.writerows(rows)
ncols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency", "Unit", "Note"]
for name, store in (("gc_content", gc_rows), ("length", len_rows), ("width", wid_rows)):
    store.sort(key=lambda r: (int(r["NCBI_ID"]), r["Note"]))
    with open(os.path.join(out_dir, name + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ncols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(store)
    print(f"{name}: {len(store)} rows")
for c in COLS[10:]:
    stats["col " + c] = sum(1 for r in rows if r[c])
print(json.dumps(stats, indent=1))
