#!/usr/bin/env python3
"""Build bugphyzz curation tables from The Microbe Directory (Shaaban et al. 2018).

Usage: sh python.sh build.py <microbe-directory.csv> <names.dmp> <out_dir>

Input is data/microbe-directory.csv from github.com/microbe-directory/microbe-directory
(MIT licence) at the commit pinned in run.sh: the 2018 release the bugphyzz
Google Sheets were copied from. The website microbe.directory no longer exists.

The table has names but no taxids. Bacteria and Archaea rows are matched to
NCBI taxids with taxoniq and the pinned ncbi-taxon-db release (current
scientific names), falling back to names.dmp from the NCBI new_taxdump
archive of the same date, because taxoniq has no synonym lookup. The
fallback compares names with case and punctuation ignored (The Microbe
Directory wrote "Frankia sp BMG5 12" for Frankia sp. BMG5.12) against
scientific names, equivalent names and synonyms, the last also without
their authority ("Prevotella copri Hayashi et al. 2007" -> "Prevotella
copri"). A name must point to exactly one current taxon; names that match
none or several are listed in unmatched.tsv and not imported.

The Microbe Directory gives no frequencies, so every row has Frequency
"unknown", as the sheet rows did; how much weight that carries is decided
at conflict resolution (waldronlab/bugphyzzExports#66). The value as
published is kept in Note.

Corrections go in overrides.tsv (action drop or add), never in the outputs.
"""
import sys, os, csv, re, collections, importlib.metadata
import taxoniq

DB_VERSION = "2026.11.1"   # https://github.com/lwaldron/taxoniq/releases/tag/ncbi-taxon-db-2026.11.1
got = importlib.metadata.version("ncbi-taxon-db")
if got != DB_VERSION:
    sys.exit(f"ncbi-taxon-db {got} installed; this build is pinned to {DB_VERSION} (see python.sh)")

md_csv, names_dmp, out_dir = sys.argv[1:4]
SOURCE = "The Microbe Directory"
stats = collections.Counter()

def norm(x):
    x = re.sub(r"[\s._\-/:,()'\[\]]+", " ", x.lower()).strip()
    return re.sub(r"\bsp\b", "sp", x)


## Fallback index from names.dmp: normalised name -> set of taxids
syn = collections.defaultdict(set)
with open(names_dmp, encoding="utf-8") as fh:
    for line in fh:
        f = line.rstrip("\t|\n").split("\t|\t")
        if len(f) < 4 or f[3] not in ("scientific name", "synonym", "equivalent name", "authority"):
            continue
        tid, name = int(f[0]), f[1]
        if f[3] != "authority":
            syn[norm(name)].add(tid)
        if f[3] in ("synonym", "authority"):
            ## drop the authority: keep "Genus species" or "Genus species subsp. x"
            m = re.match(r"^((?:Candidatus )?[A-Z][a-z]+ [a-z][a-z-]+(?: subsp\. [a-z][a-z-]+)?) (?:\(|[A-Z]|et al|ex )", name)
            if m and re.search(r"\b(1[789]|20)\d\d\b", name):   # an authority has a year
                syn[norm(m.group(1))].add(tid)


def current(tid):
    try:
        t = taxoniq.Taxon(tid)
        return t.tax_id, t.scientific_name
    except Exception:
        return None


def resolve(name):
    name = " ".join(name.split())
    for c in (name, re.sub(r" sp ", " sp. ", name)):
        try:
            t = taxoniq.Taxon(scientific_name=c)
            stats["matched by current name"] += 1
            return (t.tax_id, t.scientific_name), None
        except Exception:
            pass
    ids = syn.get(norm(name))
    if ids:
        cur = {current(i) for i in ids} - {None}
        if len(cur) > 1:
            ## a species and its own strains: take the one that contains the others
            top = [c for c in cur if all(c[0] == o[0] or c[0] in {a.tax_id for a in taxoniq.Taxon(o[0]).lineage} for o in cur)]
            if len(top) == 1:
                cur = set(top)
        if len(cur) == 1:
            stats["matched by synonym or punctuation-insensitive name"] += 1
            return cur.pop(), None
        if len(cur) > 1:
            return None, "ambiguous: " + ";".join(str(x[0]) for x in sorted(cur))
    return None, "no NCBI match"


BINARY = {   # column -> (physiology, Attribute)
    "animal_pathogen": ("animal pathogen", "animal pathogen"),
    "plant_pathogen": ("plant pathogenicity", "plant pathogenicity"),
    "antimicrobial_susceptibility": ("antimicrobial sensitivity", "antimicrobial sensitivity"),
    "biofilm_forming": ("biofilm formation", "biofilm formation"),
    "extreme_environment": ("extreme environment", "extreme environment"),
    "spore_forming": ("spore formation", "spore formation"),
}
GRAM = {"0": "gram stain negative", "1": "gram stain positive", "2": "gram stain variable"}
COGEM = {"1": "COGEM pathogenicity rating 1", "2": "COGEM pathogenicity rating 2", "3": "COGEM pathogenicity rating 3"}
NUMERIC = {"optimal_temperature": ("growth temperature", "growth temperature", "Celsius"),
           "optimal_ph": ("optimal ph", "optimal ph", "pH")}

rows = collections.defaultdict(list)
unmatched = []
seen = set()
for r in csv.DictReader(open(md_csv, newline="", encoding="utf-8")):
    if r["kingdom"] not in ("Bacteria", "Archaea"):
        continue
    stats["Bacteria and Archaea rows"] += 1
    hit, why = resolve(r["species"])
    if hit is None:
        unmatched.append({"microbe_id": r["microbe_id"], "species": r["species"], "reason": why})
        continue
    tid, cur_name = hit
    ## Two names can now be one taxon (S. oligofermentans is S. cristatus);
    ## both keep their rows, so a disagreement between them stays visible.
    if tid in seen:
        stats["names sharing a taxid with an earlier name"] += 1
    seen.add(tid)
    base = {"NCBI_ID": str(tid), "Genome_ID": "", "Accession_ID": "", "Taxon_name": cur_name,
            "Attribute_source": SOURCE, "Frequency": "unknown"}
    note = lambda col: f"{r['species']}; {col} = {r[col]}"
    for col, (phys, attr) in BINARY.items():
        if r[col] in ("0", "1"):
            rows[phys].append({**base, "Attribute": attr, "Attribute_value": "TRUE" if r[col] == "1" else "FALSE", "Note": note(col)})
    ## microbiome_location is 1 when the microbe has been found in the human
    ## microbiome; 0 is absence of evidence, not evidence of no host.
    if r["microbiome_location"] == "1":
        rows["host-associated"].append({**base, "Attribute": "host-associated", "Attribute_value": "TRUE", "Note": note("microbiome_location")})
    if r["gram_stain"] in GRAM:
        rows["gram stain"].append({**base, "Attribute": GRAM[r["gram_stain"]], "Attribute_value": "TRUE", "Note": note("gram_stain")})
    if r["pathogenicity"] in COGEM:
        rows["COGEM pathogenicity rating"].append({**base, "Attribute": COGEM[r["pathogenicity"]], "Attribute_value": "TRUE", "Note": note("pathogenicity")})
    for col, (phys, attr, unit) in NUMERIC.items():
        if r[col]:
            rows[phys].append({**base, "Attribute": attr, "Attribute_value": r[col], "Unit": unit, "Note": note(col)})

## overrides: corrections that survive a re-run
ov_path = os.environ.get("MD_OVERRIDES") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "overrides.tsv")   # MD_OVERRIDES: tests only
for o in csv.DictReader(open(ov_path, encoding="utf-8"), delimiter="\t"):
    phys = o["physiology"]
    if o["action"] == "drop":
        before = len(rows[phys])
        rows[phys] = [x for x in rows[phys] if not (x["NCBI_ID"] == o["NCBI_ID"] and x["Attribute"] == o["Attribute"])]
        if len(rows[phys]) == before:
            sys.exit(f"override matched nothing: {o}")
    elif o["action"] == "add":
        new = {"NCBI_ID": o["NCBI_ID"], "Genome_ID": "", "Accession_ID": "", "Taxon_name": o["Taxon_name"],
               "Attribute": o["Attribute"], "Attribute_value": o["Attribute_value"], "Attribute_source": SOURCE,
               "Frequency": o["Frequency"], "Note": "override: " + o["reason"]}
        if phys in ("growth temperature", "optimal ph"):
            new["Unit"] = NUMERIC["optimal_temperature" if phys == "growth temperature" else "optimal_ph"][2]
        rows[phys].append(new)
    else:
        sys.exit(f"unknown override action: {o['action']}")
    stats["overrides applied"] += 1

for phys in rows:   # identical rows from two names for one taxon
    uniq = {}
    for x in rows[phys]:
        uniq.setdefault((x["NCBI_ID"], x["Attribute"], x["Attribute_value"]), x)
    rows[phys] = list(uniq.values())

os.makedirs(out_dir, exist_ok=True)
for phys, rr in sorted(rows.items()):
    cols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency"] + (["Unit"] if phys in ("growth temperature", "optimal ph") else []) + ["Note"]
    rr.sort(key=lambda x: (int(x["NCBI_ID"]), x["Attribute"]))
    with open(os.path.join(out_dir, phys.replace(" ", "_") + ".tsv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rr)
    print(f"{phys}: {len(rr)} rows")
with open(os.path.join(out_dir, "unmatched.tsv"), "w", newline="", encoding="utf-8") as fh:
    w = csv.DictWriter(fh, fieldnames=["microbe_id", "species", "reason"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(unmatched)
for k, v in stats.items():
    print(f"{k}: {v}")
print(f"unmatched: {len(unmatched)}")
