#!/usr/bin/env python3
"""Build genome size, G+C content and coding-gene tables from the RefSeq assembly summary.

Usage: python3 -I build.py <assembly_summary_refseq.txt> <out_dir>
Needs taxoniq with the ncbi-taxon-db release pinned in DB_VERSION (see python.sh).

Keeps bacterial and archaeal assemblies that RefSeq has not flagged
(excluded_from_refseq is "na", so no metagenome-derived, partial or
contaminated assemblies) and picks one per species_taxid: a reference genome
first, then an assembly from type material, then the most complete assembly
level, then the newest release. The species taxid is first resolved against
the pinned ncbi-taxon-db (dropped if absent, moved to its species if no
longer at species rank), so the choice is per current species taxid and no
taxid appears twice; the current scientific name is the Taxon_name.

Writes genome_size.tsv, gc_content.tsv and coding_genes.tsv in the
curation-sheet layout; the assembly accession, level and flags go in Note.
"""
import sys, os, csv, collections, importlib.metadata, math
import taxoniq

DB_VERSION = "2026.11.1"   # https://github.com/lwaldron/taxoniq/releases/tag/ncbi-taxon-db-2026.11.1
got = importlib.metadata.version("ncbi-taxon-db")
if got != DB_VERSION:
    sys.exit(f"ncbi-taxon-db {got} installed; this build is pinned to {DB_VERSION} (see python.sh)")

summary, out_dir = sys.argv[1:3]
os.makedirs(out_dir, exist_ok=True)
SOURCE = "NCBI_RefSeq"
LEVEL = {"Complete Genome": 0, "Chromosome": 1, "Scaffold": 2, "Contig": 3}
TYPE = ("assembly from type material", "assembly designated as neotype", "assembly designated as reftype")

stats = collections.Counter()
resolved = {}   # species_taxid in the summary -> taxoniq Taxon at species rank, or None

def species(sp):
    if sp in resolved:
        return resolved[sp]
    try:
        t = taxoniq.Taxon(int(sp))
    except KeyError:
        stats["species taxid not in ncbi-taxon-db, dropped"] += 1; t = None
    else:
        if t.rank.name != "species":
            t = next((x for x in t.ranked_lineage if x.rank.name == "species"), None)
            stats["species taxid re-ranked, moved to its species" if t else "species taxid with no species in its lineage, dropped"] += 1
    resolved[sp] = t
    return t

def positive(x):
    try:
        v = float(x)
    except ValueError:
        return False
    return math.isfinite(v) and v > 0

best = {}
with open(summary) as fh:
    first = fh.readline()
    if not first.startswith("##"):
        sys.exit(f"{summary} does not look like an NCBI assembly summary (no '##' comment line)")
    reader = csv.DictReader(fh, delimiter="\t")
    reader.fieldnames[0] = reader.fieldnames[0].lstrip("#")
    for row in reader:
        stats["assemblies"] += 1
        if row["group"] not in ("bacteria", "archaea"):
            continue
        stats["bacteria and archaea"] += 1
        if row["excluded_from_refseq"] != "na":
            stats["flagged by RefSeq, skipped"] += 1; continue
        if row["genome_rep"] != "Full":
            stats["partial genome representation, skipped"] += 1; continue
        if row["version_status"] != "latest":
            stats["not the latest version, skipped"] += 1; continue
        stats["candidate assemblies"] += 1
        t = species(row["species_taxid"])
        if t is None:
            continue
        date = row["seq_rel_date"].replace("/", "")   # YYYY/MM/DD or "na"
        key = (row["refseq_category"] != "reference genome",
               row["relation_to_type_material"] not in TYPE,
               LEVEL[row["assembly_level"]],
               -int(date) if date.isdigit() else 0,   # newest first, "na" last
               row["assembly_accession"])
        if t.tax_id not in best or key < best[t.tax_id][0]:
            best[t.tax_id] = (key, row, t)

ncols = ["NCBI_ID", "Genome_ID", "Accession_ID", "Taxon_name", "Attribute", "Attribute_value", "Attribute_source", "Frequency", "Unit", "Note"]
tables = {"genome_size": [], "gc_content": [], "coding_genes": []}
for _, row, t in best.values():
    stats["species"] += 1
    flags = [row["assembly_level"]]
    if row["refseq_category"] == "reference genome": flags.append("reference genome")
    if row["relation_to_type_material"] != "na": flags.append(row["relation_to_type_material"])
    base = {"NCBI_ID": str(t.tax_id), "Genome_ID": "", "Accession_ID": "", "Taxon_name": t.scientific_name,
            "Attribute_source": SOURCE, "Frequency": "always", "Note": f"{row['assembly_accession']}; " + "; ".join(flags)}
    for name, col, attr, unit in (("genome_size", "genome_size", "Genome size", "base pair"),
                                  ("gc_content", "gc_percent", "gc content", "mol%"),
                                  ("coding_genes", "protein_coding_gene_count", "genes", "genes")):
        if not positive(row[col]):
            stats[f"no {col}, {name} row skipped"] += 1; continue
        value = f"{float(row[col]):g}" if col == "gc_percent" else row[col]
        tables[name].append({**base, "Attribute": attr, "Attribute_value": value, "Unit": unit})

for name, rows in tables.items():
    rows.sort(key=lambda r: int(r["NCBI_ID"]))
    with open(os.path.join(out_dir, name + ".tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ncols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    stats[f"rows {name}"] = len(rows)
stats["ncbi-taxon-db"] = DB_VERSION
for k, v in stats.items():
    print(f"{k}: {v}")
