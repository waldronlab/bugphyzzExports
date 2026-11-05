#!/usr/bin/env python3
"""Fixture test for build.py: sh python.sh test_build.py (needs the pinned ncbi-taxon-db, no network)."""
import os, subprocess, sys, tempfile, csv
here = os.path.dirname(os.path.abspath(__file__))
d = tempfile.mkdtemp(); out = os.path.join(d, "out"); summary = os.path.join(d, "assembly_summary_refseq.txt")
cols = ["assembly_accession", "bioproject", "biosample", "wgs_master", "refseq_category", "taxid", "species_taxid", "organism_name",
        "infraspecific_name", "isolate", "version_status", "assembly_level", "release_type", "genome_rep", "seq_rel_date", "asm_name",
        "asm_submitter", "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq", "relation_to_type_material",
        "asm_not_live_date", "assembly_type", "group", "genome_size", "genome_size_ungapped", "gc_percent", "replicon_count",
        "scaffold_count", "contig_count", "annotation_provider", "annotation_name", "annotation_date", "total_gene_count",
        "protein_coding_gene_count", "non_coding_gene_count", "pubmed_id"]
def row(acc, taxid, species, name, level, size, gc, pcg, cat="na", typ="na", date="2020/01/01", group="bacteria", excl="na", rep="Full"):
    r = dict.fromkeys(cols, "na")
    r.update(assembly_accession=acc, taxid=taxid, species_taxid=species, organism_name=name, assembly_level=level, genome_size=size,
             gc_percent=gc, protein_coding_gene_count=pcg, refseq_category=cat, relation_to_type_material=typ, seq_rel_date=date,
             group=group, excluded_from_refseq=excl, genome_rep=rep, version_status="latest")
    return r
rows = [
    # E. coli: the reference genome wins over a newer, complete, type-material assembly; strain taxid maps to species 562
    row("GCF_000005845.2", "511145", "562", "Escherichia coli str. K-12 substr. MG1655", "Complete Genome", "4641652", "50.800000", "4298", cat="reference genome"),
    row("GCF_000008865.2", "386585", "562", "Escherichia coli O157:H7 str. Sakai", "Complete Genome", "5594477", "50.500000", "5155", date="2024/01/01", typ="assembly from type material"),
    # B. subtilis: type material beats a complete non-type assembly; the "na" gene count drops only the coding_genes row
    row("GCF_000009045.1", "224308", "1423", "Bacillus subtilis subsp. subtilis str. 168", "Scaffold", "4215606", "43.500000", "na", typ="assembly from type material"),
    row("GCF_000227465.1", "1423", "1423", "Bacillus subtilis", "Complete Genome", "4100000", "43.000000", "4000", date="2024/01/01"),
    # skipped: metagenome-derived, partial, a virus, and a species taxid absent from ncbi-taxon-db (merged id 12)
    row("GCF_900000001.1", "1496", "1496", "Clostridioides difficile", "Contig", "4000000", "28.000000", "3500", excl="derived from metagenome"),
    row("GCF_900000002.1", "1496", "1496", "Clostridioides difficile", "Contig", "4000000", "28.000000", "3500", rep="Partial"),
    row("GCF_000840245.1", "10298", "10298", "Human alphaherpesvirus 1", "Complete Genome", "152222", "68.000000", "77", group="viral"),
    row("GCF_000000012.1", "12", "12", "Bernardetia litoralis", "Contig", "6000000", "36.000000", "5000"),
]
with open(summary, "w", newline="") as fh:
    fh.write("##  See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.\n")
    w = csv.DictWriter(fh, fieldnames=["#" + cols[0]] + cols[1:], delimiter="\t", lineterminator="\n"); w.writeheader()
    for r in rows: w.writerow({("#" + k if k == cols[0] else k): v for k, v in r.items()})
r = subprocess.run([sys.executable, "-I", os.path.join(here, "build.py"), summary, out], capture_output=True, text=True)
assert r.returncode == 0, r.stderr
def read(name): return list(csv.DictReader(open(os.path.join(out, name + ".tsv")), delimiter="\t"))
gs = read("genome_size"); gc = read("gc_content"); cg = read("coding_genes")
assert [(x["NCBI_ID"], x["Taxon_name"], x["Attribute_value"], x["Unit"]) for x in gs] == [
    ("562", "Escherichia coli", "4641652", "base pair"), ("1423", "Bacillus subtilis", "4215606", "base pair")], gs
assert gs[0]["Note"] == "GCF_000005845.2; Complete Genome; reference genome" and gs[0]["Attribute"] == "Genome size", gs[0]
assert gs[1]["Note"] == "GCF_000009045.1; Scaffold; assembly from type material", gs[1]
assert [(x["NCBI_ID"], x["Attribute_value"], x["Unit"]) for x in gc] == [("562", "50.8", "mol%"), ("1423", "43.5", "mol%")], gc
assert [(x["NCBI_ID"], x["Attribute"], x["Attribute_value"], x["Unit"]) for x in cg] == [("562", "genes", "4298", "genes")], cg
assert all(x["Attribute_source"] == "NCBI_RefSeq" and x["Frequency"] == "always" for x in gs + gc + cg)
assert "species taxid not in ncbi-taxon-db, dropped: 1" in r.stdout, r.stdout
# a file without the NCBI comment line is refused
open(summary, "w").write("x\n")
r = subprocess.run([sys.executable, "-I", os.path.join(here, "build.py"), summary, out], capture_output=True, text=True)
assert r.returncode != 0 and "does not look like" in r.stderr, r.stderr
print("build.py fixture test: OK")
