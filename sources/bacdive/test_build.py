#!/usr/bin/env python3
"""Fixture test for build.py: python3 -I test_build.py (needs no network, no taxdump download)."""
import json, gzip, os, subprocess, sys, tempfile, csv
here = os.path.dirname(os.path.abspath(__file__))
d = tempfile.mkdtemp(); raw = os.path.join(d, "raw"); tax = os.path.join(d, "tax"); out = os.path.join(d, "out")
os.makedirs(raw); os.makedirs(tax)
# minimal taxdump: 1386 Bacillus genus; 1423 B. subtilis species; 224308 B. subtilis 168 strain; 999 merged into 1423
open(os.path.join(tax, "nodes.dmp"), "w").write("1386\t|\t1\t|\tgenus\t|\n1423\t|\t1386\t|\tspecies\t|\n224308\t|\t1423\t|\tstrain\t|\n")
open(os.path.join(tax, "names.dmp"), "w").write("1386\t|\tBacillus\t|\t\t|\tscientific name\t|\n1423\t|\tBacillus subtilis\t|\t\t|\tscientific name\t|\n224308\t|\tBacillus subtilis subsp. subtilis str. 168\t|\t\t|\tscientific name\t|\n")
open(os.path.join(tax, "merged.dmp"), "w").write("999\t|\t1423\t|\n")
rec = {
  "1": {"General": {"BacDive-ID": 1, "NCBI tax id": {"NCBI tax id": 999, "Matching level": "species"}},
        "Name and taxonomic classification": {"species": "Bacillus subtilis", "strain designation": "168", "type strain": "yes"},
        "Morphology": {"cell morphology": [{"@ref": 1, "gram stain": "positive", "cell shape": "rod-shaped", "motility": "yes", "cell length": "2-3 µm"}, {"@ref": 2, "gram stain": "positive"}],
                       "colony morphology": {"@ref": 1, "type of hemolysis": "beta"}},
        "Physiology and metabolism": {"oxygen tolerance": [{"@ref": 1, "oxygen tolerance": "aerobe"}, {"@ref": 2, "oxygen tolerance": "facultative anaerobe"}],
                                      "spore formation": {"@ref": 1, "spore formation": "yes"},
                                      "halophily": [{"@ref": 1, "salt": "NaCl", "growth": "positive", "tested relation": "growth", "concentration": "0-7 %"}],
                                      "metabolite utilization": [{"@ref": 1, "metabolite": "D-glucose", "utilization activity": "+", "kind of utilization tested": "fermentation"}]},
        "Culture and growth conditions": {"culture temp": [{"@ref": 1, "growth": "positive", "type": "growth", "temperature": "15-50"}, {"@ref": 1, "growth": "positive", "type": "optimum", "temperature": "37"}]},
        "Safety information": {"risk assessment": {"@ref": 1, "biosafety level": "1", "biosafety level comment": "Risk group (German classification)"}},
        "Sequence information": {"GC content": {"@ref": 1, "GC-content": "43.5"}, "Genome sequences": [{"accession": "GCA_000009045", "NCBI tax ID": 224308}]},
        "Isolation, sampling and environmental information": {"isolation": {"@ref": 1, "sample type": "soil", "country": "Germany"}}},
  "2": {"General": {"BacDive-ID": 2, "NCBI tax id": {"NCBI tax id": 1423, "Matching level": "species"}},
        "Name and taxonomic classification": {"species": "Bacillus subtilis", "type strain": "no"},
        "Sequence information": {"GC content": {"@ref": 1, "GC-content": "95"}}}}
with gzip.open(os.path.join(raw, "000001-000050.json.gz"), "wt") as fh: json.dump(rec, fh)
# missing chunks must abort
r = subprocess.run([sys.executable, "-I", os.path.join(here, "build.py"), raw, tax, out, "100", "50"], capture_output=True, text=True)
assert r.returncode != 0 and "missing" in r.stderr, r.stderr
r = subprocess.run([sys.executable, "-I", os.path.join(here, "build.py"), raw, tax, out, "50", "50"], capture_output=True, text=True)
assert r.returncode == 0, r.stderr
rows = {x["BacDive_ID"]: x for x in csv.DictReader(gzip.open(os.path.join(out, "bacdive.csv.gz"), "rt"))}
a = rows["1"]
assert a["Taxon_name"] == "Bacillus subtilis 168" and a["Parent_NCBI_ID"] == "1423", a   # merged 999 -> 1423
assert a["NCBI_ID"] == "224308", a                      # strain taxid from the genome record
assert a["Oxygen_tolerance"] == "aerobe;facultative anaerobe" and a["Gram_stain"] == "positive", a
assert a["Culture_temperature"] == "37" and "growth positive 15-50" in a["Culture_temperature_range"], a
assert a["Halophily"] == "NaCl positive growth 0-7 %" and a["Metabolite_utiilization"] == "D-glucose + fermentation", a
assert a["Hemolysis"] == "beta" and a["Biosafety_level"] == "1" and a["Country"] == "Germany" and a["Genome_ID"] == "GCA_000009045", a
assert rows["2"]["NCBI_ID"] == "" and rows["2"]["Parent_NCBI_ID"] == "1423", rows["2"]
gc = list(csv.DictReader(open(os.path.join(out, "gc_content.tsv")), delimiter="\t"))
assert [(g["NCBI_ID"], g["Attribute_value"]) for g in gc] == [("224308", "43.5")], gc      # 95 mol% dropped
ln = list(csv.DictReader(open(os.path.join(out, "length.tsv")), delimiter="\t"))
assert ln[0]["Attribute_value"] == "2-3" and ln[0]["Unit"] == "um", ln
print("build.py fixture test: OK")
