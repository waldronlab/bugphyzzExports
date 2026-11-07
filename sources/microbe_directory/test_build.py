#!/usr/bin/env python3
"""Fixture test for build.py: sh python.sh test_build.py (needs the pinned ncbi-taxon-db, no network)."""
import os, subprocess, sys, tempfile, csv
here = os.path.dirname(os.path.abspath(__file__))
d = tempfile.mkdtemp(); out = os.path.join(d, "out")
cols = ["microbe_id", "kingdom", "phylum", "class", "order", "family", "genus", "species", "gram_stain", "gram_stain_citation",
        "microbiome_location", "microbiome_location_citation", "antimicrobial_susceptibility", "antimicrobial_susceptibility_citation",
        "optimal_temperature", "optimal_temperature_citation", "extreme_environment", "extreme_environment_citation", "biofilm_forming",
        "biofilm_forming_citation", "optimal_ph", "optimal_ph_citation", "animal_pathogen", "animal_pathogen_citation", "spore_forming",
        "spore_forming_citation", "pathogenicity", "pathogenicity_citation", "plant_pathogen", "plant_pathogen_citation"]
def row(i, species, kingdom="Bacteria", **kw):
    r = dict.fromkeys(cols, ""); r.update(microbe_id=str(i), kingdom=kingdom, species=species, **kw); return r
rows = [
    row(1, "Escherichia coli", gram_stain="0", microbiome_location="1", pathogenicity="2", optimal_temperature="37.0", spore_forming="0"),
    row(2, "Prevotella copri", gram_stain="0", microbiome_location="0"),              # authority-only synonym -> 165179
    row(3, "Vibrio sp HENC 03", optimal_ph="7.5"),                                     # punctuation -> Vibrio sp. HENC-03, 992012
    row(4, "Bacteriovorax marinus", gram_stain="0"),                                   # species 97084 and its strain 862908 -> 97084
    row(5, "Streptococcus oligofermentans", gram_stain="0"),                           # synonym of S. cristatus 45634
    row(6, "Streptococcus cristatus", gram_stain="1"),                                 # same taxid: both rows kept
    row(7, "Escherichia virus T4", kingdom="Viruses", gram_stain="0"),                 # not Bacteria or Archaea
    row(8, "Nonexistus bacterium XYZ 1", gram_stain="1"),                              # no match
    row(9, "Ambiguus nomen", gram_stain="1"),                                          # synonym of two unrelated taxa
    row(10, "Bacillus subtilis", plant_pathogen="0", animal_pathogen="", biofilm_forming="1", extreme_environment="0"),
]
md = os.path.join(d, "md.csv")
with open(md, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
names = os.path.join(d, "names.dmp")
with open(names, "w") as fh:
    for line in [("165179", "Prevotella copri Hayashi et al. 2007", "synonym"), ("992012", "Vibrio sp. HENC-03", "scientific name"),
                 ("97084", "Bacteriovorax marinus", "synonym"), ("862908", "Bacteriovorax marinus", "synonym"),
                 ("45634", "Streptococcus oligofermentans", "synonym"),
                 ("562", "Ambiguus nomen", "synonym"), ("1423", "Ambiguus nomen", "synonym")]:
        fh.write(f"{line[0]}\t|\t{line[1]}\t|\t\t|\t{line[2]}\t|\n")
ov = os.path.join(d, "overrides.tsv")
with open(ov, "w") as fh:
    fh.write("physiology\tNCBI_ID\tTaxon_name\tAttribute\tAttribute_value\tFrequency\taction\treason\n")
    fh.write("gram stain\t45634\tStreptococcus cristatus\tgram stain negative\t\t\tdrop\ttest\n")
    fh.write("optimal ph\t992012\tVibrio sp. HENC-03\toptimal ph\t7\tunknown\tadd\ttest\n")   # add listed before its drop
    fh.write("optimal ph\t992012\tVibrio sp. HENC-03\toptimal ph\t7.5\t\tdrop\ttest\n")
    fh.write("gram stain\t1423\tBacillus subtilis\tgram stain positive\tTRUE\tunknown\tadd\ttest\n")
subprocess.run([sys.executable, os.path.join(here, "build.py"), md, names, out], check=True, env={**os.environ, "MD_OVERRIDES": ov}, stdout=subprocess.DEVNULL)
read = lambda f: list(csv.DictReader(open(os.path.join(out, f)), delimiter="\t"))
gram = {(r["NCBI_ID"], r["Attribute"]) for r in read("gram_stain.tsv")}
assert gram == {("562", "gram stain negative"), ("165179", "gram stain negative"), ("97084", "gram stain negative"), ("45634", "gram stain positive"), ("1423", "gram stain positive")}, gram
assert all(r["Frequency"] == "unknown" and r["Attribute_source"] == "The Microbe Directory" for r in read("gram_stain.tsv"))
assert [r["NCBI_ID"] for r in read("host-associated.tsv")] == ["562"]                 # microbiome_location 0 gives no row
assert [(r["NCBI_ID"], r["Attribute"]) for r in read("COGEM_pathogenicity_rating.tsv")] == [("562", "COGEM pathogenicity rating 2")]
assert [(r["NCBI_ID"], r["Attribute_value"]) for r in read("spore_formation.tsv")] == [("562", "FALSE")]
assert [(r["NCBI_ID"], r["Attribute_value"], r["Unit"]) for r in read("optimal_ph.tsv")] == [("992012", "7", "pH")]
assert [(r["NCBI_ID"], r["Attribute_value"]) for r in read("growth_temperature.tsv")] == [("562", "37.0")]
assert [(r["NCBI_ID"], r["Attribute_value"]) for r in read("plant_pathogenicity.tsv")] == [("1423", "FALSE")]
assert [(r["NCBI_ID"], r["Attribute_value"]) for r in read("biofilm_formation.tsv")] == [("1423", "TRUE")]
assert [(r["NCBI_ID"], r["Attribute_value"]) for r in read("extreme_environment.tsv")] == [("1423", "FALSE")]
assert not os.path.exists(os.path.join(out, "animal_pathogen.tsv"))                 # "" gives no row
un = {r["species"]: r["reason"] for r in read("unmatched.tsv")}
assert set(un) == {"Nonexistus bacterium XYZ 1", "Ambiguus nomen"}, un
assert un["Ambiguus nomen"].startswith("ambiguous"), un
## a stale override (matches nothing) and an unknown action both stop the build
for line in ["gram stain\t562\tEscherichia coli\tgram stain positive\t\t\tdrop\tstale\n", "gram stain\t562\tEscherichia coli\tgram stain positive\t\t\treplace\tbad\n"]:
    with open(ov, "w") as fh:
        fh.write("physiology\tNCBI_ID\tTaxon_name\tAttribute\tAttribute_value\tFrequency\taction\treason\n" + line)
    r = subprocess.run([sys.executable, os.path.join(here, "build.py"), md, names, os.path.join(d, "bad")], env={**os.environ, "MD_OVERRIDES": ov}, capture_output=True, text=True)
    assert r.returncode != 0 and ("matched nothing" in r.stderr or "unknown override action" in r.stderr), r.stderr
print("ok")
