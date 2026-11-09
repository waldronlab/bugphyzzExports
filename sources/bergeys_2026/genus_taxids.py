"""Map Bergey's chapter titles to NCBI genus taxids.

Usage: sh ../ncbi_assemblies/python.sh genus_taxids.py <previous_map.json> <dois.json> <out.json>

The previous map (2026-11-04, NCBI taxdump of 2026-11-03, scientific names then
synonyms at rank genus) is kept for every title it resolved: it handles
synonyms and homonyms that taxoniq's name lookup does not. Titles it missed are
cleaned of "Incertae sedis", "Form-", "gen. nov.", quotes and the space after
an initial, then looked up with taxoniq in the pinned ncbi-taxon-db release;
SYNONYMS covers NCBI synonyms that the name lookup misses. A hit counts only
if it is a genus under Bacteria or Archaea.
"""
import json, re, sys
import taxoniq

DB_VERSION = "2026.11.1"
PROKARYOTES = {2, 2157}
# NCBI synonyms taxoniq's name lookup cannot resolve: cleaned title -> taxid
SYNONYMS = {
    "Catenococcus": "94139",  # homotypic synonym of Allocatenococcus (NCBI names page)
    # NCBI keeps these synonyms with the authority appended, which the exact
    # name lookup misses (names.dmp, new_taxdump 2026-11-01)
    "Brachymonas": "28219",  # "Brachymonas Hiraishi et al. 1995", now Allobrachymonas
    "Coenonia": "78328",  # "Coenonia Vandamme et al. 1999", now Allocoenonia
    "Candidatus Magnetobacterium": "40118",  # equivalent name of Magnetobacterium (SeqCode)
    # Not added: "Desulfoliva" gen. nov. (Galushko and Kuever 2020) is an NCBI
    # synonym of Desulfofaba (109169), which has its own chapter; the Desulfoliva
    # chapter describes only Desulfofaba fastidiosa, so mapping it would replace
    # the type species' traits. Left unmapped until NCBI accepts the genus.
}


def clean(title):
    s = re.sub(r"<[^>]+>", "", title)
    s = re.sub(r"[“”‘’\"']", " ", s)
    s = re.sub(r"(?i)\bincertae sedis\b(\s+[IVX]+(\.[IVX]+)?\s*\.?)?", " ", s)
    s = re.sub(r"(?i)\bform\s*[‐-]\s*", " ", s)
    s = re.sub(r"(?i)\b(gen|stat|sp)\s*\.\s*nov\s*\.?|\bcorrig\b|\bemendation\b", " ", s)
    s = re.sub(r"\b([A-Z]) (?=[a-z])", r"\1", s)
    s = re.sub(r"\s+", " ", s).strip()
    return re.sub(r"(?i)^candidatus\s*", "Candidatus ", s)


def lookup(name):
    try:
        t = taxoniq.Taxon(SYNONYMS[name]) if name in SYNONYMS else taxoniq.Taxon(scientific_name=name)
    except KeyError:
        return None
    if t.rank.name != "genus":
        return None
    if not {a.tax_id for a in t.lineage} & PROKARYOTES:
        return None
    return str(t.tax_id)


def main(prev_path, dois_path, out_path):
    prev = json.load(open(prev_path))
    titles = sorted({d["genus"] for d in json.load(open(dois_path))})
    out, added = {}, []
    for title in titles:
        p = prev.get(title) or {}
        if p.get("taxid"):
            out[title] = p
            continue
        name = clean(title)
        tid = lookup(name)
        out[title] = {"clean": name, "taxid": tid} if tid else {"clean": name, "taxid": None}
        if tid:
            added.append((title, name, tid))
    json.dump(out, open(out_path, "w"), ensure_ascii=False, indent=0, sort_keys=True)
    print(f"{sum(1 for v in out.values() if v.get('taxid'))} of {len(out)} titles mapped; "
          f"{len(added)} added by taxoniq {DB_VERSION}:")
    for a in added:
        print("  " + " -> ".join(a))
    print("unmapped:", sorted(t for t, v in out.items() if not v.get("taxid")))


if __name__ == "__main__":
    main(*sys.argv[1:4])
