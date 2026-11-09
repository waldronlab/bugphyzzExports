# Bergey's Manual, genus abstracts, 2026 extraction

Genus-level traits from the public abstracts of the genus chapters of Bergey's Manual of Systematics of Archaea and Bacteria (Wiley, ISBN 9781118960608), replacing the 2020 keyword scrape (bergeys.webscraper), which kept one state from "aerobic or facultatively anaerobic" and matched keywords in the wrong sentence for shape. Comparison against BacDive and Madin 2020 is in waldronlab/bugphyzz-review (`analysis/bergeys/`), summarised in waldronlab/bugphyzz#274. The export drops the 2020 rows (`Bergey's Manual` in the sheets).

## Provenance

- Chapter list and abstracts: Crossref API, `filter=isbn:9781118960608`, fetched 2026-11-04: 2,249 genus chapters (DOI prefix `10.1002/9781118960608.gbm`), all with abstracts. `dois.json` lists them; the abstract text is not committed, re-fetch it from Crossref.
- Extraction: `extract.py` with `system_prompt.txt`, model `gemini-3.8-flash` (the `gemini-flash-latest` alias on 2026-11-04, Google Generative Language API; `extract.py` now pins it), temperature 0, JSON output, 50 abstracts per call. Run 2026-11-04.
- Second pass, 2026-11-07: cell arrangement and spore shape with `system_prompt_arrangement.txt`, same model and settings, same abstracts.
- Frozen copy of the Crossref download, the batches and both runs' model outputs: the private repository [lwaldron/bergeys-bmsab-abstracts](https://github.com/lwaldron/bergeys-bmsab-abstracts), and on supermicro `~/data/bergeys_bmsab/` (`crossref_2026-11-04/`, `arrangement_2026-11-07/`, with checksums). The abstracts are Wiley's, under its text and data mining license, so they are not public; the tables here are.
- Genus to NCBI taxid: `genus_taxids.json`, made by `genus_taxids.py`. It keeps the 2026-11-04 map (NCBI taxdump of 2026-11-03, scientific names then synonyms at rank genus) for every title that map resolved, and looks up the rest with taxoniq in the pinned ncbi-taxon-db 2026.11.1 after cleaning "Incertae sedis", "Form-", "gen. nov.", quotes and the space after an initial; a hit counts only for a genus under Bacteria or Archaea. `SYNONYMS` in the script adds NCBI synonyms the lookup misses, including the ones NCBI stores with the authority appended (Brachymonas, now *Allobrachymonas*; Coenonia, now *Allocoenonia*; Candidatus Magnetobacterium). 2,094 of 2,117 distinct titles map. Of the rest, ten are section headings, eleven are names NCBI does not hold (the quoted *Desulfovibrio* split proposals, Borzia, Iyengariella, Levispirillum, Sporospirillum, Candidatus Nanopetraeus), Candidatus Haloredivivus is a clade rather than a genus in NCBI, and "Desulfoliva" gen. nov. is an NCBI synonym of *Desulfofaba*, which has its own chapter, so it is left unmapped on purpose (see the comment in `genus_taxids.py`). The check that found these is `analysis/bergeys/bergeys_synonyms.py` in waldronlab/bugphyzz-review.
- When several chapters map to one genus taxid (a 2015 original and its revision, or a genus NCBI has since merged), the newest chapter is used.
- `build_traits.py` writes the TSVs, applying `overrides.tsv`.

## Output (2026-11-09)

| Physiology | Rows | Genera |
|---|---|---|
| aerophilicity | 1,822 | 1,619 |
| gram stain | 1,624 | 1,561 |
| shape | 2,679 | 1,833 |
| motility | 1,690 | 1,517 |
| spore formation | 1,102 | 1,095 |
| growth temperature | 515 | 515 |
| gc content | 1,788 | 1,788 |
| arrangement | 1,671 | 857 |
| spore shape | 173 | 126 |
| width | 822 | 822 |
| length | 738 | 738 |

Each row carries the chapter DOI in `Note`. Frequency is the extractor's reading of the abstract: `always` for an unqualified statement, `usually` for "most", `sometimes` for "some" or "or". Dropped: 26 chapters whose title did not map to a genus taxid, 50 shape values of "other", 9 rows with a frequency outside the vocabulary. Curved rods are written as `vibrio`, added to `attributes.tsv` in bugphyzz 1.7.4. Growth temperature: only the stated optimum; abstracts that give a range without an optimum are not imported (the ranges are in bugphyzz-review `analysis/bergeys/bergeys2026_numeric.tsv`). G+C content (mol%, range as stated) is in `gc_content.tsv`. Width and length are the cell dimensions in µm as stated, a range or a single value. Spore shape is `elliptic`, `coccus` or `bacillus`; "endospore" is not a shape and is not extracted.

Compared with the 2020 scrape on genera both have: aerophilicity shares a value for 85%, gram stain 98%, shape 72% (the scrape's sentence errors), arrangement 88%, spore shape 89%; width and length ranges overlap for 95%. A hand check of 30 random genera against their abstracts found 28 right on every trait.

## Refreshing

```sh
# batches and outputs: clone lwaldron/bergeys-bmsab-abstracts (or rsync supermicro ~/data/bergeys_bmsab/) rather than re-fetching
python3 -I extract.py system_prompt.txt <batches_dir> <out_dir>   # needs ~/.pi/agent/auth.json google key
python3 -I extract.py system_prompt_arrangement.txt <batches_dir> <arrangement_out_dir>
sh ../ncbi_assemblies/python.sh genus_taxids.py genus_taxids.json dois.json genus_taxids.json   # after a taxonomy update
python3 -I build_traits.py <out_dir> genus_taxids.json . <arrangement_out_dir>
```

The batch files are the Crossref records split 50 per file.
