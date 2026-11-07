# Bergey's Manual, genus abstracts, 2026 extraction

Genus-level traits from the public abstracts of the genus chapters of Bergey's Manual of Systematics of Archaea and Bacteria (Wiley, ISBN 9781118960608), replacing the 2020 keyword scrape (bergeys.webscraper), which kept one state from "aerobic or facultatively anaerobic" and matched keywords in the wrong sentence for shape. Comparison against BacDive and Madin 2020 is in waldronlab/bugphyzz-review (`analysis/bergeys/`), summarised in waldronlab/bugphyzz#274. The export drops the 2020 rows (`Bergey's Manual` in the sheets).

## Provenance

- Chapter list and abstracts: Crossref API, `filter=isbn:9781118960608`, fetched 2026-11-04: 2,249 genus chapters (DOI prefix `10.1002/9781118960608.gbm`), all with abstracts. `dois.json` lists them; the abstract text is not committed, re-fetch it from Crossref.
- Extraction: `extract.py` with `system_prompt.txt`, model `gemini-3.8-flash` (the `gemini-flash-latest` alias on 2026-11-04, Google Generative Language API; `extract.py` now pins it), temperature 0, JSON output, 50 abstracts per call. Run 2026-11-04.
- Second pass, 2026-11-07: cell arrangement and spore shape with `system_prompt_arrangement.txt`, same model and settings, same abstracts.
- Frozen copy of the Crossref download, the batches and both runs' model outputs: the private repository [lwaldron/bergeys-bmsab-abstracts](https://github.com/lwaldron/bergeys-bmsab-abstracts), and on supermicro `~/data/bergeys_bmsab/` (`crossref_2026-11-04/`, `arrangement_2026-11-07/`, with checksums). The abstracts are Wiley's, under its text and data mining license, so they are not public; the tables here are.
- Genus to NCBI taxid: `genus_taxids.json`, made by `genus_taxids.py`. It keeps the 2026-11-04 map (NCBI taxdump of 2026-11-03, scientific names then synonyms at rank genus) for every title that map resolved, and looks up the rest with taxoniq in the pinned ncbi-taxon-db 2026.11.1 after cleaning "Incertae sedis", "Form-", "gen. nov.", quotes and the space after an initial; a hit counts only for a genus under Bacteria or Archaea. `SYNONYMS` in the script adds NCBI synonyms the lookup misses. 2,091 of 2,117 distinct titles map; most of the rest are section headings.
- When several chapters map to one genus taxid (a 2015 original and its revision, or a genus NCBI has since merged), the newest chapter is used.
- `build_traits.py` writes the TSVs, applying `overrides.tsv`.

## Output (2026-11-07)

| Physiology | Rows | Genera |
|---|---|---|
| aerophilicity | 1,819 | 1,617 |
| gram stain | 1,622 | 1,559 |
| shape | 2,675 | 1,830 |
| motility | 1,688 | 1,515 |
| spore formation | 1,100 | 1,093 |
| growth temperature | 514 | 514 |
| gc content | 1,786 | 1,786 |
| arrangement | 1,668 | 856 |
| spore shape | 173 | 126 |
| width | 820 | 820 |
| length | 736 | 736 |

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
