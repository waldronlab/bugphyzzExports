# Bergey's Manual, genus abstracts, 2026 extraction

Genus-level traits from the public abstracts of the genus chapters of Bergey's Manual of Systematics of Archaea and Bacteria (Wiley, ISBN 9781118960608), replacing the 2020 keyword scrape (bergeys.webscraper), which kept one state from "aerobic or facultatively anaerobic" and matched keywords in the wrong sentence for shape. Comparison against BacDive and Madin 2020 is in waldronlab/bugphyzz-review (`analysis/bergeys/`), summarised in waldronlab/bugphyzz#274.

## Provenance

- Chapter list and abstracts: Crossref API, `filter=isbn:9781118960608`, fetched 2026-11-04: 2,249 genus chapters (DOI prefix `10.1002/9781118960608.gbm`), all with abstracts. `dois.json` lists them; the abstract text is not committed, re-fetch it from Crossref.
- Extraction: `extract.py` with `system_prompt.txt`, model `gemini-flash-latest` (Google Generative Language API), temperature 0, JSON output, 50 abstracts per call. Run 2026-11-04.
- Genus to NCBI taxid: NCBI taxdump of 2026-11-03, scientific names then synonyms at rank genus, after stripping "Incertae Sedis", "Form-", "gen. nov." and a stray space after the initial from the chapter titles. 2,063 of 2,117 distinct titles mapped; the rest are section headings.
- When a genus has more than one chapter (2015 original plus a revision), the newest is used.
- `build_traits.py` writes the TSVs, applying `overrides.tsv`.

## Output (2026-11-04)

| Physiology | Rows | Genera |
|---|---|---|
| aerophilicity | 1,940 | 1,598 |
| gram stain | 1,738 | 1,541 |
| shape | 2,872 | 1,798 |
| motility | 1,806 | 1,502 |
| spore formation | 1,184 | 1,098 |
| growth temperature | 799 | 799 |

Each row carries the chapter DOI in `Note`. Frequency is the extractor's reading of the abstract: `always` for an unqualified statement, `usually` for "most", `sometimes` for "some" or "or". Dropped: 72 chapters whose title did not map to a genus taxid, 54 shape values of "other", 5 rows with a frequency outside the vocabulary. "vibrio" (curved rod) is written as `bacillus` with "vibrioid" in `Note`, because `vibriod cell` is not in `attributes.tsv`. Growth temperature: the stated optimum when the abstract gives one (`always`), else the midpoint of the stated range (`sometimes`, range in `Note`). G+C content was extracted too (2,083 genera, in bugphyzz-review `analysis/bergeys/bergeys2026_numeric.tsv`) but bugphyzz has no attribute for it yet.

## Refreshing

```sh
python3 -I extract.py system_prompt.txt <batches_dir> <out_dir>   # needs ~/.pi/agent/auth.json google key
python3 -I build_traits.py <out_dir> <genus_taxids.json> .
```

The batch files are the Crossref records split 50 per file; see `analysis/bergeys/` in bugphyzz-review for the helper that makes them and the taxid map.
