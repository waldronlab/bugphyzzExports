# Madin et al. 2020, condensed species traits

Madin JS, Nielsen DA, Brbic M, et al. A synthesis of bacterial and archaeal phenotypic trait data. Sci Data 7, 170 (2020). https://doi.org/10.1038/s41597-020-0497-4

- Upstream: `condensed_traits_NCBI.csv` from figshare, condensed_traits_NCBI v1.2.0, https://doi.org/10.6084/m9.figshare.13557827 (version 1, published 2021-01-12; MD5 `763522b14bbbd1e820215ea98312bc41`). 172,324 records, 21,500 NCBI species.
- Licence: CC BY 4.0, as stated on the figshare record. The same file, byte for byte, is `output/condensed_traits_NCBI.csv` at commit `396dcd4d` (2020-10-16) of https://github.com/bacteria-archaea-traits/bacteria-archaea-traits, but that repository has no licence file, so the import reads the figshare copy.
- `import.R` keeps, for multistate traits, every value held by at least 20% of a species' records (Frequency `always` if at least 90% agree, `usually` if at least 50%, else `sometimes`), so disagreement between Madin's upstream sources stays visible; for binary traits the majority only, ties dropped; for numeric traits the median per species with the record count in `Note`.
- Everything with a bugphyzz attribute is imported. `metabolism` (aerophilicity) duplicates the Madin_2020 rows pasted into the aerophilicity Google Sheet in 2022; remove those sheet rows once this file is in use. Not imported, no bugphyzz attribute yet: `doubling_h`, `range_tmp`, `isolation_source`, `pathways`, `carbon_substrates`.
- Vocabulary: gram positive/negative only (Madin has no "variable"); shape `vibrio` kept as `vibrio` (added to `attributes.tsv` in bugphyzz 1.7.4); `spiral` and `spirochete` as `spirillum`; `filament` as `filamentous`; `pleomorphic` as `pleiomorphic cell`; other shapes dropped.
- Madin pools about 25 upstream resources. Records from the five that bugphyzz imports directly (ProTraits, The Microbe Directory, KEGG, PATRIC, BacDive) are excluded so they are not counted twice; they are older snapshots of the same tables and add nothing. That removes 29,500 of 172,324 records. The upstream resource(s) behind each kept row are in the `Data_source` column (GOLD, Engqvist, GenBank, FAPROTAX, Pasteur, IJSEM/Fierer, ...), so conflict resolution can weigh them separately.

## Output (2026-11-05, after excluding directly imported resources)

| Physiology | Species |
|---|---|
| aerophilicity | 9,658 |
| gram stain | 11,545 |
| spore formation | 5,173 |
| motility | 6,554 |
| shape | 6,992 |
| genome size | 14,149 |
| coding genes | 3,539 |
| growth temperature (optimum) | 5,639 |
| optimal ph | 3,228 |
| gc content | 3,537 |

Refresh: `Rscript import.R` (downloads the figshare file and checks its MD5; set `MADIN_CSV=<local file>` to use a local copy).
