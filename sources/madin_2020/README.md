# Madin et al. 2020, condensed species traits

Madin JS, Nielsen DA, Brbic M, et al. A synthesis of bacterial and archaeal phenotypic trait data. Sci Data 7, 170 (2020). https://doi.org/10.1038/s41597-020-0497-4

- Upstream: https://github.com/bacteria-archaea-traits/bacteria-archaea-traits, `output/condensed_traits_NCBI.csv`, commit `396dcd4d` (2020-10-16). 172,324 records, 21,500 NCBI species.
- `import.R` takes the majority value per species for discrete traits (Frequency `always` if at least 90% of records agree, `usually` if at least 50%, else `sometimes`) and the median per species for numeric traits.
- Everything with a bugphyzz attribute is imported. `metabolism` (aerophilicity) duplicates the Madin_2020 rows pasted into the aerophilicity Google Sheet in 2022; remove those sheet rows once this file is in use. Not imported, no bugphyzz attribute yet: `gc_content`, `doubling_h`, `range_tmp`, `isolation_source`, `pathways`, `carbon_substrates`.
- Vocabulary: gram positive/negative only (Madin has no "variable"); shape `vibrio` written as `bacillus`; `spiral` and `spirochete` as `spirillum`; `filament` as `filamentous`; `pleomorphic` as `pleiomorphic cell`; other shapes dropped.
- About half of Madin's records come from GOLD and PATRIC, which overlap other bugphyzz sources; conflict resolution in the export handles the overlap.

## Output (2026-11-04)

| Physiology | Species |
|---|---|
| aerophilicity | 10,886 |
| gram stain | 12,246 |
| spore formation | 6,884 |
| motility | 7,524 |
| shape | 8,023 |
| genome size | 14,861 |
| coding genes | 3,608 |
| growth temperature (optimum) | 7,023 |
| optimal ph | 3,928 |

Refresh: `Rscript import.R` (reads the pinned upstream CSV; set `MADIN_CSV=<local file>` to use a local copy).
