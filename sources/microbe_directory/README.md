# The Microbe Directory

Eleven traits for 3,658 bacterial and archaeal taxa from The Microbe Directory, imported from its published 2018 table. This replaces the rows copied by hand into the Google Sheets, which carried copying errors, and the source's own errors are corrected in `overrides.tsv`.

Citation: Shaaban H, Westfall DA, Mohammad R, et al. The Microbe Directory: an annotated, searchable inventory of microbes' characteristics. Gates Open Res. 2018;2:3. doi:10.12688/gatesopenres.12772.1. Data under the MIT licence of the repository below. The website (microbe.directory) is gone.

## Provenance

- Upstream: `data/microbe-directory.csv` in https://github.com/microbe-directory/microbe-directory at commit `02fbf2633dfab0c418701197b6f869fa038d6776` (unchanged since 2018). 7,678 rows, 4,077 of them Bacteria or Archaea; viruses, viroids and eukaryotes are skipped. This is the version the bugphyzz sheets were copied from: matched by name, the sheet values agree with it, not with the later MicrobeDirectory 2.0 (github.com/dcdanko/MD2).
- Taxids: the table has names only. `build.py` matches them with taxoniq against ncbi-taxon-db 2026.11.1, then, for names that are not current, against `names.dmp` from NCBI's `new_taxdump_2026-11-01.zip` (the dump that release was built from), because taxoniq has no synonym lookup. The fallback ignores case and punctuation (the table writes `Frankia sp BMG5 12`), uses scientific names, equivalent names and synonyms, and strips the authority from old names (`Prevotella copri Hayashi et al. 2007`). A name must resolve to one current taxon; when it hits a species and that species' own strains, the species is taken. 2,939 names matched as current names and 912 through the fallback; 223 matched nothing and 3 were ambiguous (`unmatched.tsv`). 28 names share a taxid with another name (*Streptococcus oligofermentans* is now *S. cristatus*), and both keep their rows.
- Values, from the 0/1 codes:

  | Column | Physiology | Mapping |
  |---|---|---|
  | `gram_stain` | gram stain | 0 negative, 1 positive, 2 variable |
  | `animal_pathogen`, `plant_pathogen`, `antimicrobial_susceptibility`, `biofilm_forming`, `extreme_environment`, `spore_forming` | animal pathogen, plant pathogenicity, antimicrobial sensitivity, biofilm formation, extreme environment, spore formation | 1 TRUE, 0 FALSE |
  | `pathogenicity` | COGEM pathogenicity rating | 1, 2, 3 |
  | `microbiome_location` | host-associated | 1 TRUE; 0 gives no row, since not found in the human microbiome is not evidence of no host |
  | `optimal_temperature`, `optimal_ph` | growth temperature (Celsius), optimal ph | as published |

- Frequency is `unknown` on every row, as on the sheet rows: the table gives no frequencies. How much weight that carries is for conflict resolution (#66). `Note` keeps the name and code as published.
- `overrides.tsv` has two actions. `drop` removes rows matching `NCBI_ID` and `Attribute`, and `Attribute_value` too when it is filled in; a drop that matches nothing stops the build, so stale overrides show up. `add` appends a row. All drops run before any add. Current entries: Gram stain for three *Parabacteroides* and *Moraxella catarrhalis* (published positive; waldronlab/bugphyzz#260), the negative Gram stain that came in under *S. oligofermentans*, and five impossible numbers: 121 °C for *Burkholderia gladioli* (waldronlab/bugphyzz#261), −60 °C and −10 °C optima, and two negative pH values.

## What changed from the sheet rows

The sheets held 14,023 Microbe Directory rows. Compared with them, this import:

- Fixes the plant pathogenicity sheet, where 619 of 1,451 rows carried the taxid of a different taxon than the name (*Leptospira interrogans* under a *Burkholderia* id, plant pathogen TRUE).
- Replaces antimicrobial sensitivity: 588 sheet rows had no value in the source, and only 180 of the rest agreed with it.
- Replaces host-associated: the sheet had 1,884 FALSE rows (from `microbiome_location` 0, and some from 1) and no TRUE rows; here 722 TRUE rows. The 464 habitat rows (`colonization of host`) are not imported; host-associated carries the same information.
- Adds spore formation (1,982 rows), which never reached a sheet, and 187 to 967 taxa per trait that the sheet copy missed, mostly names it could not match.

Once this is on `devel`, delete the rows with `Attribute_source` "The Microbe Directory" from the 12 sheets (11 physiologies plus habitat), otherwise they are counted twice.

## Output

| File | Rows |
|---|---|
| `gram_stain.tsv` | 3,168 |
| `extreme_environment.tsv` | 2,782 |
| `plant_pathogenicity.tsv` | 2,212 |
| `animal_pathogen.tsv` | 2,038 |
| `spore_formation.tsv` | 1,982 |
| `growth_temperature.tsv` | 1,846 |
| `COGEM_pathogenicity_rating.tsv` | 1,451 |
| `optimal_ph.tsv` | 1,131 |
| `antimicrobial_sensitivity.tsv` | 1,086 |
| `host-associated.tsv` | 722 |
| `biofilm_formation.tsv` | 606 |

## Refreshing

```sh
sh run.sh <scratch_dir>      # downloads the table and the taxdump (about 170 MB, not committed), then builds
sh python.sh test_build.py   # fixture test, no network
```
