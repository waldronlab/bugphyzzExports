# bugphyzzExports news

Changes to the data pipeline and the exported files, newest first. Released data versions are on Zenodo (concept DOI 10.5281/zenodo.10980653); the `devel` branch holds the weekly export.

## 2026-11 (unreleased, on `devel`)

- **Export stops on unmapped attribute values** (#65): right after import, every discrete value must be in bugphyzz's `attributes.tsv`, after `physiologies()` has applied bugphyzz's synonym table (`attribute_synonyms.tsv`, waldronlab/bugphyzz#279). Unmapped values used to be dropped with only a line in `log_file.msg`, which lost all 14,992 BacDive shape rows. Habitat, disease association and antimicrobial resistance are imported once, up front, so the check runs before ASR.
- **Genus, species and strain only** (#64): the export no longer carries the 5,792 rows above genus (family to superkingdom). All of them were ASR on internal nodes of the LTP tree, and the deep ones are wrong with near-certain scores: Bacteria, Bacillota, Bacilli and Clostridia gram-negative, Bacteria at 58.6 °C. The `mixed` GMT files lose their family-to-phylum signatures.
- **Antimicrobial resistance names kept in the vocabulary** (#83): the AMR rewrite lowercased names and invented a few, so 1,372 PATRIC rows (kanamycin A, erythromycin A, beta-lactam, carbapenems, carbapenem-resistant *Acinetobacter*) were dropped by `filterData()`. Case is now kept, carbapenem-resistant *Acinetobacter* becomes resistance to carbapenems, and the rewritten names are checked against the vocabulary, which gains `resistance to beta-lactam` and `resistance to carbapenems` in bugphyzz 1.7.9.
- **Programmatic sources under version control** (`sources/`, #69, #70): a manifest plus one TSV per physiology for each scripted source, with the script and pinned upstream version, read by `bugphyzz::physiologies()` from bugphyzz 1.7.4. Hand curation stays in the Google Sheets.
  - `bergeys_2026`: sentence-level LLM extraction of all 2,249 Bergey's Manual genus-chapter abstracts (aerophilicity, Gram stain, shape, motility, spore formation, growth temperature optimum, G+C content), source name `Bergeys_2026`; the 2020 keyword scrape stays in the sheets for now (waldronlab/bugphyzz#274).
  - `madin_2020`: Madin et al. 2020 species traits in full (ten traits), excluding records Madin took from resources bugphyzz imports directly (#73, #77); upstream resource kept per row in `Data_source`.
  - `bacdive`: every BacDive strain from the public API (99,393 strains; #76, #78), in the layout of the former sheet, plus G+C content, cell length and width. Needs bugphyzz 1.7.5, which splits multi-valued fields.
  - `microbe_directory`: The Microbe Directory 2018 table (eleven traits, 3,658 taxa), matched to current taxids with taxoniq and NCBI synonyms. Replaces the hand-copied sheet rows, which had the taxid of another taxon on 619 of 1,451 plant pathogenicity rows, antimicrobial sensitivity values not in the source, and host-associated FALSE for every taxon. Source errors (*Parabacteroides* and *Moraxella* Gram positive, 121 °C for *Burkholderia gladioli*) are corrected in `overrides.tsv` (waldronlab/bugphyzz#260, #261). Adds spore formation.
  - `ncbi_assemblies`: genome size, G+C content and protein-coding gene count of one RefSeq assembly per bacterial or archaeal species (77,673 species; #74), source name `NCBI_RefSeq`; taxids checked against a pinned `ncbi-taxon-db` release with taxoniq.
- **New attribute `gc content`** (range, mol%) from Bergey's, Madin and BacDive (#71).
- **Validation thresholds applied at export** (#63, #67): imputed (ASR) annotations are kept only for attribute values with MCC >= 0.7 (discrete) or R2 >= 0.5 (numeric), and `rarely`/`never` states are dropped, so the exported files match the manuscript and `importBugphyzz()` defaults. The multistate file went from 172 MB to about 15 MB.
- **Weekly export restored** (#62, #68): the action had failed since 2026-03-29 because the Bioconductor devel image stopped shipping `remotes`; first clean run 2026-11-04.

## 2024-06-27, v1.0.3 (Zenodo 10.5281/zenodo.12574596)

- Last release before the changes above; the version used by `importBugphyzz()` by default.
