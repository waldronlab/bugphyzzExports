# Wikipedia introductions, 2026 extraction (preliminary)

Genus and species traits from the introductory paragraphs of English Wikipedia articles, as carried by taxoniq (waldronlab/bugphyzz#296). The method is the one used for the Bergey's abstracts (`../bergeys_2026/`): an LLM reads each text and fills a fixed schema, and `build_traits.py` turns the result into curation-sheet tables.

Status: a first run for review. The source needs the `Wikipedia_2026` row in `bugphyzz/inst/extdata/attribute_sources.tsv` (branch `wikipedia-2026` of waldronlab/bugphyzz), with confidence `medium` and evidence `nas`. The export comparison below was run with these.

## Provenance

- Text: the `description` index of a taxoniq build with the Wikipedia stage, on supermicro at `~/taxoniq-upstream` (taxoniq v1.0.4, commit 60b3788, built 2026-10-30/31 from the NCBI taxdump of 2026-10-30; `wikipedia_extracts.json` md5 `dc63bc25d44c8b62f7cd7a867922e24c`, 232,207 records). The taxonomy-only release wheels (ncbi-taxon-db 2026.11.1) do not carry this index. Wikidata maps each article to an NCBI taxid; taxoniq stores the taxid, the article title and the extract.
- `collect.py` keeps the taxa under Bacteria or Archaea: 10,092 with an extract (7,603 species, 1,939 genera, 30 subspecies, 520 above genus), after skipping 4,140 taxids that the 2026-10-30 taxonomy no longer has and 18 empty extracts. It strips the HTML and writes 50 taxa per batch.
- Extraction: `extract.py` with `system_prompt.txt`, model `gemini-3.8-flash` (Google Generative Language API), temperature 0, JSON output. Run 2026-11-10, 202 batches; three batches came back one taxon short and were re-run.
- Second pass, 2026-11-10: halophily class, temperature class and pH class with `system_prompt_extremophily.txt`, same model and settings, on the 2,186 taxa whose text has a word for one of them (`select_batches.py`, 44 batches). Only the class words count; nothing is inferred from a habitat or a stated temperature.
- `system_prompt.txt` gained the habitat terms "skin" and "intertidal zone" after the run. The run returned them as "other: skin" and "other: tidal flat", which `build_traits.py` maps (`OTHER_HABITAT`).
- `build_traits.py` writes the TSVs, applying `overrides.tsv` (empty so far).
- The text, the batches and the model output are not committed, because the extracts are CC BY-SA. They are on supermicro and on Levi's laptop under `~/data/wikipedia_2026/run_2026-11-10/`.

## Licence

Wikipedia text is CC BY-SA 4.0. The tables here hold values read from the text and a link to the article, and none of the text. Once this source is in the export, the bugphyzz data files are distributed under CC BY-SA 4.0 instead of CC BY 4.0 (export file headers, the repository README and `.zenodo.json`). The software stays Artistic-2.0.

## Output (2026-11-10)

| Physiology | Rows | Taxa | Genus rows | Species rows |
|---|---|---|---|---|
| gram stain | 5,654 | 5,632 | 873 | 4,771 |
| habitat | 6,232 | 5,236 | 504 | 5,713 |
| aerophilicity | 3,365 | 3,321 | 546 | 2,812 |
| shape | 2,974 | 2,804 | 413 | 2,552 |
| motility | 2,641 | 2,628 | 347 | 2,290 |
| spore formation | 1,506 | 1,505 | 239 | 1,263 |
| host-associated | 821 | 814 | 146 | 667 |
| animal pathogen | 792 | 774 | 112 | 668 |
| arrangement | 274 | 193 | 87 | 186 |
| growth temperature | 189 | 189 | 11 | 177 |
| length | 130 | 130 | 17 | 113 |
| width | 125 | 125 | 20 | 105 |
| optimal ph | 110 | 110 | 5 | 105 |
| plant pathogenicity | 106 | 105 | 15 | 89 |
| gc content | 36 | 36 | 8 | 28 |
| spore shape | 22 | 21 | 4 | 17 |
| temperature class | 696 | 685 | 124 | 571 |
| halophily class | 562 | 554 | 67 | 493 |
| ph class | 254 | 250 | 39 | 214 |

Each row carries the article URL and the Wikidata item in `Note`. Frequency is the extractor's reading of the text, as for Bergey's.

What `build_traits.py` drops: taxa above genus (513 values); texts the model marked as not about the taxon (208, for example an article about a disease attached to the taxid of its agent); values whose trait has no supporting word in the text (176, of which 123 host-associated), a guard against the model answering from its own knowledge; 350 habitats outside the fixed vocabulary (sponges 25, coral 10, algae 9; left out on purpose, see below); 106 growth temperatures without a stated optimum; 20 shapes of "other".

A check of 30 random genus and species extractions against their text, made on the first half of the run, found every extracted value stated in the text, and one habitat ("tidal flat") left out.

## Halophily, temperature and pH classes

Three physiologies that are new to bugphyzz, added in bugphyzz 1.7.18 (`spreadsheet_links.tsv`, `attributes.tsv`). `physiologies()` returns them and the export includes them. Imputed values are kept only for classes that pass the cross-validation; `ph class` has values for fewer than 1% of the tips of the tree and is not imputed.

- `halophily class`: halophilic (441), halotolerant (120), non-halophilic (1). "Slight", "moderate" or "extreme" is in `Note`. BacDive's `halophily` is a different thing, the NaCl range of growth.
- `temperature class`: thermophile (284), mesophile (174), psychrophile (116), psychrotolerant (71), hyperthermophile (34), thermotolerant (17).
- `ph class`: alkaliphile (109), acidophile (81), alkalitolerant (37), neutrophile (22), acidotolerant (5).

The classes line up with the measured values other sources hold for the same taxa. Median growth temperature: psychrophile 17 °C (n 108), psychrotolerant 23 (67), mesophile 30 (137), thermotolerant 40 (14), thermophile 55 (231), hyperthermophile 85 (29). Median optimal pH: acidophile 4.1 (n 30), neutrophile 7.0 (9), alkalitolerant 8.3 (12), alkaliphile 9.0 (30). Against `extreme environment` in The Microbe Directory, counting the -phile classes other than mesophile and neutrophile as extreme: 149 extreme in both, 27 in neither, 12 extreme only in Wikipedia, 18 only in The Microbe Directory.

## Habitat vocabulary

Wikipedia habitats are written only with terms bugphyzz already has; the other 350 values (225 distinct, mostly used once) are dropped and not added to the vocabulary. Two were already there under another name or group: "tidal flat" is written as `intertidal zone`, and `skin` moved from the isolation-site group to habitat as well.

## Agreement with the other sources

`analysis/compare_sources.R`, results in `analysis/results/`. Strain records count for their species, and only taxa where both sources give a single state enter the percentage. The same measure for the other sources against BacDive is the yardstick.

| Attribute | Wikipedia vs BacDive (species) | Madin_2020 (medium) vs BacDive | BacMap or Microbe Directory (low) vs BacDive | Wikipedia vs Bergeys_2026 (genus) |
|---|---|---|---|---|
| aerophilicity | 95.5%, kappa 0.90, n 1,746 | 94.8%, 0.88 | BacMap 83.5%, 0.73 | 97.7%, 0.96, n 302 |
| gram stain | 98.7%, 0.97, n 3,231 | 97.8%, 0.96 | TMD 97.9%, 0.96; BacMap 95.9%, 0.91 | 98.8%, 0.97, n 484 |
| shape | 95.2%, 0.77, n 1,433 | 97.6%, 0.89 | BacMap 92.0%, 0.74 | 90.3%, 0.81, n 134 |
| motility | 95.5%, 0.91, n 1,442 | 97.0%, 0.94 | | 96.5%, 0.93, n 200 |
| spore formation | 99.4%, 0.98, n 673 | 99.6%, 0.99 | TMD 97.7%, 0.94 | 98.0%, 0.95, n 152 |

Numeric values against BacDive: growth temperature within 5 °C for 92% of 163 species (Spearman 0.92), G+C within 2 mol% for 86% of 22, length within 1 µm for all 22, width within 0.25 µm for 95% of 22. Optimal pH is within 0.5 of Madin_2020 for 39 of 40.

Wikipedia agrees with BacDive about as well as Madin_2020 does and better than BacMap. The two are not independent: many of the species articles are one-sentence stubs written from the original description, which is also where BacDive and Madin take much of their data. The 350 taxa where Wikipedia shares no state with BacDive or Bergeys_2026 are in `analysis/results/disagreements_with_high_confidence.tsv`; in the ones read, the model had copied the article faithfully and the article and BacDive differ (motile by gliding or twitching against BacDive's "non-motile", "facultative anaerobe" against a BacDive strain grown aerobically).

Pathogenicity and host association have no high-confidence reference. BacDive records animal pathogens only as TRUE, and all 201 shared species are TRUE in Wikipedia. Against The Microbe Directory, 289 shared taxa: 191 TRUE in both, 20 FALSE in both, 75 TRUE in Wikipedia and FALSE in The Microbe Directory, 3 the reverse. The 75 include *Neisseria meningitidis*, *Edwardsiella tarda*, *Borrelia turicatae* and the agent of chancroid, so here the error is in The Microbe Directory. Plant pathogenicity agrees with The Microbe Directory for 39 of 40.

## What it adds

`analysis/added.R` (curated data) and `analysis/compare_exports.R` (two full runs of `export_bugphyzz.R` on supermicro, with and without this source, bugphyzz 1.7.17 and taxPPro 1.0.1). The export runs predate the three classes, which the export does not include, and the `skin` and `intertidal zone` habitat terms (70 rows), so the habitat counts of the export are 70 rows short.

Curated data: 7,968 taxa get at least one row; 386 of them have no curated row of any physiology in bugphyzz (332 genera, 52 species, 2 below species) and 67 more are species known only through their strains. Of 26,485 rows, 15,475 are for a taxon with no curated value of that physiology at its own rank. Habitat is the largest part: 6,132 of 6,232 habitat rows are a term bugphyzz does not have for that taxon, and 4,541 of 5,236 taxa have no habitat at all. The three classes are new for every taxon.

Export, taxa with a value per attribute (baseline, with Wikipedia): habitat 3,875 to 8,441; motility 14,821 to 15,630; shape 16,873 to 17,680; host-associated 1,467 to 1,945; animal pathogen 2,776 to 3,174; gram stain 29,864 to 29,993; arrangement 2,028 to 2,140; length 6,625 to 6,735; width 6,589 to 6,691; aerophilicity 29,024 to 29,102; plant pathogenicity 2,205 to 2,268; optimal pH 4,004 to 4,062; spore formation 25,923 to 25,963. In total 7,784 new taxon by attribute pairs and 61 new taxa, nearly all curated values rather than ASR.

For aerophilicity, gram stain and spore formation, which have ASR in the export, most Wikipedia taxa already had a predicted value, so the effect is to replace predictions with stated values: 2,341 more curated gram stain rows and 1,078 fewer ASR rows. Of 146,059 taxon by attribute pairs present in both runs, 602 change their top state (0.41%). Where both runs predict by ASR the top state is the same for 99.1% (aerophilicity), 99.8% (gram stain) and 99.8% (spore formation) of taxa. The validation thresholds in `validation_summary.tsv` were not recomputed with this source.

At confidence `medium` Wikipedia outranks The Microbe Directory and BacMap and ties with Madin_2020 and ProTraits. taxPPro breaks a tie by source name, and Wikipedia wins it. Of the 602 changed pairs, 258 had a predicted or pooled value replaced by a stated one, 119 are ASR predictions that shift, and 225 are one curated value replacing another:

| Attribute | Was The Microbe Directory | Was Madin_2020 | Was ProTraits |
|---|---|---|---|
| animal pathogen | 76 (73 FALSE to TRUE, the cases above where Wikipedia is right) | | |
| plant pathogenicity | 1 | | |
| shape | | 38 | 7 |
| aerophilicity | | 29 | 6 |
| gram stain | | 23 | 4 |
| motility | | 17 | |
| spore formation | | 4 | |
| arrangement | | | 6 |
| host-associated | | | 2 |

In another 12 pairs the top state changes while the source stays Madin_2020 or ProTraits. 77 curated growth temperatures and 9 optimal pH values change (median 2.5 °C and 0.25 pH). No value from BacDive or Bergeys_2026 changes. In a run at confidence `low` the Madin_2020 and ProTraits columns were empty and 438 pairs changed (0.30%); the Microbe Directory column was the same.

## Known limits

- Coverage is what Wikidata links to an NCBI taxid; 4,140 linked taxids are no longer in the taxonomy and are skipped rather than remapped through merged taxids.
- An article can be edited by anyone and the extract is a snapshot of unknown date within the taxoniq build; the Note links the live article, not the revision.
- `physiologies()` leaves `Rank` empty for about half the rows because bugphyzz's internal `ranksParents` table predates these taxids (the same holds for Madin_2020). The export looks ranks up again, so nothing is lost there.
- "Curved rod" is written as `vibrio`, as in Bergeys_2026, and disagrees with `bacillus` in BacDive or Bergeys_2026 for 19 taxa.

## Refreshing

```sh
# on a machine with the Wikipedia-enabled taxoniq build
python3 -I collect.py <dir containing the ncbi_taxon_db package> <scratch_dir>
python3 -I extract.py system_prompt.txt <scratch_dir>/batches <scratch_dir>/out 8   # needs ~/.pi/agent/auth.json google key
python3 -I select_batches.py <scratch_dir>/records.jsonl <scratch_dir>/batches_extremophily 'halo|salt|salin|thermo|psychro|cryo|mesoph|acido|acid-|alkal|neutroph|cold|heat'
python3 -I extract.py system_prompt_extremophily.txt <scratch_dir>/batches_extremophily <scratch_dir>/out_extremophily 8
python3 -I build_traits.py <scratch_dir>/out <scratch_dir>/records.jsonl . <scratch_dir>/out_extremophily
# comparison; physiologies.rds is saveRDS(bugphyzz::physiologies("all")) without this source
Rscript analysis/compare_sources.R physiologies.rds . <scratch_dir>/records.jsonl analysis/results
Rscript analysis/added.R physiologies.rds . <scratch_dir>/records.jsonl analysis/results
Rscript analysis/compare_exports.R <baseline export dir> <export dir with this source> analysis/results
```
