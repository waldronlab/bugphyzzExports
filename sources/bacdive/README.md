# BacDive, strain records from the public API

Replaces the BacDive Google Sheet (a 2023 snapshot of 97,334 strains, 5,757 of them with an NCBI id) with every strain record in BacDive, each carrying the NCBI taxid BacDive assigns.

Citation: Schober I, Koblitz J, Sardà Carbasse J, et al. BacDive in 2025: the core database for prokaryotic strain data. Nucleic Acids Res. 2025;53(D1):D748-D756. doi:10.1093/nar/gkae959. License: CC BY 4.0 (https://bacdive.dsmz.de/about).

## Provenance

- API: https://api.bacdive.dsmz.de (no registration since February 2026). `fetch.py` calls `GET /fetch/<id;id;...>` for ids 1 to 180,000 in blocks of 50 (the API answers 404 for a block where no id exists; ids have long gaps) and stores the raw JSON gzipped in a scratch directory that is not committed. Four copies run in parallel on disjoint ranges (`run.sh`); a full fetch takes about 30 minutes.
- Fetched 5 November 2026; highest id with a record 176,393; 99,393 strains.
- NCBI ids: BacDive reports a taxid and whether it matches the strain or the species. A strain-level taxid goes in `NCBI_ID` (Rank strain), taken from the strain record or, as the former sheet did, from the strain's genome assembly record; the species taxid, from BacDive or from the NCBI taxdump of 3 November 2026 by species name, goes in `Parent_NCBI_ID`. Merged taxids are resolved. 4,854 strains have a strain-level id (1,387 of them from the genome record) and 97,213 a species-level one (the sheet: 5,757 and 94,553).
- Multi-valued fields such as oxygen tolerance ("anaerobe;facultative anaerobe" from two references) keep every value; bugphyzz 1.7.5 splits them into one row each. The sheet kept only the first.
- `build.py` writes the tables below. Nothing is edited by hand.

## Output

`bacdive.csv.gz`: one row per strain in the column layout of the former sheet, so `bugphyzz:::.getBacDive()` reads it unchanged. Multi-valued fields are joined with `;` in the same text forms the sheet used (for example `NaCl positive growth 0-5 %` for halophily, `glucose + fermentation` for metabolite utilization). `Culture_temperature` is the reported optimum when there is one, otherwise the first temperature with positive growth; every reported temperature entry is kept in `Culture_temperature_range`.

Narrow tables in the curation-sheet layout, read through `manifest.tsv`, for attributes the sheet never carried:

| File | Rows | Taxa | Notes |
|---|---|---|---|
| `gc_content.tsv` | 15,570 | | mol%, range as reported, outside 20 to 80 dropped |
| `length.tsv` | 5,932 | | cell length, um |
| `width.tsv` | 5,810 | | cell width, um |

Rows in the narrow tables use the strain taxid where BacDive has one, otherwise the species taxid.

## Refreshing

```sh
sh run.sh <scratch_dir> <ncbi_taxdump_dir>
```

Then update the date and counts above from the build output.
