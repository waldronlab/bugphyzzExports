# Curation sheet snapshots

The bugphyzz hand-curation Google Sheets, snapshotted on 2026-11-24 by
`snapshot.sh` (waldronlab/bugphyzz#305, #113). The sheets are manual copies
of rows from the original sources named in `Attribute_source`, with
information added during curation; the sheets themselves are an intermediate
artifact, not a source. Each file is the CSV the sheet serves through
Google's publish-to-the-web link, as listed in bugphyzz
`inst/extdata/spreadsheet_links.tsv` and `spreadsheet_customlinks.tsv` at
commit e33fb4a. One file per physiology, spaces as underscores.

`snapshot.sh` is a one-time migration tool. Once bugphyzz reads these files
instead of the sheets, do not re-run it: it would overwrite the frozen
snapshot with whatever the sheets serve at that moment.

Rules:

- Rows copied from a cited source (the `Attribute_source` column) are not
  edited. A wrong row gets a drop in `../sheet_overrides.tsv` (ADR 0005);
  editing the copy would attribute to the source a claim it never made.
- These files are deliberately not in `../manifest.tsv`. bugphyzz will read
  them as the sheet data itself, with the `link` column of
  `spreadsheet_links.tsv` and `spreadsheet_customlinks.tsv` pointing here
  (waldronlab/bugphyzz#306; until that merges, bugphyzz and the export still
  read the live sheets). Listing them in the manifest would read every row
  twice.
- The sheets stay view-only with their publish-to-the-web CSV links
  unchanged, because bugphyzz 1.8.x reads those links until Bioconductor
  3.24 is retired (ADR 0007).

Format notes:

- The CSVs come from Google with CRLF line endings and no final newline. The
  script adds the final newline and `.gitattributes` has git store the files
  with LF; `read.csv()` gives the same result either way.
- Three files have a line break inside a quoted field, kept as the sheets
  have it: `hydrogen_gas_producing.csv` and `lactate_producing.csv` (taxon
  "Butyrate-producing bacterium T1-815") and
  `mutation_rate_per_site_per_year.csv` (five "Yersinia pseudotuberculosis"
  records). Those records span two lines in the file.

Licence: the curated data is CC BY 4.0, as bugphyzz data releases were before
Wikipedia-derived annotations were added (see the repository README); the
`Attribute_source` column names the underlying source of each row.

| File | Google Sheet |
|---|---|
| `acetate_producing.csv` | https://docs.google.com/spreadsheets/d/1PGW4C8fkyZYUJKA8yAW7ZjWEB6QXYdoTb6wrNxatK2c |
| `aerophilicity.csv` | https://docs.google.com/spreadsheets/d/1JgUzttAZf5Qu2obZ_2xXqsnEr_I1XKn_ZzNvrE6ssO4 |
| `antimicrobial_resistance.csv` | https://docs.google.com/spreadsheets/d/1ldVVQr6zT9g_Bf6WsJtAjke8yyP-gLp8llMYh6bipag |
| `arrangement.csv` | https://docs.google.com/spreadsheets/d/1qkNPdVDjpkdvZ-UtWeIrmWA3aYCInIQrseDDYqXiH-w |
| `butyrate_producing.csv` | https://docs.google.com/spreadsheets/d/1uPOIyeFD5ClPR2W2XyT6jAdD8FkgGuNBSl_8F2CT_uI |
| `coding_genes.csv` | https://docs.google.com/spreadsheets/d/1MVQ3o4fRYA4EOjp9mNAF0AvyQ0SzzdTMoZXruAIHmb4 |
| `disease_association.csv` | https://docs.google.com/spreadsheets/d/1MAtkI-UEc7T33BYvLBAkw4PDmEmNBroMVUi4KwN8cCc |
| `fatty_acid_composition.csv` | https://docs.google.com/spreadsheets/d/1u_zgX9dDuk2hIUMlfHqC2pDAmc3ktGf6Vj2EbzRiV8M |
| `genome_size.csv` | https://docs.google.com/spreadsheets/d/1LR0WnLf-trrGHLObcr4t2rs_4W4MHS21SZ8KoujF1k4 |
| `gram_stain.csv` | https://docs.google.com/spreadsheets/d/1ygWBlBczB24fyL0DTPmTModCuEDi-1j1_rsYBZ5JEOc |
| `growth_medium.csv` | https://docs.google.com/spreadsheets/d/1r7GIZWjBbmC-rRZYbHUwM31LGJ6YGexdZozYOvFZXVs |
| `growth_temperature.csv` | https://docs.google.com/spreadsheets/d/1yJjKSd59tjpq1Oaifd7S3jFBWFe2ZU-kmxV4HtYJznc |
| `habitat.csv` | https://docs.google.com/spreadsheets/d/16WVR910o0p_CyiqYxSa1hkDTdtsAM4Oy8TlbZNIr7ys |
| `health_associated.csv` | https://docs.google.com/spreadsheets/d/1oKvcKub-2NnbwYsIyXPETt7nUl89RwTk5Vwzk1WX0j0 |
| `host-associated.csv` | https://docs.google.com/spreadsheets/d/1SHk5e6ASkTDTg5CGKVEryBeHRPawXdtOoH_3Kge9Z48 |
| `hydrogen_gas_producing.csv` | https://docs.google.com/spreadsheets/d/14CNe6W0OmR5iY9lQU2IPFtYZsLd8a77ehR3lVJkx1A8 |
| `isolation_site.csv` | https://docs.google.com/spreadsheets/d/1QmwZcfoXhBdWPxTt6qa7hhMckHGV7c0lU2pTPpcSVfI |
| `lactate_producing.csv` | https://docs.google.com/spreadsheets/d/1Am4XPHK9zyEuPV1h1je2TmsSHLitJJOLPcuhN6-v_hk |
| `length.csv` | https://docs.google.com/spreadsheets/d/1gtD-kPb3XEQF0HAFcMMx3Deunj5GxVchPHXWHF-ZK9w |
| `mutation_rate_per_site_per_generation.csv` | https://docs.google.com/spreadsheets/d/1YIlMJCwqtzU1v2hTrPOPbjhxjmj_t-ol6WvMLkkz1bs |
| `mutation_rate_per_site_per_year.csv` | https://docs.google.com/spreadsheets/d/1K9KUScp1VeEkeXjJmKft_Rdr4xqlyD9Cx7eQpijzWMQ |
| `shape.csv` | https://docs.google.com/spreadsheets/d/1OpjfL-EIosX9iUg8nRxv03qpAx8iKo0FpcK_ZQokr-k |
| `sphingolipid_producing.csv` | https://docs.google.com/spreadsheets/d/1mVeLtc94SykCBBHhxxPhrbUgFopZ4OMWHrM_c5-grjU |
| `spore_shape.csv` | https://docs.google.com/spreadsheets/d/1XPtnoqeSfJFwDC3QHsXuNrTSJpgIZgUpv4ZmSFtT38g |
| `width.csv` | https://docs.google.com/spreadsheets/d/1Wc14Gp4nygdR6XNnZEJ5mZMJtCP6pniqa9GQOwzriNo |
