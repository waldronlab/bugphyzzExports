# 0007. Hand curation moves from Google Sheets to `sources/sheets/`

- **Status:** Accepted
- **Date:** 2026-11-24
- **Deciders:** Levi Waldron
- **Amends:** 0003

## Context and Problem Statement

ADR 0003 kept hand curation in Google Sheets while moving scripted sources into this repository. The sheets are no longer edited (corrections go in `sheet_overrides.tsv`, ADR 0005), but `bugphyzz::physiologies()` still reads all of them live from Google on every call. That leaves two problems:

1. Nothing records what a sheet said when an export ran, so an export cannot be rebuilt.
2. A value cannot link back to the row it came from. A live sheet has no stable row address, but a line in a tagged file does, e.g. `https://github.com/waldronlab/bugphyzzExports/blob/vX.Y.Z/sources/sheets/aerophilicity.csv?plain=1#L123`. BugSigDB Explorer needs this for clickable provenance.

## Decision

1. Every curation sheet is snapshotted as-is into `sources/sheets/`, one CSV per physiology, by `sources/sheets/snapshot.sh`, a one-time migration tool. The sheets are manual copies of rows from cited original sources with curator additions, so the snapshot, like the sheets, is an intermediate artifact; the `Attribute_source` column carries the provenance. The files stay CSV so the reader in bugphyzz stays `read.csv()`.
2. `physiologies()` reads these files instead of Google (waldronlab/bugphyzz#306). This reverses item 4 of ADR 0003 ("Hand curation stays in Google Sheets") and restates item 2: `physiologies()` binds the manifest files to the rows of `sources/sheets/`, not to live sheet rows. `sources/sheets/` is not listed in `manifest.tsv`, or every row would be read twice.
3. Rows copied from a cited source are not edited in the files, for the reason in ADR 0005, which stays in force: an edited copy would misattribute the source. Corrections remain drops in `sheet_overrides.tsv`.
4. Rows that do not come from a cited source are corrected, and new curation is added, by pull requests editing the files, reviewed like any other change.
5. The Google Sheets are frozen (view-only, with a note pointing here) with their publish-to-the-web CSV links unchanged and no new rows or columns, because bugphyzz 1.8.x (Bioconductor 3.24) reads those links until 3.24 is retired. The `source_link` column of `spreadsheet_links.tsv` keeps the historical sheet URL for that reason and for provenance.

## Alternatives Considered

- **Keep reading the sheets live.** Rejected for the two problems above.
- **Convert to TSV like the programmatic sources.** Rejected: the snapshot's value is being a verbatim capture of what the sheets' CSV links serve, verifiable by direct comparison, and read by the same `read.csv()` call bugphyzz already uses, so the migration changes where the data is read, not how. A format conversion would make the files a transformation to be trusted rather than a capture to be checked, for no gain: these files are not in `manifest.tsv`, so no code benefits from them matching the TSV layout.
- **List the files in `manifest.tsv`.** Rejected: `physiologies()` reads them as the sheet data itself, so a manifest row would duplicate every row.
- **Edit corrections directly into the snapshot files.** Rejected in ADR 0005; a reorganisation of rows by original source is a separate, later decision.

## Consequences

- An export can be rebuilt from a tag, and every sheet-derived value has a stable line address for provenance links.
- Curation changes become pull requests: reviewed, dated, and attributed in git history.
- Sheet data reaches users through BiocFileCache like the manifest files, so a change shows up only after `forceDownload = TRUE` (waldronlab/bugphyzz#306).
- The frozen sheets linger until Bioconductor 3.24 is retired, and a drift between them and the snapshot would affect only bugphyzz 1.8.x.

References: waldronlab/bugphyzz#305, waldronlab/bugphyzz#306, #113, ADR 0003, ADR 0005.
