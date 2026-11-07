# 0005. Corrections to sheet rows are a versioned drop list applied at export

- **Status:** Proposed
- **Date:** 2026-11-07
- **Deciders:** Levi Waldron

## Context and Problem Statement

Many Google Sheet rows are copies of external sources (BacMap, ProTraits, KEGG, ...), and some are wrong: *Parabacteroides* Gram positive, a plasmid length given as a genome size. Editing the sheet loses the record of what the source said and why it was changed, and changing a value makes the sheet attribute to the source a claim it never made. ADR 0003 covers corrections to scripted sources (`overrides.tsv`) but not to sheet rows.

## Decision

1. Corrections to sheet rows go in `sources/sheet_overrides.tsv`: one row per correction, with physiology, NCBI_ID, Taxon_name, Attribute, optional Attribute_value, Attribute_source and a reason that names the issue.
2. The only action is drop. A correct value that no source gives is not added under another source's name.
3. The sheets are not edited, and `physiologies()` still returns the rows. The export drops them right after import (`inst/scripts/sheet_overrides.R`) and logs each drop.
4. A row that matches nothing stops the export, so corrections made obsolete by a sheet or source change are removed rather than silently ignored.

## Alternatives Considered

- **Edit the sheets.** Rejected: no review, no reason recorded with the data, and a changed value misattributes the source.
- **Apply the drops in `physiologies()`.** Rejected: ADR 0002 keeps every source's claim visible there.
- **Allow adds.** Deferred until a case needs one; it would need a curator source in `attribute_sources.tsv`.

## Consequences

- Every correction is reviewed in a pull request and dated in git.
- `physiologies()` and the export differ by exactly the rows in this file.

References: waldronlab/bugphyzz#260, waldronlab/bugphyzz#262, ADR 0002, ADR 0003.
