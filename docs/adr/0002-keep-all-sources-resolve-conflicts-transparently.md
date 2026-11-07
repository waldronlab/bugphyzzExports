# 0002. Keep all sources and resolve conflicts transparently

- **Status:** Accepted
- **Date:** 2026-11-05
- **Deciders:** Levi Waldron

## Context and Problem Statement

Sources disagree. Re-extracting Bergey's in 2026, importing Madin et al. 2020 in full and rebuilding BacDive from its API put three or more sources on many taxa, sometimes with different values. Dropping a source whenever it disagrees with another hides the disagreement and loses the records it alone holds. Madin 2020 also aggregates five resources that bugphyzz already imports directly (ProTraits, The Microbe Directory, KEGG, PATRIC, BacDive), so part of it repeats them.

## Decision

1. Every source is imported, and every value it holds stays in the curated data with `Attribute_source` (and, for aggregated sources, the upstream `Data_source`) recorded per row.
2. Disagreements are resolved at conflict resolution, by rules that are written down and logged, not by deleting rows upstream (#66).
3. The one exclusion is a true secondary duplicate: a Madin record whose upstream source is also imported directly (29,469 records that add no taxa or values), or an earlier extraction of a text that a later extraction of the same text supersedes (the 2020 Bergey's scrape, superseded by `Bergeys_2026`; waldronlab/bugphyzz#274).

## Alternatives Considered

- **Drop the lower-quality source for an attribute.** Rejected: errors are per record, not per source, and dropped rows can't be audited.
- **Keep Madin's duplicates.** Rejected: older snapshots of data already present, double-counting the same evidence in conflict resolution.

## Consequences

- Conflict resolution rules carry the weight and must be explicit and logged (#66); until then, ties between sources are broken alphabetically by source name.
- Users can see every source's claim for a taxon in `physiologies()`.

References: #66, #73, #77.
