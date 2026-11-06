# 0003. Programmatic sources live in `sources/`, hand curation in Google Sheets

- **Status:** Proposed
- **Date:** 2026-11-05
- **Deciders:** Levi Waldron, Claude Code

## Context and Problem Statement

Some curation data comes from scripts (the Bergey's extraction, Madin 2020, the BacDive API) and some from curators. Script output used to be pasted into the same Google Sheets as hand curation. A re-run then overwrote hand fixes, and nothing recorded how the data was made or let anyone repeat it. The output is too large for the Bioconductor package tarball.

## Decision

1. Each programmatic source has a directory under `sources/` with the script that produced it, its pinned input version, one TSV per physiology in the sheet column layout, and a README with provenance and licence.
2. `sources/manifest.tsv` lists the files; `bugphyzz::physiologies()` reads them from the `devel` branch through BiocFileCache and binds them to the sheet rows.
3. Output TSVs are never edited by hand. Corrections go in the source's `overrides.tsv`, applied on every re-run.
4. Hand curation stays in Google Sheets.
5. Taxonomy lookups use a dated taxoniq `ncbi-taxon-db` release, recorded in the README.

## Alternatives Considered

- **Keep pasting into Sheets.** Rejected for the reasons above.
- **Ship the data in the bugphyzz package.** Rejected: size, and each data change would need a package release.

## Consequences

- Sources can be re-run and diffed; fixes survive re-runs.
- `physiologies()` needs network access to this repository (falls back to sheets only, with a warning).

References: #69, #70, #78, waldronlab/bugphyzz#275.
