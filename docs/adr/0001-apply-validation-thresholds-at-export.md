# 0001. Apply validation thresholds at export

- **Status:** Proposed
- **Date:** 2026-11-05
- **Deciders:** Levi Waldron, Claude Code

## Context and Problem Statement

Annotations imputed by ancestral state reconstruction (ASR) are only as good as their cross-validation in waldronlab/taxPProValidation. Until 2026-11, three thresholds were in use: `importBugphyzz()` defaulted to 0.8, the documentation said 0.5, and the manuscript used MCC 0.7 for discrete and R2 0.5 for numeric attributes. The exported files carried every ASR row, so the Zenodo release did not match the paper's counts and users got different data depending on how they loaded it.

## Decision

1. The export drops ASR rows whose attribute value fails MCC 0.7 (discrete) or R2 0.5 (numeric) at rank `all` in `validation_summary.tsv`, as shipped in the installed bugphyzz package.
2. ASR rows with no validation row are dropped, whatever the threshold.
3. Imputed "rarely" and "never" states are dropped.
4. `importBugphyzz()` uses the same values as defaults (`mccThreshold`, `r2Threshold`) and can only filter further.

## Alternatives Considered

- **Export everything and filter only in `importBugphyzz()`.** Rejected: the GMT and CSV files are used directly by people who never call the package, and the Zenodo files must match the paper.
- **A per-attribute threshold.** Not adopted now; a threshold picked after seeing a result is hard to defend (waldronlab/bugphyzz-manuscript#15).

## Consequences

- Release files are about a fifth of their former size and match the manuscript's counts.
- A new validation run changes the export: v1.0.5 removed aerobic, spirillum and COGEM rating 3 predictions and added G+C content.
- A new numeric attribute gets no ASR rows until taxPProValidation has a row for it.

References: #63, #67, waldronlab/bugphyzz#263, waldronlab/bugphyzz#271.
