# 0006. Impute whole attributes, judged by distance-matched kappa

- **Status:** Proposed
- **Date:** 2026-11-09
- **Deciders:** Levi Waldron
- **Amends:** 0001

## Context and Problem Statement

ADR 0001 keeps an imputed (ASR) value when that value's cross-validated MCC is at least 0.7. Judging each value separately has three problems:

1. **Skewed imputations.** Within one attribute, some states are imputed and others are not. In v1.1.0 every imputed aerophilicity value is "anaerobic", although held-out aerobic predictions were 92% correct, so users counting anaerobes see an inflated share.
2. **Unrepresentative test.** Cross-validation scores held-out curated taxa, which usually have a curated relative nearby. When few taxa are curated, the taxa actually imputed are much farther from one. For COGEM pathogenicity rating (525 curated LTP tips), the median Nearest Sequenced Taxon Index (NSTI, the distance on the tree to the nearest curated tip) is 0.073 for held-out curated tips and 0.151 for the imputed taxa, and held-out accuracy falls from 0.95 to 0.74 as NSTI grows. Its validation MCC of 0.74 (per value, unweighted, as in ADR 0001) therefore overstated the accuracy of its 22,000 imputations. On the 10,018 taxa that also had a curated biosafety level, imputed COGEM ratings agreed with it for 69%, while predicting rating 1 for every taxon would agree for 92%.
3. **No credit for near misses.** Mistaking a facultative anaerobe for an aerobe counts as the same error as mistaking an anaerobe for an aerobe.

## Decision

1. An attribute's states are imputed all together or not at all.
2. The decision uses one cross-validation score per attribute, Cohen's kappa:
   - unweighted for binary attributes and for attributes whose states have no order;
   - linearly weighted for attributes whose states form an ordered scale, so an error between adjacent states counts as half an error. The order is declared in bugphyzz's `attributes.tsv` before validation.

   Unlike MCC, kappa is lowered when the predicted states are distributed differently from the curated ones.
3. Held-out tips are weighted so that their NSTI distribution matches that of the Living Tree Project (LTP) tips to be imputed;  see "Covariate Shift Adaptation by Importance Weighted Cross Validation" (doi: 10.5555/1314498.1390324). The weight of a held-out tip is the share of tips to be imputed in its NSTI bin divided by the share of held-out tips in that bin, with bins at 0.02, 0.05, 0.1 and 0.2 and no clipping.  Bin shares are pooled over all folds.
4. An attribute is imputed if the lower 95% confidence bound of its score is at least 0.7 (R² at least 0.5 for numeric attributes, with the same weighting). The bound is a two-sided t-interval over the per-fold scores, with k − 1 degrees of freedom for k folds. The folds share one tree, so the interval may be too narrow. Imputations of attributes that fail are written to `bugphyzz_asr_dropped.csv.gz`.
5. Imputed rows keep their `Score`: the posterior probability of the state for discrete attributes, and a step function of NSTI for numeric attributes (0.9 up to NSTI 0.05, 0.5 up to 0.15, 0.1 up to 0.2, 0 beyond). The export drops imputed rows with Score below 0.5. In `importBugphyzz()`:
   - `asrScoreThreshold` (default 0.5, the export's own cut-off) filters imputed rows further, trading coverage within an attribute for precision; `excludeRarely` is deprecated in favour of it;
   - `kappaThreshold` (default 0.7) replaces the per-value `mccThreshold`, and `r2Threshold` (default 0.5) now applies per attribute;
   - thresholds below the export's defaults make `importBugphyzz()` load the imputations from `bugphyzz_asr_dropped.csv.gz` automatically, so `includeDroppedAsr` is removed.

   `bugphyzz_asr_dropped.csv.gz` holds every imputation left out of the main files: those of attributes that fail decision 4 and those with Score below 0.5. Per-attribute kappa and R² come from `validation_summary.tsv` in bugphyzz, as MCC does now.
6. The release notes report, for each imputed attribute, the precision of each state on the held-out folds, distance-matched, with the number of held-out predictions.

Ordered scales declared now: aerophilicity (aerobic, facultatively anaerobic, anaerobic), gram stain (negative, variable, positive), biosafety level (1, 2, 3), COGEM pathogenicity rating (1, 2, 3).

## Evidence

Five-fold held-out ASR with the export's pipeline (`fitMk` ER, `ancr`), 2026-11-09. Four random splits gave aerophilicity lower bounds of 0.783 to 0.786; the table shows the first split. Scores are distance-matched kappa, linearly weighted for ordered attributes; the lower 95% bound is in parentheses. The last column is the distance-matched MCC (multiclass R_K for more than two states) for comparison.

| Attribute | Curated LTP tips | Kappa | Decision | MCC |
|---|---|---|---|---|
| gram stain | 14,193 | 0.947 (0.939) | impute | 0.942 (0.934) |
| spore formation | 6,930 | 0.879 (0.875) | impute | 0.879 (0.875) |
| aerophilicity | 13,464 | 0.799 (0.786) | impute | 0.716 (0.700) |
| COGEM pathogenicity rating | 525 | 0.652 (0.551) | none | 0.670 (0.568) |
| motility | 10,934 | 0.619 (0.602) | none | 0.620 (0.603) |
| shape | 11,717 | 0.564 (0.530) | none | 0.568 (0.536) |
| arrangement | 374 | 0.529 (0.450) | none | 0.567 (0.489) |
| biosafety level | 10,709 | 0.532 (0.490) | none | 0.533 (0.489) |
| hemolysis | 421 | 0.207 (0.049) | none | 0.261 (0.113) |

Without the ordinal weighting, aerophilicity's lower bound is 0.699 to 0.703. Any threshold between 0.61 and 0.78 gives the decisions above.

The effective sample size of the weighted held-out tips is 97% to 98% of their number for aerophilicity and gram stain, 56% for COGEM pathogenicity rating, 67% for arrangement and 22% for hemolysis (largest weight 6.3). No fold had an empty bin.

## Alternatives Considered

- **Keep per-value MCC (ADR 0001).** Rejected for the skew and for the distance problem above.
- **Attribute-level MCC (R_K) for nominal and binary attributes, weighted kappa for ordered ones.** Not adopted: it gives the same decisions on the evidence above, mixes two scores under one threshold, and does not penalise a skewed mix of predicted states.
- **Minimum curated coverage (for example 5% of LTP tips).** Not adopted: the percentage would be chosen after seeing where COGEM falls, and distance-matching excludes COGEM without a new number.
- **Drop individual imputed rows beyond a fixed NSTI.** Not adopted: for the attributes that pass, accuracy changes little with NSTI (aerophilicity 0.87 at NSTI ≤ 0.02, 0.84 above 0.2).
- **A stricter default Score threshold to clean up weak minority states.** Not adopted as default: at Score ≥ 0.95, precision of imputed facultative anaerobes rises only from 0.59 to 0.68 (counting any curated state as correct) and of gram-variable from 0.22 to 0.46, while removing about half of anaerobic and aerobic imputations. Users can set `asrScoreThreshold`.

## Consequences

- Aerophilicity gains aerobic and facultatively anaerobic imputations (about 10,000 rows instead of about 1,200); gram stain gains "variable". COGEM pathogenicity rating loses all imputations (about 22,000 rows).
- Weak minority states are imputed with their attribute: facultatively anaerobic (precision about 0.59) and gram stain variable (about 0.22, 27 held-out predictions).
- The 0.7 threshold is kept from ADR 0001, where it was set for MCC. It is a default: on the evidence above and given the ordinal weighting, the decisions do not depend on it, and users can choose another (decision 5). For binary attributes |kappa| ≤ |MCC|, with equality when the predicted and curated distributions of states match.
- waldronlab/taxPProValidation must report distance-matched kappa (weighted for ordered attributes) per attribute and fold. Until it does, the export cannot apply this rule from `validation_summary.tsv`.
- `importBugphyzz()` arguments change as in decision 5. `excludeRarely` is in Bioconductor 3.23 and goes through the deprecation cycle; `mccThreshold` and `includeDroppedAsr` have not been released and are removed without one.
- Release counts no longer match the manuscript.

References: ADR 0001, #63, #67.
