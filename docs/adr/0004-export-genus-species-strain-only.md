# 0004. Export genus, species and strain only

- **Status:** Accepted
- **Date:** 2026-11-05
- **Deciders:** Levi Waldron

## Context and Problem Statement

ASR also assigns states to internal nodes of the LTP tree that map to families, orders and higher ranks. No curated annotation sits above genus, and the deep nodes take the root's state with near-certain scores (waldronlab/taxPPro#16), so those rows look confident and carry no information.

## Decision

The export keeps rows at genus, species and strain and drops the rest (5,792 rows in the 2026-11-05 data), logging the count dropped per attribute.

## Alternatives Considered

- **Keep higher ranks with a warning.** Rejected: users of the GMT files never see warnings.
- **Fix the scores at deep nodes first.** Deferred to waldronlab/taxPPro#16; this decision can be revisited if that work makes them meaningful.

## Consequences

- Signatures above genus are no longer published.
- If taxPPro#16 produces calibrated scores at higher ranks, revisit with a new ADR.

References: #64, #85.
