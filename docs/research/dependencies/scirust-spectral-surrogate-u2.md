# SciRust spectral-surrogate dependency pin for NoiseLab U2

Purpose: freeze the exact SciRust revision that satisfies the generic phase-randomized spectral-surrogate prerequisite declared by `docs/research/fluctuation-universality-u2-stage0.md`.

## Pinned dependency

- Repository: `Memorithm/scirust`
- Branch qualified for consumption: `master`
- Immutable commit: `0e2eaccac631b689f97c242c47bad11d433847d9`
- Merge provenance: SciRust pull request `#1422`, merged into `master`
- Required implementation path: `scirust-signal/src/surrogate.rs`
- Required API family: deterministic phase-randomized real-valued spectral surrogate with explicit seed and fail-closed input validation

This pin replaces the earlier preregistration-time candidate reference to PR #1421. A pull-request head is not an allowed scientific dependency for U2.

## Qualification boundary

This document records dependency provenance only. It does **not** execute Stage U2, inspect U2 outcomes, alter preregistered thresholds, change source-family selection, or claim that the SciRust primitive proves any NoiseLab hypothesis.

Before Stage U2 execution, NoiseLab code that consumes the surrogate primitive must itself be reviewed and qualified against this exact SciRust revision (or a newly preregistered replacement revision). If SciRust changes after this commit, U2 must not silently float to the newer revision.

## Frozen U2 invariants retained

The existing U2 preregistration remains authoritative for:

- the four source families;
- scales `1, 2, 4, 8, 16`;
- descriptor definitions and convergence score;
- shuffled empirical-marginal null;
- phase-randomized spectral null;
- minimum `199` surrogate repetitions;
- explicit seeds and one-sided `+1` empirical p-value correction;
- `alpha = 0.05` decision boundaries;
- preservation of all negative/protocol-failure outcomes;
- prohibition on using unrelated protected TDI final holdouts.

No scientific result is recorded by this pin.
