# U2 archived-surrogate score replay

This document specifies the narrow integrity contract implemented by
`verify_archived_surrogate_scores`.

The verifier operates on an already captured U2 bundle. It does not regenerate
source trajectories or surrogate transformations. From the frozen panel
configuration it reconstructs the preregistered job identities, validates the
complete surrogate-array index, decodes the exact retained IEEE-754 binary64
left/right arrays and recomputes the declared multiscale convergence score.
Every retained score must match the replayed score at the exact `f64::to_bits()`
level.

A successful pair must retain its complete two-null score set. A pair that
failed the existing U2 numerical/protocol path may have no retained scores; the
replay verifier does not invent missing observations or reinterpret that failure.
Duplicate job identities, malformed indexes, incorrect file sizes, non-finite
values, malformed completion markers and score-bit disagreement fail closed.

This is reproducibility and byte-level computation-integrity evidence only. It
does not recompute or authorize U2 p-values or panel decisions, does not change
any frozen source parameter, seed, scale, alpha or surrogate count, and does not
support a universality, resonance, common-mechanism or cosmological-origin claim.
The synthetic replay regressions are software fixtures, not scientific data.

The broader capture format, retention policy and integrity limitations are
documented in [`u2-input-capture.md`](u2-input-capture.md).
