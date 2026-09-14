# U2 execution-contract audit

Audited NoiseLab revision: `a9e4c1bcaed803ea1ace9f34a01e87f9ecee1996`.

## Finding

The panel called `U2Readiness::preregistered().validate()`, but its actual public
`StageU2PanelConfig` was only checked for an alpha in `(0, 1]`. A caller could
change alpha, scales or seed roots while retaining `StageU2PanelMode::Scientific`
and receive `scientific_claim_permitted=true`. Validating a fresh default
readiness object does not validate the caller's configuration.

## Correction

`run_stage_u2_panel` now rejects configuration drift before generating residuals
or inspecting outcomes, in both scientific and non-scientific smoke modes:

- scales exactly `[1, 2, 4, 8, 16]`, including order and count;
- alpha exactly the existing `0.05` binary floating-point value;
- source and surrogate seed roots exactly the constants already chosen by the
  first executable panel implementation.

The alpha and scales come from the frozen U2 preregistration. The seed roots are
implementation-frozen choices, not numbers retroactively attributed to the
original preregistration. This correction preserves them; it does not search for
better seeds or change already generated observations.

Invalid non-finite/out-of-range alpha retains its prior error category. Other
drift reports the precise field through `PreregistrationDrift`. Public fields
remain for compatibility; modifying them no longer labels a custom experiment as
an execution of this canonical panel. A custom experiment needs its own protocol
and must use a separately identified execution path.

## Validation and scope

New regressions cover both modes, alternate finite thresholds, a one-ULP alpha
change, empty/reordered/duplicated/changed scales, both seed roots, and non-finite
or out-of-range alpha. Invalid cases exercise the public runner and return before
source generation. Canonical-mode validation does not execute scientific loads.
Existing deterministic smoke and negative-result paths remain unchanged.

No U0/U1 endpoint, result, descriptor, decision label, source definition, SciRust
pin, surrogate count, FieldLab interface, or holdout is changed. No full U2 panel
is run by this audit. Tests are accepted only from actual execution, not from
this document.

The legacy `scientific_claim_permitted` flag identifies selection of the frozen
scientific workload, not successful rejection of a null hypothesis. Protocol
failures, negative results and all six pairs must still be retained and reviewed.
It does not certify universality, novelty, a physical origin, or a cosmological
origin. The original exploratory multiplicity policy remains unchanged.
