# Bistable Langevin / Kramers calibration — Stage 0 v2

Status: PREREGISTERED PROTOCOL — executable Stage 0 v2 runner available. This document replaces the invalidated convention-mismatched preregistration. The panel can be executed via `examples/bistable_kramers_stage0_report.rs` (`src/bistable_stage0.rs`). Smoke (`NOISELAB_BISTABLE_SMOKE` / default) is not scientific evidence; full load requires `NOISELAB_BISTABLE_FULL=1`. No Stage 0 v2 scientific results.md is claimed by shipping the runner alone.

## Provenance and implementation pin

Historical pin named in the original freeze text: NoiseLab commit `d96b042fceb47210078aad441be548babe6a78da` (not present in the current clone history). The executable model convention remains that of `src/langevin.rs`.

### Pin-sync note (requalification)

The runner and frozen factors in `src/preregistered.rs` / `src/bistable_stage0.rs` are requalified against current `main` (`src/langevin.rs`, SciRust pin in `Cargo.toml`) without forging outcomes under the missing historical SHA. Re-run the scientific panel under a retained producing commit before any calibration claim. Past references to `d96b042…` remain historically accurate as the original intended freeze identity, not as a reproducible object in this repository.

The model convention is fixed to the implementation, not inferred from the invalidated Stage 0:

```text
dx = (a*x - b*x^3 + A*cos(2*pi*f*t + phi))*dt + sqrt(2*D)*dW
U(x) = b*x^4/4 - a*x^2/2
```

with `a > 0`, `b > 0`. Therefore:

```text
xmin = +/-sqrt(a/b)
DeltaU = a^2/(4b)
r0 = a*sqrt(2)/(2*pi)
```

for the unforced symmetric potential and the unit-mobility weak-noise overdamped Kramers prefactor implemented by NoiseLab.

The previously invalidated document remains immutable provenance and MUST NOT be edited into agreement with this protocol.

## Research question

Under a fixed weak/subthreshold periodic forcing regime that has not been used by the existing executable calibration test, does the mean phase-locked well-occupancy response exhibit a reproducible interior maximum whose location is compatible with the independently computed Kramers half-period matching control?

This is a calibration of a known mechanism. It is not a search for a universal optimum and cannot establish novelty from a numerical maximum.

## Null hypothesis H0

On the preregistered noise grid, the mean coherent switching response has no uncertainty-separated interior maximum whose location is compatible with the independently calculated weak-noise Kramers rate-matching region. Boundary maxima, monotone responses, broad plateaus, incompatible peak locations and statistically unresolved maxima all fail to reject H0 and must be retained.

## Frozen model parameters

Use a parameter point distinct from the existing engineering regression test:

```text
a = 1.2
b = 1.0
forcing_amplitude A = 0.15
forcing_frequency_hz f = 0.04
phase_rad = 0.0
```

These values are fixed before any Stage 0 v2 sweep. The implementation's static subthreshold check must return true; otherwise the protocol fails closed and no sweep is interpreted.

The independent analytic control is computed before simulation from:

```text
DeltaU = a^2/(4b)
r0 = a*sqrt(2)/(2*pi)
rK(D) = r0 * exp(-DeltaU/D)
```

and the half-period matching equation

```text
rK(D_K) = 2*f.
```

For the frozen parameters above, the algebraic control gives `DeltaU = 0.36`, `r0 ~= 0.27009489484713184`, and `D_K ~= 0.2958709422522067`. These are analytic preregistration values, not empirical outcomes.

## Frozen noise grid

Use the following multiplicative factors around the analytic control, with no adaptive insertion after observing responses:

```text
[0.25, 0.40, 0.63, 1.00, 1.58, 2.50, 4.00] * D_K
```

The exact floating-point `D` values used by the runner must be derived once from the analytic `D_K` above and recorded in the evidence artifact before trajectories are generated.

Also run the deterministic `D = 0` control separately. It is not eligible to become the interior-peak candidate.

## Frozen simulation budget

Use:

```text
dt = 0.02 s
80 forcing periods total
20 forcing periods burn-in
16 independent seeds
initial state = -sqrt(a/b)
```

The 16 seeds are fixed to:

```text
[101, 211, 307, 401, 503, 601, 701, 809,
 907, 1009, 1103, 1201, 1301, 1409, 1511, 1601]
```

The same seed set is used at every `D` so comparisons are paired by seed. No seed may be dropped because its trajectory weakens the apparent effect; simulation failures are recorded as protocol failures.

## Primary metric

The primary response is the existing NoiseLab `coherent_switching_response`: map each post-burn-in state to well occupancy `-1/+1` and measure the magnitude of its cosine/sine projection at the declared forcing frequency.

The experimental unit for uncertainty is the independent seed/trajectory. Individual time samples are never treated as independent replicates.

For each `D`, report:

- mean coherent switching amplitude across the 16 seeds;
- sample standard deviation across seeds;
- standard error `s/sqrt(16)`;
- all per-seed amplitudes in the evidence artifact.

## Controls

The same code path must execute:

1. `D = 0` with all declared seeds, verifying seed independence of the deterministic trajectory.
2. `A = 0` over the identical positive-`D` grid to distinguish periodic phase locking from noise-only switching.
3. The full frozen weak/subthreshold grid above without adaptive retuning.
4. A separate falsification regime with `forcing_frequency_hz = 0.20` while holding `a`, `b`, `A`, phase, `dt`, periods, burn-in, seeds and grid factors fixed. Its Kramers control is recomputed analytically. This regime is not pooled with the primary calibration and exists to test degradation of simple adiabatic rate matching.

No colored-noise, non-Gaussian-noise, FHN, attention or KV result is part of this protocol.

## Decision rules

`BEST_SAMPLED_POINT` may identify the finite-grid maximizer without implying resonance.

`INTERIOR_RESPONSE_PEAK` requires:

- the maximizer to be one of grid indices 1 through 5, never an endpoint;
- its mean to exceed both immediate neighboring grid means;
- its mean to exceed each scan-edge mean by at least two pooled standard errors, using seed-level standard errors;
- deterministic and `A = 0` controls to show no implementation failure that explains the maximum.

`MECHANISM_COMPATIBLE` additionally requires the interior peak coordinate to lie within a factor of 2 of the independently calculated `D_K`. This deliberately coarse tolerance reflects the weak-noise asymptotic nature of the Kramers control and the finite logarithmic grid. It may not be widened after outcome inspection.

If the primary weak/subthreshold regime lacks an `INTERIOR_RESPONSE_PEAK`, H0 is not rejected. If an interior peak exists but fails the factor-of-2 Kramers compatibility rule, report an empirical peak without Kramers mechanism compatibility. If the falsification regime happens to agree as well or better, report that observation rather than suppressing it.

`CONFIRMED` is outside this Stage 0 and requires a separate, later protocol.

## Anti-leakage and evidence

Before generating trajectories, the runner must record:

- NoiseLab revision;
- SciRust revision resolved by Cargo;
- compiler/toolchain;
- full model/run parameters;
- analytic `DeltaU`, `r0`, `D_K` and derived positive-`D` grid;
- exact seed list;
- source hash of this preregistration.

The primary grid, seeds, model parameters, burn-in, metric and decision thresholds may not change after any Stage 0 v2 response is inspected. Any exploratory follow-up receives a new protocol identifier and cannot overwrite this evidence.

Negative, equivalent and inconclusive results are retained.

## Executable runner

Library entry points:

- `BistableStage0V2::primary()` — outcome-blind materialization of frozen model, run, seeds, analytic `D_K`, and positive-`D` grid;
- `run_bistable_kramers_stage0` — panel including `D = 0`, primary forced grid, `A = 0` control, falsification preflight, and protocol-language classification;
- `classify_bistable_stage0` — fail-closed decision helper (`H0NotRejected`, interior peak, mechanism compatibility, control/protocol failure).

Report executable: `cargo run --example bistable_kramers_stage0_report`.

Environment (same `report_mode` contract as FHN/U2):

- `NOISELAB_BISTABLE_FULL` — scientific budget when `1`/`true`;
- `NOISELAB_BISTABLE_SMOKE` — explicit smoke when `1`/`true`;
- unset both → non-scientific smoke;
- both enabled → error before any trajectory.

Smoke abbreviates to `4` total / `1` burn-in forcing periods and the first `2` preregistered seeds. The positive-`D` multiplicative factors and analytic `D_K` are not retuned. The falsification regime at `f = 0.20` remains fail-closed when the analytic Kramers control has no positive finite solution (current frozen parameters).

## Transfer gates

Only after this calibration has a complete evidence artifact and its controls have been evaluated may its methodological lessons be considered for the planned FitzHugh-Nagumo coherence-resonance stage. No empirical result from this protocol is promoted to ProofLab as `PROVED`; formal proof status remains exclusively determined by the configured Lean kernel.
