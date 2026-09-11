# Bistable Langevin / Kramers calibration preregistration

Status: Stage 0 preregistration only. This document fixes the non-final protocol before implementation or inspection of experimental outcomes.

## Research question

Can NoiseLab recover the classical weak-noise, weak-periodic-forcing rate-matching regime of a symmetric overdamped bistable Langevin system, and can it correctly retain negative or inconclusive outcomes when the assumptions of the approximation are violated?

This calibration is not a search for a universal optimal noise level and is not evidence of a novel stochastic-resonance law.

## Model family

Use an overdamped one-dimensional Langevin system with a symmetric double-well potential,

```text
dx = -U'(x) dt + A cos(2 pi f t + phi) dt + sqrt(2 D) dW_t
U(x) = a x^4 / 4 - b x^2 / 2,   a > 0, b > 0.
```

The deterministic minima are at `x = +/-sqrt(b/a)` and the barrier is at `x = 0`. Under the stated convention, the zero-forcing barrier height is

```text
DeltaU = b^2 / (4 a).
```

Implementation must reuse SciRust RNG, integration, stochastic-process, signal and statistics primitives when suitable public APIs exist. Any missing generic numerical primitive belongs in SciRust first; NoiseLab should not create a competing generic implementation.

## Null hypothesis H0

Within a preregistered weak-forcing/adiabatic calibration regime, the measured response as a function of `D` has no reproducible interior optimum whose location is materially more consistent with the independently calculated Kramers half-period rate-matching prediction than with matched monotone/no-peak controls and grid-resolution uncertainty.

Failure to reject H0, absence of an interior peak, or disagreement with the rate-matching prediction must be retained as a valid negative or inconclusive calibration result.

## Independent oracle

For the unforced symmetric potential, use the overdamped Kramers approximation with a prefactor computed from the local curvatures at the minimum and barrier under the exact model convention used by the implementation. The oracle is calculated from model parameters, never fitted to the measured response sweep.

Define the predicted noise intensity `D_K` by solving

```text
r_K(D_K) ~= 2 f
```

when a positive solution exists under the preregistered parameter regime. The implementation must report the exact prefactor convention and units. If assumptions needed for the approximation are not met, the oracle is marked inapplicable rather than coerced into a prediction.

## Baselines and controls

At minimum compare:

1. `D = 0` deterministic forcing control.
2. A low-noise regime below the predicted matching region.
3. A high-noise regime above the predicted matching region.
4. The full preregistered logarithmic `D` grid, evaluated without adaptive retuning from observed outcomes.
5. A matched no-periodic-forcing (`A = 0`) control to distinguish stochastic resonance from noise-only switching/coherence effects.
6. A deliberately non-adiabatic or stronger-forcing regime in which the simple Kramers matching approximation is expected to degrade; this is a falsification/control regime, not a failure to hide.

Later colored-noise or non-Gaussian variants require separate preregistration and cannot be pooled with this white-noise calibration.

## Primary response metric

The primary metric is the response amplitude at the known forcing frequency, estimated with a fixed preregistered spectral estimator from the stationary analysis window. The exact estimator, windowing rule, transient discard and normalization must be fixed in code before outcome inspection.

Secondary diagnostics may include transition count, phase locking, residence-time distribution and signal-to-noise ratio, but they cannot replace the primary metric after results are observed.

## Replication and uncertainty

Use independent deterministic seeds generated through the pinned SciRust RNG API. The number of seeds, simulation duration, burn-in, sampling interval and `D` grid are fixed before running the measurement sweep.

Report per-`D` mean response and uncertainty across seeds. Candidate peak claims require a confidence interval or bootstrap interval whose resampling unit is the independent seed/trajectory, not individual correlated time samples.

No single-seed maximum can be classified as resonance.

## Decision rules

A sampled point may receive `BEST_SAMPLED_POINT` solely because it maximizes the declared finite-grid objective.

`INTERIOR_RESPONSE_PEAK` requires all of the following:

- the maximizer is not a boundary point of the preregistered `D` grid;
- its mean response exceeds both adjacent grid points by a preregistered minimum prominence;
- uncertainty analysis does not make the apparent peak indistinguishable from both neighbors under the declared criterion;
- the same conclusion is not produced by an implementation defect detected by deterministic/no-forcing controls.

`MECHANISM_COMPATIBLE` additionally requires that the measured peak location is consistent, within a preregistered tolerance combining simulation uncertainty and grid resolution, with the independently computed Kramers rate-matching region in the weak-forcing/adiabatic calibration regime.

A peak outside that region is not relabeled to preserve a positive result. Strong-forcing/non-adiabatic disagreement is retained and reported.

`CONFIRMED` is reserved for a later isolated confirmatory protocol and is out of scope for this Stage 0 document.

## Anti-leakage / provenance

- Record exact NoiseLab and SciRust revisions, compiler/toolchain, model parameters, seeds and configuration hashes.
- Do not change the `D` grid, prominence threshold, primary metric or acceptance tolerance after viewing the calibration outcomes.
- Exploratory follow-up must use a new protocol identifier and cannot overwrite this preregistration.
- Preserve negative, equivalent and inconclusive runs.
- Do not transfer any empirical claim to ProofLab as `PROVED`; only formalized mathematical obligations accepted by the configured Lean kernel can receive that status.

## Transfer gates

Only after this white-noise bistable calibration behaves consistently with its controls should NoiseLab proceed to FitzHugh-Nagumo coherence resonance. RoPE/FLAT-ATTENTION, KVLab or other complex-system perturbation studies remain downstream and cannot be used to tune this calibration.

## External methodological context

The classical rate-matching control follows the stochastic-resonance/Kramers literature already cited in `docs/RESONANT_OPERATING_POINTS.md`. Recent 2026 work on bistable and excitable systems shows that colored noise, interaction order and regime choice can shift or suppress response maxima; those results motivate explicit falsification regimes here but are not treated as evidence for NoiseLab outcomes.
