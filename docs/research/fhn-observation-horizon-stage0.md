# FitzHugh–Nagumo observation-horizon robustness — Stage 0

## Status

Preregistered robustness control. No outcome has been observed under this protocol and no novelty claim is authorized.

This control extends the existing FitzHugh–Nagumo coherence-resonance calibration without changing its model, noise grid, seeds, spike definition, uncertainty rule, or acceptance interval. Its sole purpose is to test whether the sampled coherence-resonance decision is stable to observation duration rather than being an artefact of a single finite horizon.

## Frozen dependency and parent protocol

Parent protocol: `docs/research/coherence-resonance.md`.

The parent protocol pins SciRust revision:

```text
f57d598bf03e5dfb16ec6423e4e43105a77540d3
```

This robustness control must use the same dependency revision and the same FHN implementation as the parent protocol. A dependency or implementation change invalidates direct comparison and requires a new Stage 0.

## Research question

Does the parent protocol's qualitative decision — presence or absence of an uncertainty-separated sampled interior minimum of `CV_ISI` — remain unchanged when the post-burn-in observation horizon is varied while every other declared experimental factor is held fixed?

This is a finite-horizon robustness question. It is not a test of asymptotic stationarity and it does not establish a universal observation time.

## Null hypothesis

H0: under the declared parameterization, the parent decision is not robust to observation horizon; at least one preregistered horizon changes the decision class or moves an accepted minimum outside the parent's accepted coordinate interval.

A negative or horizon-sensitive outcome must be retained. It must not be repaired by adding horizons, moving the noise grid, changing seeds, weakening the uncertainty rule, or discarding inconvenient replicates.

## Frozen experimental factors

Use the parent values unchanged:

```text
epsilon = 0.01
a = 1.05
dt = 0.001
burn_in_steps = 10_000
spike_threshold = 0.0
min_spikes per realization = 5
noise amplitudes = [0.02, 0.03, 0.05, 0.075, 0.10, 0.20, 0.40, 0.70]
seeds = [11, 23, 37, 41, 53, 67, 79, 97]
uncertainty weight = 2.0
accepted coordinate interval = [0.03, 0.20]
```

The only independent variable introduced by this control is total trajectory length.

## Preregistered horizons

Evaluate exactly three total-step budgets:

```text
40_000 steps
80_000 steps
160_000 steps
```

The parent protocol's `80_000`-step run is the reference horizon. All three horizons retain the same `10_000` burn-in steps, so the analyzed post-burn-in windows are 30_000, 70_000, and 150_000 steps respectively.

No additional horizon may be introduced after reading outcomes from these three runs.

## Metrics and decision rule

For every horizon and every noise amplitude:

1. run all eight declared seeds;
2. reject that horizon's calibration if any replicate has fewer than `min_spikes` after burn-in;
3. compute each replicate's `CV_ISI` exactly as in the parent protocol;
4. compute mean CV, sample standard deviation, standard error, and the parent's `w*SE` conservative edge-separation rule;
5. classify the horizon as one of:
   - `INTERIOR_MINIMUM_ACCEPTED`;
   - `NO_SEPARATED_INTERIOR_MINIMUM`;
   - `INSUFFICIENT_SPIKES`;
   - `PROTOCOL_FAILURE`.

The robustness control passes only if all three horizons yield the same decision class and, when that class is `INTERIOR_MINIMUM_ACCEPTED`, every selected coordinate lies inside `[0.03, 0.20]`.

Exact equality of selected coordinates across horizons is not required. The control concerns the preregistered decision class and broad accepted region, not a claim that a discrete optimum is invariant to sample duration.

## Interpretation boundary

A passing result supports only the statement that the parent calibration's qualitative decision is stable across these three finite observation horizons for the declared numerical system.

It does not establish:

- asymptotic convergence;
- ergodicity;
- independence from initial conditions in other FHN regimes;
- a universal observation duration;
- transfer to inverse stochastic resonance, networks, attention systems, or physical hardware;
- a new resonance mechanism.

## Motivation and external baseline

Yamakou, Krüger, and Schulz-Baldes (2026), *Analysis of inverse stochastic resonance: Effects of neural excitability and timescale separation*, arXiv:2608.03454, distinguish finite-time Monte Carlo behavior from asymptotic behavior in a bistable FitzHugh–Nagumo regime. That work studies a different phenomenon and parameter regime; it is used here only to motivate an explicit finite-horizon robustness control, not as evidence for NoiseLab's coherence-resonance outcome.

Primary source: https://arxiv.org/abs/2608.03454

## Provenance rule

This document must be committed before any execution performed specifically for this horizon comparison. Results, including negative or inconclusive outcomes, must be recorded without editing these frozen horizons, seeds, grid coordinates, burn-in, decision classes, or acceptance interval.
