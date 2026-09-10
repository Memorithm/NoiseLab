# FitzHugh–Nagumo coherence-resonance calibration

## Status

Calibration experiment. No novelty claim.

This experiment is a truth/control step in NoiseLab's progression from systems with established resonance behavior toward systems whose useful perturbation regimes are unknown.

## Model

Use the canonical excitable FitzHugh–Nagumo form

```text
epsilon * dx/dt = x - x^3/3 - y
dy/dt = x + a + sigma * xi(t)
```

with `epsilon > 0`, `a > 1`, and Gaussian white noise applied to the slow variable. The implementation uses the equivalent stochastic increment `sigma * sqrt(dt) * N(0,1)` in a fixed-step Euler–Maruyama update of `y`; deterministic Gaussian draws come from SciRust `scirust_sim::SplitMix64`.

The deterministic equilibrium is

```text
x* = -a
y* = a^3/3 - a
```

and is used as the initial state for every realization.

The local Euler/Euler–Maruyama stepping code is experiment-specific. NoiseLab does not present it as a general SDE solver or as a replacement for a future validated SciRust SDE primitive.

## Research question

Does this declared excitable system exhibit a sampled intermediate noise amplitude at which noise-induced spikes are more temporally regular than in both lower-noise and higher-noise regimes?

The target quantity is the coefficient of variation of inter-spike intervals,

```text
CV_ISI = sample_stddev(ISI) / mean(ISI)
```

where a spike is defined by an upward crossing of the declared activator threshold. Crossing times are linearly interpolated between samples.

Lower `CV_ISI` means more regular spike timing.

## Predeclared calibration configuration

The executable regression uses:

```text
epsilon = 0.01
a = 1.05
dt = 0.001
steps = 80_000
burn_in_steps = 10_000
spike_threshold = 0.0
min_spikes per realization = 5
noise amplitudes = [0.02, 0.03, 0.05, 0.075, 0.10, 0.20, 0.40, 0.70]
seeds = [11, 23, 37, 41, 53, 67, 79, 97]
uncertainty weight = 2.0
```

These values are calibration parameters, not claimed physical constants or universal optima.

## Decision rule

For every noise amplitude:

1. run every declared seed;
2. reject the calibration if any realization has fewer than `min_spikes` after burn-in rather than silently dropping that seed;
3. compute each realization's `CV_ISI`;
4. compute the mean CV and its sample standard deviation across seeds.

The best sampled interior point is retained only if:

- its mean CV is strictly below both immediate neighboring sampled means; and
- its uncertainty-penalized upper value `mean + w*SE` lies below both scan-edge uncertainty-penalized lower values `mean - w*SE`.

This `w*SE` construction is a conservative ranking heuristic. It is not described as a formal confidence interval, p-value, or proof.

The regression accepts the selected coordinate only if it lies in the broad predeclared interval `[0.03, 0.20]`. The purpose of this interval is to detect gross implementation/protocol regressions without encoding a single expected optimum into the algorithm.

## Null / negative outcomes

The calibration is not promoted if:

- there is no strict sampled interior CV minimum;
- the candidate minimum is not separated from both scan edges by the declared uncertainty heuristic;
- the low-noise or high-noise edge is at least as regular under the declared metric;
- an input, state, spike train, or noise sweep is malformed;
- a replicate does not provide enough spikes for the declared CV measurement.

A negative result is a valid NoiseLab result and must not be repaired by moving the sweep, changing seeds, weakening the uncertainty rule, or relaxing the metric after observing CI output.

## Interpretation boundary

A passing calibration supports the narrow statement that, for this declared numerical system and protocol, the sampled inter-spike regularity has a reproducible intermediate optimum consistent with the established coherence-resonance phenomenon.

It does **not** establish:

- a universal noise optimum;
- a universal "perfect moment" for arbitrary systems;
- a new noise family;
- an analytic law relating the optimum to arbitrary system parameters;
- applicability to attention, RoPE, FLAT-ATTENTION, physical hardware, or other targets without independent experiments.

## References

- A. S. Pikovsky and J. Kurths, "Coherence Resonance in a Noise-Driven Excitable System," *Physical Review Letters* 78, 775–778 (1997), DOI: 10.1103/PhysRevLett.78.775.
- B. Lindner, J. García-Ojalvo, A. Neiman, and L. Schimansky-Geier, "Effects of noise in excitable systems," *Physics Reports* 392, 321–424 (2004), DOI: 10.1016/j.physrep.2003.10.015.

Primary/reference literature motivates the calibration target. Executed NoiseLab evidence remains authoritative for this repository's specific implementation and parameterization.
