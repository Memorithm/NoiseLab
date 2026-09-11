# Noise as Information

NoiseLab treats "noise" as an observational classification, not as a guarantee that a fluctuation is useless.

The working hypothesis for this research axis is deliberately narrower than the claim that all noise is informative:

> A measurable fraction of the apparent noise in some systems can carry information about unmodelled, hidden or unresolved degrees of freedom.

This hypothesis is falsifiable and system-dependent. A negative result is a valid result.

## Why this belongs in NoiseLab

A common processing pipeline starts from

```text
x(t) = s(t) + n(t)
```

and attempts to suppress `n(t)`. That decomposition can be useful, but in a real system the residual can instead be a projection of omitted state:

```text
x(t) = F(z_1(t), z_2(t), ..., z_k(t))
```

If a declared noise-like component is statistically dependent on a hidden or system state `Z`, removing it before characterization can destroy evidence about `Z` even while conventional signal-to-noise metrics improve.

NoiseLab therefore adopts the rule:

> Characterize before filtering.

Filtering is an intervention whose information loss or gain must be measured, not assumed.

## Research categories

Experiments should distinguish at least these roles:

1. **Nuisance noise** — fluctuations that reduce a declared estimation or control objective and show no reproducible relation to the variables under study.
2. **Informative noise** — fluctuations statistically dependent on a declared hidden/system variable.
3. **Functional noise** — perturbations whose presence improves a declared system objective in a bounded operating region.
4. **Structural apparent noise** — unresolved deterministic or chaotic dynamics that appear stochastic at the chosen observation scale.
5. **Measurement noise** — fluctuations introduced by the measurement chain, which may carry information about the instrument rather than the target system.

These categories are not mutually exclusive and must not be assigned from appearance alone.

## Stage 0 executable diagnostic

`src/information.rs` provides a deliberately small first diagnostic:

- empirical entropy of a discrete hidden/system state;
- equal-width histogram mutual-information estimate between a real-valued noise-like component and that state;
- before/after audit for denoising, smoothing, thresholding or another declared transformation;
- signed information change, so an apparent information gain is not silently clipped away;
- known-answer controls where one bit is present, absent, retained or destroyed.

The current estimator is a histogram plug-in estimator. It is not an exact continuous mutual-information estimator and is sensitive to sample count, bin count and range. Experiments using it must preregister the binning rule and include bin-count sensitivity checks.

## Required experimental protocol

For a candidate noise-like component `N`, hidden/system state `Z` and transformation `T`, record at minimum:

```text
H(Z)
I(N; Z)
I(T(N); Z)
delta_I = I(T(N); Z) - I(N; Z)
```

alongside the conventional task objective, raw data and transformation parameters.

A denoiser is not considered information-preserving merely because it improves SNR or visual smoothness. Conversely, a positive mutual-information estimate is not sufficient to establish mechanism or causality.

Every confirmatory experiment should include:

- a clean/no-perturbation or identity control where applicable;
- a negative independence control;
- repeated seeds or independent trajectories where stochastic generation is involved;
- uncertainty or resampling analysis appropriate to the estimator;
- a preregistered observation window and binning/estimation rule;
- preservation of raw measurements before any filtering;
- explicit separation of exploratory and confirmatory analysis.

## Next research steps

The Stage 0 histogram diagnostic is only a calibration layer. Subsequent work should compare stronger estimators and system classes without replacing this known-answer baseline:

1. lagged and conditional information for dynamical systems;
2. permutation/surrogate baselines to detect finite-sample estimator bias;
3. information retention across frequency-selective filters;
4. hidden-state inference from oscillator, bistable, FitzHugh-Nagumo and laser residuals;
5. attention/RoPE/FLAT and KV-state experiments where perturbations may expose otherwise hidden internal state;
6. transfer of reusable, domain-independent information-theory primitives to SciRust once their API and numerical behavior are qualified.

No result from this axis should be described as proving that "noise is information" in general. The admissible claim is always tied to the declared system, state variable, estimator and protocol.
