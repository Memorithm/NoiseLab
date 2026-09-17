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

## Stage 0.1 finite-sample permutation null

`permutation_mutual_information_null` adds a deterministic association null without replacing the Stage 0 known-answer baseline. For a fixed observation vector and discrete hidden-state vector it:

- quantizes the observation once under the declared equal-width bin count;
- preserves the observation values and the exact hidden-state label multiset;
- independently permutes the hidden-state labels for every surrogate with SciRust `SplitMix64` and unbiased Fisher-Yates indices;
- retains every surrogate mutual-information score rather than only a summary;
- reports the surrogate mean, the number of surrogate scores at or above the observed score, and the finite one-sided p-value `(1 + exceedances) / (1 + permutations)`;
- rejects zero permutations and requests above the explicit safety ceiling.

`bins`, `permutations` and `seed` are protocol inputs and must be preregistered. A small permutation p-value rejects only the declared exchangeable-label null under this histogram estimator. It does **not** establish causality, mechanism, a universal information property, or an information-preserving denoiser. Temporal data with autocorrelation generally require a structure-preserving null rather than unrestricted label permutation; this Stage 0.1 primitive must not be used to erase that distinction.

## Stage 0.2 structure-preserving cyclic-shift null

`cyclic_shift_mutual_information_null` supplies the first bounded time-series surrogate for cases where unrestricted label exchangeability is not defensible. The observation series stays fixed while each preregistered non-zero offset circularly rotates the **entire hidden-state sequence**, preserving its label marginal and circular lag organization. Duplicate, zero and out-of-range offsets fail closed; the exact offset list and every surrogate MI score are retained in declaration order.

The returned `(1 + exceedances) / (1 + shifts)` quantity is named a **corrected tail fraction**, not an unconditional p-value. Circular wrapping and the chosen offsets must be scientifically justified and frozen before outcomes are inspected. Periodic structure can legitimately preserve high MI under some shifts—for example a half-period binary phase inversion—and that remains null evidence rather than being forced toward independence. This primitive does not establish causality or solve nonstationarity, boundary effects, lag selection or conditional information.

## Stage 0.3 spectrum-matched observation null

`spectral_surrogate_mutual_information_null` adds a second structure-preserving control for power-of-two real-valued observation records. It holds the hidden-state sequence fixed and uses the repository's pinned SciRust `phase_randomized_surrogate` through `spectral_phase_null`, preserving Fourier magnitudes (subject to the upstream floating-point contract) while randomizing phase relationships. Every derived surrogate seed and every surrogate MI score is retained.

The observed series defines the equal-width histogram range once; every phase-randomized surrogate is quantized against those same bin edges. Values outside the observed range are clamped to the edge bins instead of moving the bins draw-by-draw. `bins`, surrogate count and root seed are protocol inputs and must be frozen before outcomes are inspected. The finite `(1 + exceedances) / (1 + surrogates)` value is reported only as a corrected tail fraction. It is not automatically an exact p-value and does not establish causality, a common physical mechanism, universality, or denoising benefit. This control specifically asks whether the observed state association survives a null that preserves the observation's second-order spectral magnitude structure while disrupting phase alignment.

## Stage 0.4 non-circular lagged association diagnostic

`lagged_histogram_mutual_information_bits` evaluates a caller-declared lag set without circular wrapping. The complete observation series is quantized once, so bin edges remain fixed across lags; each lag then uses only its overlapping observation/state pairs. By convention a positive lag pairs `N[t]` with `Z[t + lag]`, a negative lag pairs `N[t - lag]` with `Z[t]`, and lag zero is contemporaneous association. Duplicate lags and lags leaving no overlap fail closed, and declaration order plus aligned sample count are retained.

`conditional_histogram_mutual_information_bits` adds the first bounded conditional-information calibration primitive. It estimates empirical `I(N; Z | C)` with the same frozen equal-width quantization of the real-valued observation and exact discrete labels for the declared state and conditioner. Known-answer controls include a state-coded signal whose association vanishes when conditioning on the state itself and an XOR construction whose marginal MI is zero while conditional MI is one bit. A positive residual estimate means only that association remains under the declared `C` and finite histogram estimator; it is not causal evidence, does not prove `C` sufficient, and does not identify a physical mechanism. Binning and the conditioning variable must be frozen before outcome inspection.

This diagnostic is descriptive and does not infer causal direction. Searching a large lag range after inspecting outcomes creates a selection problem; confirmatory use must preregister the lag set or declare a separate multiplicity/selection procedure. Different lags also use different aligned sample counts, which must remain visible rather than being silently treated as equal evidence.

## Required experimental protocol

For a candidate noise-like component `N`, hidden/system state `Z` and transformation `T`, record at minimum:

```text
H(Z)
I(N; Z)
I(T(N); Z)
delta_I = I(T(N); Z) - I(N; Z)
```

alongside the conventional task objective, raw data and transformation parameters.

When the unrestricted permutation null is scientifically admissible, also retain its `bins`, `permutations`, `seed`, complete surrogate score sequence, exceedance count and corrected p-value. For a cyclic-shift null, retain the preregistered offset list, complete surrogate scores, exceedance count and corrected tail fraction, and justify circular wrapping. If neither exchangeability nor circular shifting is defensible, declare another structure-preserving null instead of silently applying either primitive.

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

The Stage 0 histogram diagnostic and Stage 0.1 unrestricted permutation null are calibration layers. Subsequent work should compare stronger estimators and system classes without replacing these known-answer baselines:

1. conditional-information experiments on preregistered dynamical-system state/condition variables, building on the implemented calibration primitive and non-circular lagged association diagnostic;
2. richer structure-preserving surrogate families beyond the Stage 0.2 cyclic-shift and Stage 0.3 spectrum-matched controls, especially for nonstationary and noncircular records;
3. information retention across frequency-selective filters;
4. hidden-state inference from oscillator, bistable, FitzHugh-Nagumo and laser residuals;
5. attention/RoPE/FLAT and KV-state experiments where perturbations may expose otherwise hidden internal state;
6. transfer of reusable, domain-independent information-theory primitives to SciRust once their API and numerical behavior are qualified.

No result from this axis should be described as proving that "noise is information" in general. The admissible claim is always tied to the declared system, state variable, estimator and protocol.
