# SciRust reuse contract

NoiseLab is a research bench, not a second scientific-computing platform. Generic numerical infrastructure should be reused from SciRust whenever the upstream contract is suitable.

## Pinned upstream state

Initial audited revision:

`Memorithm/scirust@f57d598bf03e5dfb16ec6423e4e43105a77540d3`

The pin is deliberate. Research artifacts must record the exact SciRust revision used to generate evidence. Upgrading the pin is a scientific change when it can alter generated samples, spectra, trajectories or metrics.

## Reused now

### Deterministic stochastic source

`scirust-sim::SplitMix64`

NoiseLab uses the upstream explicit-seed stream and its Gaussian transform. The upstream implementation is validated against the published SplitMix64 output sequence and keeps ambient randomness out of the simulation path.

### Correlated stochastic control

`scirust-sim::stochastic::ou_path`

The Ornstein-Uhlenbeck path uses SciRust's exact Gaussian transition law rather than introducing an Euler-Maruyama approximation in NoiseLab. This gives NoiseLab a correlated stochastic control with an analytic stationary variance oracle.

### Spectral analysis

`scirust-signal`

The initial bridge reuses Hann windows, Welch PSD, PSD centroid, PSD spread and PSD flatness. Welch metadata (`segments`, `dropped_samples`) is retained in NoiseLab results so spectral estimates cannot silently hide how much data contributed.

## Candidate reuse next

Before implementing any equivalent locally, audit these SciRust capabilities:

- `scirust-signal`: autocorrelation, entropy, kurtosis, skewness, filtering, STFT/spectrogram and noise characterization;
- `scirust-sim`: `System`, deterministic RK4, adaptive Dormand-Prince, symplectic second-order integration and ready-made nonlinear systems;
- `scirust-stats`: probability distributions and statistical tests needed for matched-control experiments;
- SciRust reproducibility/provenance substrate where it can be depended on without importing unrelated application semantics;
- tensor and attention primitives only through narrow adapters when NoiseLab reaches tensor perturbation experiments.

## Upstream rule

If NoiseLab needs a generic mathematical primitive that would also benefit SciRust, implement or improve it in SciRust first, validate it there, then consume the pinned revision from NoiseLab. NoiseLab should own only noise-specific ontology, perturbation composition, experiment protocols, response signatures, novelty tests and resonance-search semantics.

## Licensing

NoiseLab and SciRust use PolyForm Noncommercial 1.0.0. Upstream source remains authored and versioned in SciRust; NoiseLab consumes it as a dependency rather than vendoring copies, preserving provenance and avoiding divergence.
