# NoiseLab

NoiseLab is Memorithm's Rust research bench for the controlled generation, characterization, injection and experimental study of noise and stochastic perturbations.

The project has four complementary goals:

1. build a reproducible taxonomy and implementation of known noise/process families without conflating physical noise, stochastic processes, random fields, numerical errors and procedural textures;
2. search for falsifiable new perturbation families and for system-dependent operating points where a perturbation improves a declared objective rather than merely degrading it;
3. test whether an apparently noisy component carries measurable information about hidden, unresolved or unmodelled system state, and quantify how transformations such as denoising preserve or destroy that information;
4. test whether fluctuation series from different mechanisms approach common statistical representations under controlled coarse graining, with surrogate nulls that separate non-trivial temporal/scaling structure from generic averaging effects.

A measured optimum is evidence for one declared system, objective and protocol. It is not automatically a universal law, a new noise family or a proof. Likewise, measured dependence between a noise-like component and system state is evidence for that declared experiment, not a claim that all noise is informative. Multiscale convergence is treated only as a candidate universality signal until it survives declared null models, independent source families and confirmatory validation.

## Scientific direction

NoiseLab will study perturbations over targets including dynamical systems, numerical algorithms, attention mechanisms, RoPE, FLAT-ATTENTION and KV-state experiments. Mature mathematical conjectures may be transferred to ProofLab; information/attention experiments may be transferred to TDI; reusable mathematical primitives belong upstream in SciRust rather than being duplicated here.

The bench now separates four related research views:

- **Noise as Perturbation** — characterize and inject reproducible stochastic processes;
- **Noise as Functional Mechanism** — search for bounded operating regions where perturbations improve a declared objective;
- **Noise as Information** — test whether fluctuations contain information about declared hidden/system variables and measure what filtering removes;
- **Fluctuation Universality** — test whether independent mechanisms converge toward shared multiscale statistical representations beyond matched surrogate controls.

For the information axis, NoiseLab adopts a simple rule: **characterize before filtering**. Raw observations must be preserved, and filtering is treated as an intervention whose information loss or gain is measured rather than assumed. See `docs/research/noise-as-information.md`.

For fluctuation universality, NoiseLab uses controlled block coarse graining, fixed multiscale descriptors, scaling exponents and independently shuffled temporal-null surrogates. A passing surrogate gate is labeled only as a candidate shared scaling structure, never as proof of a common physical or cosmological origin. See `docs/research/fluctuation-universality.md`.

A central research line is **resonant operating-point discovery**. Given a system state `x`, perturbation family `N(theta)`, intervention location `I`, control parameters `u` and a declared utility `J`, NoiseLab will estimate bounded optima of the form

```text
(theta*, u*, t*) = argmax E[J(system | x, N(theta), I, u, t)]
```

subject to explicit safety/resource constraints and comparison against a no-perturbation control. There is no assumption that a universal "perfect moment" exists.

## SciRust foundation

NoiseLab intentionally reuses pinned SciRust primitives instead of creating competing implementations.

Current pin:

```text
Memorithm/scirust@f57d598bf03e5dfb16ec6423e4e43105a77540d3
```

The first bridge uses:

- `scirust-sim::SplitMix64` for explicit-seed deterministic pseudo-random streams and universality-test surrogates;
- `scirust-sim::stochastic::ou_path` for an exact-transition correlated stochastic control process;
- `scirust-signal` for Hann windows, Welch PSD, spectral centroid, spread and flatness.

See `docs/SCIRUST_REUSE.md`.

## Build

```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo fmt --all -- --check
```

## Licensing

NoiseLab is licensed under the PolyForm Noncommercial License 1.0.0. See `LICENSE.md`.

Pinned dependencies retain their own license terms; depending on an upstream Memorithm crate does not relicense that upstream code as part of NoiseLab.
