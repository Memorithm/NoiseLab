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

## Executable research reports

The U2 universality panel and FHN observation-horizon panel have executable
reporters. They default to **non-scientific smoke**, reject conflicting FULL/SMOKE
flags and emit explicit workload metadata before their TSV tables.

The `Research report smoke contracts` workflow builds both reporters at an exact
source revision, checks invalid flag rejection, runs each smoke workload twice,
requires identical reports and retains diagnostics plus the resolved dependency
graph for 14 days. This checks execution and report contracts, not scientific
universality or FHN robustness. See [report automation](docs/REPORT_AUTOMATION.md).

U2 can additionally retain the exact post-burn-in residual arrays before pair
analysis through `NOISELAB_U2_CAPTURE_DIR`: lossless binary64 inputs, extraction
provenance, descriptors, the preregistered surrogate job manifest, every realized
surrogate array and every retained surrogate score. Array files bind exact
pair/null/repetition/seed identities and are persisted before pair analysis;
score export revalidates those frozen job identities after analysis. The CI
repeats and compares these snapshots, checks no-overwrite behavior and seals the
combined U2/FHN artifact directory with a SHA-256 manifest. Integrity is not
authentication or a scientific verdict; full pre-extraction trajectories remain
outside the capture bundle. Archived surrogate arrays can also be replayed
directly without regenerating source trajectories or null transformations:
`verify_archived_surrogate_scores` validates the frozen job/index binding,
decodes the retained binary64 arrays, recomputes the declared multiscale score
and requires exact IEEE-754 bit equality with each archived score row.
`replay_u2_statistics_from_scores` then provides a separate core primitive that
revalidates the exact preregistered score ordering and independently recomputes
the one-sided `+1` p-values and frozen U2 decision from those retained scores and
an observed convergence score, without regenerating surrogate arrays. Neither
primitive binds the replay to a stored pair-result artifact or authorizes a
scientific claim; that end-to-end evidence boundary remains separate. See
[input capture and bundles](docs/research/u2-input-capture.md).

```bash
cargo build --examples
python3 scripts/check_smoke_reports.py --output-dir /tmp/noiselab-smoke-new-run
```

Use a new output directory for each run. Existing output files are not overwritten.
Full scientific loads remain subject to their frozen protocols and separate
result review; smoke artifacts must not be promoted to scientific evidence.

## SciRust foundation

NoiseLab intentionally reuses pinned SciRust primitives instead of creating competing implementations.

Current pin used by both `scirust-sim` and `scirust-signal` in `Cargo.toml`:

```text
Memorithm/scirust@0e2eaccac631b689f97c242c47bad11d433847d9
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
