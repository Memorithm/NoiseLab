# Fluctuation universality — Stage U2 preregistration

Status: **preregistered; SciRust spectral-surrogate dependency pinned; NoiseLab consumption qualified against the pinned revision; outcome-blind readiness gate merged; U2 not executed by this update**.

This document freezes the Stage U2 cross-mechanism protocol before any U2 outcomes are inspected. It extends the already-recorded negative/diagnostic U1 work. It does not revise U0/U1 results and does not assert universality.

## Question

After the U0/U1 calibration program, do fluctuation observables from mechanistically distinct NoiseLab systems approach a common multiscale representation more strongly than expected from controls that preserve one-point statistics and, separately, Fourier magnitudes?

## Hypotheses

`H0` — any apparent cross-mechanism convergence is explained by empirical marginals, generic averaging, and/or second-order spectral structure captured by the declared surrogate controls.

`H1` — at least one preregistered cross-mechanism pair exhibits positive multiscale convergence that separates from both the shuffled-null distribution and the phase-randomized spectral-null distribution under the fixed decision rule below.

Rejecting `H0` at U2 yields only `CrossMechanismCandidate`. It is not evidence of a universal law, a common physical source, cosmological ancestry, or a new kind of noise.

## Frozen source families

U2 uses only existing NoiseLab mechanisms whose scalar observable can be generated with preserved provenance:

1. driven damped oscillator residual;
2. semiconductor-laser relaxation residual;
3. bistable Langevin trajectory residual;
4. FitzHugh–Nagumo fluctuation residual.

Attention/RoPE perturbation residuals remain excluded from this first U2 panel. They may enter a later preregistration only after a scalar observable and baseline are independently frozen.

Each source series must record the producing NoiseLab commit, all model parameters, seed, sample cadence, burn-in/window rule and observable-extraction rule. No source may be selected or discarded because of its U2 score.

## Fixed representation

Reuse the U0 multiscale descriptor and scales without retuning:

- non-overlapping block means at block sizes `1, 2, 4, 8, 16`;
- variance after raw-series standardization;
- lag-1 Pearson autocorrelation;
- population excess kurtosis;
- mean-centered sign-change rate;
- pairwise descriptor distance exactly as defined by the U0 protocol;
- primary convergence score `distance(scale_first) - distance(scale_last)`.

Any future descriptor change requires a new preregistration and cannot be back-applied to this U2 panel.

## Controls

Every observed cross-mechanism pair is compared against two independently seeded null families.

### Null A — shuffled empirical marginal

Use the existing Fisher–Yates/`SplitMix64` shuffled control from U0. It preserves each empirical one-point distribution and destroys temporal ordering.

### Null B — phase-randomized spectral surrogate

Use the reusable SciRust spectral-surrogate primitive only from the immutable revision pinned by NoiseLab. The required contract is:

- preserve DC and Nyquist bins when present;
- preserve Fourier magnitudes up to transform roundoff;
- randomize only nontrivial conjugate-pair phases from an explicit seed;
- preserve real-valued reconstruction;
- reject non-finite or otherwise unsupported input rather than silently modifying it.

The preregistration-time candidate PR #1421 is no longer a permitted dependency reference. The qualified dependency is pinned in `docs/research/dependencies/scirust-spectral-surrogate-u2.md` to SciRust `master` commit `0e2eaccac631b689f97c242c47bad11d433847d9`, which contains `scirust-signal/src/surrogate.rs` after PR #1422 merged. NoiseLab PR #26 qualified the local spectral-null adapter against that exact revision, and PR #27 added the outcome-blind `U2Readiness` gate that fixes the revision, four source families, six unordered pairs, and at least 199 surrogates per null family before execution. If the dependency revision changes, a new explicit pin and requalification are required before execution. Neither adapter qualification nor readiness validation is U2 outcome evidence.

The phase-randomized null tests whether an apparent U2 signal exceeds what can be explained by the observed power spectrum/second-order Fourier structure. It does not preserve all nonlinear temporal structure.

## Repetitions and seeds

For each cross-mechanism pair and each null family:

- use at least `199` surrogate repetitions;
- use explicit recorded seeds;
- use the standard one-sided `+1` empirical p-value correction;
- keep the same observed source series for both null families;
- do not regenerate a source series because a surrogate result is inconvenient.

Higher repetition counts may be preregistered in a later confirmatory stage but cannot replace this panel after outcomes are observed.

## Primary decision rule

For each cross-mechanism pair, compute the observed convergence score and the two one-sided empirical p-values:

- `p_shuffle` against Null A;
- `p_phase` against Null B.

With `alpha = 0.05`, classify:

- `NoObservedConvergence` if observed score `<= 0`;
- `CompatibleWithMarginalNull` if score `> 0` and `p_shuffle > alpha`;
- `SpectrumExplainedCandidate` if score `> 0`, `p_shuffle <= alpha`, and `p_phase > alpha`;
- `CrossMechanismCandidate` only if score `> 0`, `p_shuffle <= alpha`, and `p_phase <= alpha`.

No multiple-comparison correction is silently invented in this exploratory U2 panel. Therefore `CrossMechanismCandidate` remains exploratory. A confirmatory family-level claim requires a separate preregistration with an explicit multiplicity policy and untouched source families.

## Mandatory negative-result preservation

The evidence record must retain:

- non-positive convergence;
- shuffled-null compatibility;
- phase-null compatibility;
- numerical/protocol rejection;
- pairs for which a source observable cannot be generated under the frozen extraction rule;
- results that contradict the broader universality conjecture.

No failed pair may be removed from the reported U2 matrix.

## Anti-leakage and provenance

- Do not inspect U2 outcomes while choosing source parameters, observable rules, descriptor scales, surrogate counts or decision thresholds.
- Do not use final/confirmatory holdouts from unrelated TDI programs.
- Keep raw series, derived descriptors, surrogate seeds, source revisions and SciRust dependency revision addressable in the result artifact.
- Treat SciRust as the provider of the generic spectral-surrogate primitive; do not duplicate FFT/RNG logic in NoiseLab.
- If the SciRust primitive changes after the pinned revision, requalification is required and the exact consumed revision must remain recorded.
- Any U2 executor must call the outcome-blind readiness gate before generating the first U2 pair result; this documentation update itself does not execute U2.

## Stop condition

Stage U2 is complete when all six unordered pairs among the four frozen source families have either:

1. a reproducible result under both null families, or
2. a recorded protocol/implementation failure that prevents a valid comparison.

The stage must stop without retuning after that matrix is complete. Any mechanism investigation motivated by a candidate belongs to a separately preregistered follow-up.
