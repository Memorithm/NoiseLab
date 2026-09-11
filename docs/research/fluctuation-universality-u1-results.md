# Fluctuation universality — Stage U1 observed results

## Evidence identity

This document records the first observed Stage U1 panel after the protocol and executable panel had already been frozen and merged.

Authoritative NoiseLab revision:

```text
2c0ab4aedda07ec602949eb3d6ec6f8b3edb27d2
```

Authoritative post-merge GitHub Actions run:

```text
Fluctuation U1 Report run 34592385867
```

The run executed `StageU1PanelConfig::default()` on `main` and completed successfully. The general NoiseLab CI on the same revision also completed Format, Check, Clippy and Test successfully.

The values below were not used to alter the Stage U1 endpoint, source parameters, scale grid, surrogate count or evidence labels.

## Results

| Pair | Purpose | Fine distance | Terminal distance | Convergence | Scaling exponent gap | Empirical p | Evidence |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| GaussianWhiteA / GaussianWhiteB | WithinIidControl | 0.040912028339 | 0.082719417414 | -0.041807389075 | 0.010599751841 | 0.105 | NoObservedConvergence |
| OrnsteinUhlenbeckFast / OrnsteinUhlenbeckSlow | WithinOuFamily | 0.198860757298 | 0.733764684576 | -0.534903927278 | 0.146223250376 | 1.000 | NoObservedConvergence |
| GaussianWhiteA / OrnsteinUhlenbeckSlow | DistinctTemporalControl | 1.090557449894 | 2.840212037136 | -1.749654587242 | 0.967194650354 | 1.000 | NoObservedConvergence |
| OrnsteinUhlenbeckSlow / DoubleWellLangevin | CrossMechanismExploratory | 0.345334111361 | 0.362798866391 | -0.017464755030 | 0.020581321097 | 0.930 | NoObservedConvergence |

The exact same rows were observed in the preceding PR run and in the post-merge `main` run, confirming deterministic reproduction under the frozen seeds.

## Interpretation

### Same-family IID control

The two independent Gaussian white-noise realizations were already close at the finest scale. Their descriptor distance increased slightly rather than shrinking under the tested coarse-graining range. This is not evidence against Gaussian finite-sample equivalence; it shows why Stage U1 retained absolute fine and terminal distances instead of interpreting a convergence score alone as class membership.

### OU correlation-time change

Fast and slow Ornstein-Uhlenbeck processes became substantially farther apart in the selected descriptor as block size increased. Under this scale grid, correlation time is not erased quickly enough to produce a shared coarse representation according to the Stage U0 distance.

This is useful calibration: the pipeline does not automatically collapse processes merely because they share a named stochastic mechanism.

### IID versus correlated control

Gaussian white noise and slow OU noise diverged strongly with coarse graining. This is the desired qualitative behavior of a negative temporal-structure control: the current representation remains sensitive to temporal organization rather than reducing all finite-variance stochastic processes to a superficial common class over the tested scales.

### OU versus bistable Langevin

The exploratory cross-mechanism pair is the most relevant to the broad conjecture in this first panel. Its fine and terminal distances were similar, but the terminal distance was slightly larger:

```text
0.345334111361 -> 0.362798866391
```

The resulting convergence score was negative and the shuffled-null empirical p-value was `0.93`. Stage U1 therefore provides **no observed convergence** for this pair.

The variance-scaling-exponent gap was small (`0.020581321097`), which is itself an important guardrail: a similar single scaling exponent is not sufficient evidence of a common universality class when the preregistered multivariate descriptor does not converge.

## Stage U1 conclusion

Stage U1 does **not** support the fluctuation-universality conjecture under the current source panel, scale range, descriptor and shuffled-surrogate null.

This is a negative result, not a reason to rewrite the endpoint. It rules out one simple version of the conjecture:

> the tested Gaussian, OU and bistable-Langevin fluctuations do not spontaneously converge toward the same selected multiscale representation over block sizes 1 through 16.

It does **not** prove that no fluctuation universality exists. In particular, Stage U1 does not yet test:

- other observables from the same dynamical systems;
- substantially wider scale ranges;
- critical or near-critical regimes;
- spectral/phase-preserving null models;
- information-theoretic invariants;
- oscillator, semiconductor-laser, FitzHugh-Nagumo, attention or real-world physical datasets;
- genuine renormalization transformations adapted to each system rather than arithmetic block means alone.

## Next protocol requirements

Any Stage U2 extension must be specified before its outcomes are inspected. The next protocol should add at least:

1. phase-randomized or otherwise spectrum-matched surrogates so spectral similarity cannot masquerade as deeper universality;
2. explicit observables from additional existing NoiseLab mechanisms, especially driven oscillator, semiconductor laser and FitzHugh-Nagumo;
3. a wider but preregistered scale grid where sample support remains adequate;
4. an information-theoretic descriptor linked to the existing `Noise as Information` axis;
5. multiplicity handling for any confirmatory family of pairwise tests;
6. at least one protected source family not used to design the Stage U2 representation.

The Stage U1 values above remain immutable evidence and must not be discarded if a later representation produces a positive result.
