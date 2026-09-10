# Resonant operating-point research

## Question

NoiseLab studies whether a system has a state-, time- and perturbation-dependent operating point at which a controlled perturbation improves a declared objective.

The general experimental object is

```text
z = (t_onset, amplitude, frequency, phase, ...)
G(z) = E[J(perturbed system at z) - J(matched control)] - cost(z)
z* = argmax G(z), subject to declared constraints
```

There is no assumption that one universal "perfect moment" exists for arbitrary systems. An optimum is always conditional on the target system, current state, intervention, utility, horizon, constraints and experimental distribution.

## Evidence levels

NoiseLab keeps these statements separate:

1. `BEST_SAMPLED_POINT`: one candidate maximizes the declared score over the finite candidate set.
2. `INTERIOR_RESPONSE_PEAK`: a scalar sweep contains a sampled local maximum above both immediate neighbors by a declared minimum prominence.
3. `MECHANISM_COMPATIBLE`: the measured peak is consistent with a preregistered physical/mathematical mechanism and its controls.
4. `CONFIRMED`: an isolated confirmatory protocol reproduces the preregistered effect.
5. `PROVED`: reserved for a mathematical statement accepted through an appropriate formal proof workflow; numerical search cannot assign this status.

## Paired search

`search_operating_points` evaluates the unperturbed control once per seed and reuses the same seed for every candidate. This paired/common-random-numbers design reduces variance from stochastic background differences without claiming to remove model uncertainty.

Candidates are ranked with

```text
conservative_score = mean_uplift - w * standard_error
```

where `w` is explicitly configured. This is a ranking penalty, not automatically a confidence interval: a formal interval requires its own distributional or resampling assumptions.

The initial candidate coordinates are onset time, amplitude, frequency and phase. Later work may add duration, correlation time, spectral exponent, injection site and structured/algebraic perturbation parameters.

## Analytic controls

### Linear damped oscillator

For

```text
x'' + 2 zeta omega_n x' + omega_n^2 x = (F/m) cos(omega t)
```

the steady-state displacement-amplitude maximum occurs at

```text
omega_r = omega_n * sqrt(1 - 2 zeta^2)
```

when `zeta < 1/sqrt(2)`. For greater damping there is no non-zero-frequency displacement-amplitude peak of that form.

This is a deliberately narrow control. It is not a general formula for nonlinear or stochastic systems.

Reference: MIT OpenCourseWare, Differential Equations, resonance material; see also standard damped-driven oscillator derivations.

- https://ocw.mit.edu/courses/18-03sc-differential-equations-fall-2011/

### Classical stochastic-resonance rate matching

For a symmetric bistable system in the classical weak-noise/adiabatic picture, the Kramers-like escape rate can be approximated as

```text
r_K(D) = r_0 exp(-DeltaU / D)
```

under the convention that `D` and `DeltaU` share the same energy-like scale. A common stochastic-resonance matching heuristic asks for roughly one transition per half forcing period,

```text
r_K(D*) ~= 2 f_signal.
```

Under the simplified rate equation above this gives

```text
D* = DeltaU / ln(r_0 / (2 f_signal))
```

only when `2 f_signal < r_0`.

This matching rule is not universal. Strong forcing, non-adiabatic response, asymmetry, colored noise and modulation-dependent prefactors can move or invalidate the simple prediction. NoiseLab therefore treats it as an oracle/control formula to test against numerical experiments, not as a law imposed on the data.

Primary references:

- L. Gammaitoni, P. Hänggi, P. Jung, F. Marchesoni, "Stochastic resonance", Reviews of Modern Physics 70, 223 (1998), DOI 10.1103/RevModPhys.70.223. https://doi.org/10.1103/RevModPhys.70.223
- L. Gammaitoni, F. Marchesoni, S. Santucci, "Stochastic Resonance as a Bona Fide Resonance", Physical Review Letters 74, 1052 (1995), DOI 10.1103/PhysRevLett.74.1052. https://doi.org/10.1103/PhysRevLett.74.1052
- D. Ryvkine, M. I. Dykman, "Noise-induced escape of periodically modulated systems: From weak to strong modulation", Physical Review E 72, 011110 (2005), DOI 10.1103/PhysRevE.72.011110. https://doi.org/10.1103/PhysRevE.72.011110

### Coherence resonance

An excitable system can exhibit maximum regularity at an intermediate noise amplitude even without a periodic input. This is a different mechanism from classical signal-plus-noise stochastic resonance and needs a different objective, for example inter-event regularity rather than input/output spectral gain.

Reference:

- A. S. Pikovsky, J. Kurths, "Coherence Resonance in a Noise-Driven Excitable System", Physical Review Letters 78, 775 (1997), DOI 10.1103/PhysRevLett.78.775. https://doi.org/10.1103/PhysRevLett.78.775

## Online optimization direction

For systems where a static grid is insufficient, NoiseLab will compare model-based prediction against model-free online extremum-seeking approaches. Extremum seeking is relevant because it estimates and tracks an input that optimizes a measured output without requiring a complete plant model, but it has its own excitation and convergence assumptions.

Recent reference used for this direction:

- "Retrospective Cost-Based Extremum Seeking Control with Vanishing Perturbation for Online Output Minimization", American Control Conference 2024, DOI 10.23919/ACC60939.2024.10644461. https://doi.org/10.23919/ACC60939.2024.10644461

## SciRust calibration targets

The pinned SciRust revision already supplies useful systems with analytic oracles:

- `scirust_sim::laser::SemiconductorLaser::relaxation_frequency()` gives a small-perturbation relaxation-oscillation frequency;
- `scirust_sim::apd::Apd::snr()` has a genuine intermediate optimal avalanche gain because signal amplification competes with excess avalanche noise;
- `scirust_sim::mechanics::SpringMassDamper` and `Pendulum` provide linear and nonlinear oscillator targets;
- `scirust_sim::electrical` supplies RLC and Van der Pol systems for electrical and nonlinear resonance studies.

The first NoiseLab operating-point tests use the SciRust APD optimum as a non-resonant control: an interior optimum must not be mislabeled as resonance merely because an optimizer finds it.

## Next experimental ladder

1. Recover the known linear oscillator resonance frequency from blind numerical sweeps.
2. Recover SciRust laser relaxation frequency from driven-response sweeps.
3. Show that the APD optimum is detected as an operating-point optimum but not automatically classified as resonance.
4. Implement a bistable Langevin system and test the Kramers half-period prediction across weak/strong forcing regimes.
5. Implement coherence-resonance controls on an excitable system.
6. Generalize the candidate space to `(time, amplitude, frequency, phase, correlation time, injection site)` and compare grid, adaptive refinement and extremum-seeking search.
7. Only after these controls pass, transfer the machinery to attention/RoPE/FLAT-ATTENTION and KV-state experiments.
