# BANC V888 sparse-dynamics perturbation bootstrap for NoiseLab

Status: research bootstrap. No claim of stochastic resonance, biological robustness or model improvement is made by this document.

## Scope

NoiseLab will use only FlyWire BANC v888-derived topology artifacts for this line. Raw V888 data stays external and is not committed here.

Codex currently identifies BANC v888 as Female Adult Fly Brain and Nerve Cord, snapshot 2026-05-20, with 158,262 neurons and 3,037,361 aggregated connections.

Sources:
- https://codex.flywire.ai/?dataset=banc
- https://codex.flywire.ai/faq

SciRust owns reusable graph/event primitives. NoiseLab owns controlled perturbation protocols and noise characterization.

## Scientific question

How do controlled stochastic, structural and temporal perturbations affect sparse recurrent systems, and are any observed benefits specific to V888-derived structure rather than generic sparsity or degree statistics?

## Bootstrap sequence

### NL-V888-0 — perturbation contract

Freeze perturbation representation:
- target set;
- time/window;
- seed;
- amplitude/probability;
- correlation model;
- intervention layer;
- whether perturbation is observational, structural or dynamical.

All perturbations must be replayable.

### NL-V888-1 — spike/event timing perturbations

Implement:
- event-time jitter;
- event dropout;
- burst injection;
- Poisson drive;
- renewal-process controls;
- correlated event noise.

Keep expected event counts explicit.

### NL-V888-2 — structural perturbations

Implement bounded:
- random edge dropout;
- weight perturbation;
- delay perturbation;
- node silencing;
- reciprocal-pair disruption;
- inter-module edge disruption.

Structural interventions are distinct from stochastic temporal noise.

### NL-V888-3 — matched topology panel

Every V888-derived run has:
- random sparse;
- degree-matched;
- reciprocity/module-matched control where available.

This prevents attributing generic sparse-network effects to the connectome prior.

### NL-V888-4 — information-retention diagnostics

Measure whether perturbation residuals contain information about:
- hidden system state;
- future task deficit;
- impending activity transition;
- recovery success/failure.

Use existing Noise-as-Information methodology and nulls where assumptions apply.

### NL-V888-5 — stochastic-resonance search

Search bounded noise levels only after competent no-noise baselines exist.

Required controls:
- zero noise;
- monotonic degradation control;
- random topology;
- matched topology controls;
- multiple seeds.

An interior optimum is local evidence for the declared task/system only.

### NL-V888-6 — perturbation-aware SML experiments

Inject perturbations into:
- SML graph routing;
- SBG state;
- bounded memory;
- recurrent depth controller.

Keep model semantics and perturbation machinery separate.

### NL-V888-7 — FLAT sparse-router perturbations

Test:
- candidate-mask dropout;
- mask false positives;
- mask false negatives;
- score/value perturbations after routing.

Measure both quality and admitted-key density.

### NL-V888-8 — event scheduler sensitivity

Compare whether fixed-step versus event-driven implementations respond differently to identical physical-model perturbations. Any implementation-induced difference is a correctness problem before it is a scientific observation.

### NL-V888-9 — robustness transfer

Test whether perturbation operating points transfer across:
- topology sizes;
- activity densities;
- task horizons;
- random seeds;
- V888-derived versus matched controls.

Do not infer universality from one transferable parameter.

## Outputs

Record:
- topology fingerprint;
- perturbation specification;
- event statistics before/after;
- active edges;
- task metric;
- information diagnostics;
- seed;
- exact software SHAs.

## Promotion

NoiseLab does not promote topology/model decisions directly. Stable stochastic primitives may move to SciRust. A reproducible model-level effect returns to TDI/SML for independent evaluation.
