# FLAT-ATTENTION RoPE injection control

## Status

Integration/control experiment. No benefit, resonance, robustness, or novelty claim.

## Dependency provenance

NoiseLab consumes the CPU scalar oracle from:

```text
repository = https://github.com/Memorithm/FLAT-ATTENTION.git
revision   = 4529a2079434965e13e90ddd2e98ecc88ee0cb3a
oracle     = forward_reference_grouped_rope
```

The revision is pinned in `Cargo.toml`. This experiment does not enable FLAT's optional WGPU backend.

SciRust remains pinned separately and supplies `scirust_sim::SplitMix64` for deterministic Gaussian draws.

## Purpose

Before searching for useful/noise-enhanced attention regimes, establish that NoiseLab can perturb the real FLAT RoPE/GQA reference contract without replacing or approximating its attention semantics.

The clean and perturbed paths both execute FLAT's deterministic online-softmax oracle. Q and K are raw projection values rotated by FLAT internally; V is not rotated.

## Injection sites

This first control exposes three additive Gaussian sites:

- `Query`: raw Q before FLAT's internal RoPE;
- `Key`: raw K before FLAT's internal RoPE;
- `Value`: V, which does not participate in QK score construction.

For a selected tensor element `x`, the intervention is

```text
x' = x + sigma * N(0, 1)
```

using a declared deterministic SciRust seed.

## Observables

For the perturbed pass relative to the exact clean control, NoiseLab reports:

- RMS context/output delta;
- maximum absolute context/output delta;
- RMS FLAT log-sum-exp delta;
- maximum absolute FLAT log-sum-exp delta.

These are perturbation-response metrics, not utility scores.

## Required invariants

The executable regression requires:

1. `sigma = 0` is exactly identical to the control for Q, K and V injection;
2. same site, inputs, sigma and seed produce exactly the same response record;
3. Q or K perturbation changes the score normalizer on a non-degenerate test case;
4. V perturbation changes the context output but leaves FLAT's log-sum-exp exactly unchanged, because attention scores/weights are functions of Q and K, not V;
5. negative or non-finite perturbation standard deviation fails closed.

If any invariant fails, no higher-level resonance or beneficial-noise experiment may rely on this adapter until the cause is understood.

## Next experiments after this control

Once qualified, NoiseLab can add preregistered sweeps over perturbation amplitude and injection site using task-specific utility functions. RoPE-parameter perturbations, logits, cache-state perturbations, WGPU parity, and FLAT/KVLab integration remain separate future experiments. A numerical optimum in any such sweep is not automatically a resonance mechanism.
