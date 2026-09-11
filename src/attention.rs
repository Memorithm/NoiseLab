//! Perturbation controls for the real FLAT-ATTENTION RoPE/GQA oracle.
//!
//! This module is intentionally an integration/control layer. It does not
//! implement another attention mechanism and does not claim that perturbing an
//! attention tensor is beneficial. The clean and perturbed executions both go
//! through FLAT-ATTENTION's deterministic scalar
//! [`forward_reference_grouped_rope`] oracle.

use flat_attention::{
    forward_reference_grouped_rope, FlatAttentionConfig, FlatAttentionError, FlatAttentionOutput,
    GroupedAttentionShape, RotaryEmbeddingConfig,
};
use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Raw attention tensor receiving an additive perturbation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AttentionNoiseSite {
    /// Perturb raw query projection values before RoPE.
    Query,
    /// Perturb raw key projection values before RoPE.
    Key,
    /// Perturb value projection values. Values are never rotated by RoPE.
    Value,
}

/// Difference between one perturbed FLAT oracle execution and its clean control.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AttentionPerturbationResponse {
    /// Tensor that received the perturbation.
    pub site: AttentionNoiseSite,
    /// Standard deviation of the additive Gaussian perturbation.
    pub noise_stddev: f64,
    /// Deterministic SciRust RNG seed.
    pub seed: u64,
    /// RMS change of the context/output tensor.
    pub output_rms_delta: f64,
    /// Maximum absolute change of the context/output tensor.
    pub output_max_abs_delta: f64,
    /// RMS change of FLAT's per-query log-sum-exp statistic.
    pub lse_rms_delta: f64,
    /// Maximum absolute change of FLAT's per-query log-sum-exp statistic.
    pub lse_max_abs_delta: f64,
}

/// Failure while constructing or evaluating an attention perturbation control.
#[derive(Debug)]
pub enum AttentionExperimentError {
    /// Gaussian standard deviation was negative, NaN, or infinite.
    InvalidNoiseStddev(f64),
    /// FLAT-ATTENTION rejected the shape, inputs, configuration, or execution.
    Flat(FlatAttentionError),
}

impl Display for AttentionExperimentError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidNoiseStddev(value) => write!(
                formatter,
                "attention noise standard deviation must be finite and non-negative, got {value}"
            ),
            Self::Flat(error) => write!(formatter, "FLAT-ATTENTION oracle rejected experiment: {error}"),
        }
    }
}

impl Error for AttentionExperimentError {}

impl From<FlatAttentionError> for AttentionExperimentError {
    fn from(value: FlatAttentionError) -> Self {
        Self::Flat(value)
    }
}

/// Execute FLAT-ATTENTION's deterministic head-local RoPE/GQA oracle unchanged.
pub fn flat_rope_control(
    q: &[f32],
    k: &[f32],
    v: &[f32],
    shape: GroupedAttentionShape,
    config: FlatAttentionConfig,
    rotary: RotaryEmbeddingConfig,
) -> Result<FlatAttentionOutput, FlatAttentionError> {
    forward_reference_grouped_rope(q, k, v, shape, config, rotary)
}

/// Compare one additive Gaussian tensor perturbation with the exact clean FLAT
/// oracle execution.
///
/// The same Q/K/V inputs, shape, attention configuration and RoPE configuration
/// are used on both sides. Only the declared [`AttentionNoiseSite`] is changed.
/// `noise_stddev == 0` skips RNG draws and therefore provides an exact
/// bit-identical control path.
pub fn evaluate_flat_rope_gaussian_perturbation(
    q: &[f32],
    k: &[f32],
    v: &[f32],
    shape: GroupedAttentionShape,
    config: FlatAttentionConfig,
    rotary: RotaryEmbeddingConfig,
    site: AttentionNoiseSite,
    noise_stddev: f64,
    seed: u64,
) -> Result<AttentionPerturbationResponse, AttentionExperimentError> {
    if !noise_stddev.is_finite() || noise_stddev < 0.0 {
        return Err(AttentionExperimentError::InvalidNoiseStddev(noise_stddev));
    }

    let control = flat_rope_control(q, k, v, shape, config, rotary)?;
    let mut perturbed_q = q.to_vec();
    let mut perturbed_k = k.to_vec();
    let mut perturbed_v = v.to_vec();

    if noise_stddev > 0.0 {
        let target = match site {
            AttentionNoiseSite::Query => &mut perturbed_q,
            AttentionNoiseSite::Key => &mut perturbed_k,
            AttentionNoiseSite::Value => &mut perturbed_v,
        };
        add_gaussian_in_place(target, noise_stddev, seed);
    }

    let perturbed = flat_rope_control(
        &perturbed_q,
        &perturbed_k,
        &perturbed_v,
        shape,
        config,
        rotary,
    )?;
    let (output_rms_delta, output_max_abs_delta) = delta_metrics(&control.output, &perturbed.output);
    let (lse_rms_delta, lse_max_abs_delta) = delta_metrics(&control.lse, &perturbed.lse);

    Ok(AttentionPerturbationResponse {
        site,
        noise_stddev,
        seed,
        output_rms_delta,
        output_max_abs_delta,
        lse_rms_delta,
        lse_max_abs_delta,
    })
}

fn add_gaussian_in_place(values: &mut [f32], stddev: f64, seed: u64) {
    let mut rng = SplitMix64::new(seed);
    for value in values {
        *value += (stddev * rng.next_gaussian()) as f32;
    }
}

fn delta_metrics(control: &[f32], perturbed: &[f32]) -> (f64, f64) {
    debug_assert_eq!(control.len(), perturbed.len());
    if control.is_empty() {
        return (0.0, 0.0);
    }
    let mut squared = 0.0f64;
    let mut maximum = 0.0f64;
    for (&left, &right) in control.iter().zip(perturbed) {
        let delta = f64::from(right) - f64::from(left);
        squared += delta * delta;
        maximum = maximum.max(delta.abs());
    }
    ((squared / control.len() as f64).sqrt(), maximum)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case() -> (
        Vec<f32>,
        Vec<f32>,
        Vec<f32>,
        GroupedAttentionShape,
        FlatAttentionConfig,
        RotaryEmbeddingConfig,
    ) {
        let shape = GroupedAttentionShape {
            batch: 1,
            q_heads: 2,
            kv_heads: 1,
            seq_len: 4,
            head_dim: 4,
        };
        let q = (0..32)
            .map(|index| {
                let x = index as f32 + 1.0;
                (0.17 * x).sin() + 0.03 * x
            })
            .collect();
        let k = (0..16)
            .map(|index| {
                let x = index as f32 + 1.0;
                (0.23 * x).cos() - 0.02 * x
            })
            .collect();
        let v = (0..16)
            .map(|index| {
                let x = index as f32 + 1.0;
                (0.11 * x).sin() + (0.07 * x).cos()
            })
            .collect();
        (
            q,
            k,
            v,
            shape,
            FlatAttentionConfig {
                causal: false,
                softmax_scale: None,
            },
            RotaryEmbeddingConfig {
                theta: 10_000.0,
                position_offset: 7,
            },
        )
    }

    #[test]
    fn zero_noise_is_an_exact_control_for_every_site() {
        let (q, k, v, shape, config, rotary) = case();
        for site in [
            AttentionNoiseSite::Query,
            AttentionNoiseSite::Key,
            AttentionNoiseSite::Value,
        ] {
            let response = evaluate_flat_rope_gaussian_perturbation(
                &q, &k, &v, shape, config, rotary, site, 0.0, 42,
            )
            .unwrap();
            assert_eq!(response.output_rms_delta, 0.0);
            assert_eq!(response.output_max_abs_delta, 0.0);
            assert_eq!(response.lse_rms_delta, 0.0);
            assert_eq!(response.lse_max_abs_delta, 0.0);
        }
    }

    #[test]
    fn value_noise_changes_context_but_not_attention_normalizer() {
        let (q, k, v, shape, config, rotary) = case();
        let response = evaluate_flat_rope_gaussian_perturbation(
            &q,
            &k,
            &v,
            shape,
            config,
            rotary,
            AttentionNoiseSite::Value,
            0.15,
            11,
        )
        .unwrap();
        assert!(response.output_rms_delta > 0.0);
        assert!(response.output_max_abs_delta > 0.0);
        assert_eq!(response.lse_rms_delta, 0.0);
        assert_eq!(response.lse_max_abs_delta, 0.0);
    }

    #[test]
    fn query_and_key_noise_change_score_normalizers_on_nondegenerate_case() {
        let (q, k, v, shape, config, rotary) = case();
        for site in [AttentionNoiseSite::Query, AttentionNoiseSite::Key] {
            let response = evaluate_flat_rope_gaussian_perturbation(
                &q, &k, &v, shape, config, rotary, site, 0.15, 23,
            )
            .unwrap();
            assert!(response.output_rms_delta > 0.0);
            assert!(response.lse_rms_delta > 0.0);
        }
    }

    #[test]
    fn perturbation_response_is_exactly_reproducible_for_same_seed() {
        let (q, k, v, shape, config, rotary) = case();
        let first = evaluate_flat_rope_gaussian_perturbation(
            &q,
            &k,
            &v,
            shape,
            config,
            rotary,
            AttentionNoiseSite::Query,
            0.2,
            99,
        )
        .unwrap();
        let second = evaluate_flat_rope_gaussian_perturbation(
            &q,
            &k,
            &v,
            shape,
            config,
            rotary,
            AttentionNoiseSite::Query,
            0.2,
            99,
        )
        .unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn invalid_noise_standard_deviation_fails_closed() {
        let (q, k, v, shape, config, rotary) = case();
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                evaluate_flat_rope_gaussian_perturbation(
                    &q,
                    &k,
                    &v,
                    shape,
                    config,
                    rotary,
                    AttentionNoiseSite::Query,
                    invalid,
                    1,
                ),
                Err(AttentionExperimentError::InvalidNoiseStddev(_))
            ));
        }
    }
}
