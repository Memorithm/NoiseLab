//! Stage 0.7 calibration: attention/FLAT perturbation response ↔ discrete-state MI.
//!
//! This module ties the FLAT-ATTENTION RoPE/GQA perturbation controls in
//! [`crate::attention`] to the Stage 0 equal-width histogram mutual-information
//! estimator, with an optional Stage 0.5 identity-filter audit. It is an
//! exploratory **known-answer** calibration panel: it does **not** claim that
//! attention noise is informative in production models, establish causality,
//! or discover structure in natural LLM activations.
//!
//! ## Protocol
//!
//! Under a frozen tiny FLAT config (small batch/heads/seq/dim) and frozen seeds:
//!
//! 1. For each preregistered observation family, run many deterministic FLAT
//!    clean/perturbed oracle pairs and concatenate a per-query residual feature
//!    series `R` (mean-abs output delta, or abs LSE delta).
//! 2. Declare a balanced discrete label `Z` (block-alternating bit).
//! 3. Unit-normalize `R`, then build the **positive** observation
//!    `N = R_unit + A · (2Z − 1)` that injects one labelled bit by design.
//! 4. Build a **negative** control that keeps the same `N` but Fisher–Yates
//!    shuffles `Z` under a frozen seed.
//! 5. Optionally audit Stage 0.5 identity filtering on the positive pair.
//!
//! ## KV-state handling
//!
//! FLAT's scalar `forward_reference_grouped_rope` oracle does **not** expose a
//! persistent KV-cache API. True KV-state information experiments are therefore
//! **deferred**. The panel's `KvAnalogueValueSiteOutputDelta` family uses
//! Value-site perturbation as an explicit **KV-analogue** (V participates in
//! the context write path; it does not rotate under RoPE and does not change
//! attention LSE), without inventing a fake KV cache surface.
//!
//! Classification is fail-closed against frozen MI thresholds.

use crate::attention::{
    evaluate_flat_rope_gaussian_perturbation_pair, per_query_abs_lse_delta,
    per_query_mean_abs_output_delta, AttentionExperimentError, AttentionNoiseSite,
    AttentionPerturbationSpec,
};
use crate::filtered_information::{
    filtered_histogram_mutual_information_bits, FilteredInformationError,
    FilteredMutualInformationAudit, FrequencySelectiveFilterSpec,
};
use crate::information::{histogram_mutual_information_bits, InformationError};
use crate::u2_readiness::U2_SCIRUST_REVISION;
use flat_attention::{FlatAttentionConfig, GroupedAttentionShape, RotaryEmbeddingConfig};
use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Module identity retained in every panel provenance record.
pub const ATTENTION_STATE_MODULE: &str = "attention_state_information";

/// Retained residual length for Stage 0.7.
pub const ATTENTION_STATE_SAMPLES: usize = 2_048;
/// Frozen equal-width bin count for the Stage 0 histogram estimator.
pub const ATTENTION_STATE_BINS: usize = 8;
/// Root from which per-family data and shuffle seeds are derived.
pub const ATTENTION_STATE_SEED_ROOT: u64 = 0x4130_375f_5354_4154;
/// Additive amplitude of the injected labelled bit after unit-variance
/// normalization of the base residual `R`.
pub const ATTENTION_STATE_INJECTION_AMPLITUDE: f64 = 2.0;
/// Block length for the declared balanced binary label `Z`.
pub const ATTENTION_STATE_LABEL_BLOCK: usize = 64;
/// Minimum positive-control MI (bits) required to pass known-answer.
pub const ATTENTION_STATE_POSITIVE_MI_MIN_BITS: f64 = 0.25;
/// Maximum negative-control MI (bits) allowed to pass known-answer.
pub const ATTENTION_STATE_NEGATIVE_MI_MAX_BITS: f64 = 0.08;
/// Maximum absolute identity-filter delta (bits) allowed when audited.
pub const ATTENTION_STATE_IDENTITY_ABS_DELTA_MAX_BITS: f64 = 1.0e-12;
/// Frozen Gaussian perturbation stddev applied at the declared attention site.
pub const ATTENTION_STATE_NOISE_STDDEV: f64 = 0.15;
/// FLAT revision pin retained in provenance (must match Cargo.toml).
pub const ATTENTION_STATE_FLAT_REVISION: &str = "4529a2079434965e13e90ddd2e98ecc88ee0cb3a";

const FLAT_BATCH: usize = 1;
const FLAT_Q_HEADS: usize = 2;
const FLAT_KV_HEADS: usize = 1;
const FLAT_SEQ_LEN: usize = 4;
const FLAT_HEAD_DIM: usize = 4;
const FLAT_ROTARY_THETA: f32 = 10_000.0;
const FLAT_POSITION_OFFSET: usize = 7;

/// Preregistered observation families for the Stage 0.7 known-answer panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AttentionStateFamily {
    /// Value-site mean-abs per-query output delta (explicit KV-analogue).
    KvAnalogueValueSiteOutputDelta,
    /// Query-site absolute per-query LSE delta (attention-score normalizer path).
    QuerySiteLseDelta,
    /// Key-site mean-abs per-query output delta (RoPE-rotated key path).
    KeySiteOutputDelta,
}

impl AttentionStateFamily {
    /// Stable protocol label for provenance.
    pub fn label(self) -> &'static str {
        match self {
            Self::KvAnalogueValueSiteOutputDelta => "kv_analogue_value_site_output_delta",
            Self::QuerySiteLseDelta => "query_site_lse_delta",
            Self::KeySiteOutputDelta => "key_site_output_delta",
        }
    }

    fn noise_site(self) -> AttentionNoiseSite {
        match self {
            Self::KvAnalogueValueSiteOutputDelta => AttentionNoiseSite::Value,
            Self::QuerySiteLseDelta => AttentionNoiseSite::Query,
            Self::KeySiteOutputDelta => AttentionNoiseSite::Key,
        }
    }

    fn feature_rule(self) -> &'static str {
        match self {
            Self::KvAnalogueValueSiteOutputDelta => {
                "concatenated per-query mean-abs context/output delta under Value-site Gaussian perturbation (KV-analogue; no persistent KV cache API)"
            }
            Self::QuerySiteLseDelta => {
                "concatenated per-query abs LSE delta under Query-site Gaussian perturbation"
            }
            Self::KeySiteOutputDelta => {
                "concatenated per-query mean-abs context/output delta under Key-site Gaussian perturbation"
            }
        }
    }
}

/// Declaration order of the Stage 0.7 panel (frozen).
pub const ATTENTION_STATE_PANEL_FAMILIES: [AttentionStateFamily; 3] = [
    AttentionStateFamily::KvAnalogueValueSiteOutputDelta,
    AttentionStateFamily::QuerySiteLseDelta,
    AttentionStateFamily::KeySiteOutputDelta,
];

/// Known-answer classification for one family under frozen thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttentionStateKnownAnswerClass {
    /// Positive MI above threshold and negative MI below threshold (and
    /// identity delta within bound when audited).
    KnownAnswerPassed,
    /// Estimators ran but controls missed frozen thresholds.
    KnownAnswerFailed,
}

/// Control arm identity retained in provenance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttentionStateControlKind {
    /// `N = unit(R) + A·(2Z−1)` with declared balanced `Z`.
    PositiveInjectedBit,
    /// Same `N`, Fisher–Yates-shuffled `Z` under a frozen seed.
    NegativeShuffledLabel,
}

impl AttentionStateControlKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::PositiveInjectedBit => "positive_injected_bit",
            Self::NegativeShuffledLabel => "negative_shuffled_label",
        }
    }
}

/// Fail-closed errors for the Stage 0.7 panel.
#[derive(Debug, Clone, PartialEq)]
pub enum AttentionStateInformationError {
    Information(InformationError),
    Filtered(FilteredInformationError),
    Attention(String),
    DegenerateResidual {
        family: AttentionStateFamily,
    },
    LengthMismatch {
        family: AttentionStateFamily,
        got: usize,
        expected: usize,
    },
    ProtocolBins {
        bins: usize,
    },
    NonFiniteResidual {
        family: AttentionStateFamily,
        index: usize,
    },
}

impl Display for AttentionStateInformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Information(error) => write!(f, "attention-state information error: {error}"),
            Self::Filtered(error) => write!(f, "attention-state filtered audit error: {error}"),
            Self::Attention(message) => write!(f, "attention-state FLAT error: {message}"),
            Self::DegenerateResidual { family } => {
                write!(f, "{} residual has zero empirical variance", family.label())
            }
            Self::LengthMismatch {
                family,
                got,
                expected,
            } => write!(
                f,
                "{} residual length {got} != expected {expected}",
                family.label()
            ),
            Self::ProtocolBins { bins } => write!(
                f,
                "Stage 0.7 requires at least 2 histogram bins, got {bins}"
            ),
            Self::NonFiniteResidual { family, index } => write!(
                f,
                "{} residual became non-finite at index {index}",
                family.label()
            ),
        }
    }
}

impl Error for AttentionStateInformationError {}

impl From<InformationError> for AttentionStateInformationError {
    fn from(value: InformationError) -> Self {
        Self::Information(value)
    }
}

impl From<FilteredInformationError> for AttentionStateInformationError {
    fn from(value: FilteredInformationError) -> Self {
        Self::Filtered(value)
    }
}

impl From<AttentionExperimentError> for AttentionStateInformationError {
    fn from(value: AttentionExperimentError) -> Self {
        Self::Attention(value.to_string())
    }
}

/// Provenance attached to one family calibration result.
#[derive(Debug, Clone, PartialEq)]
pub struct AttentionStateProvenance {
    pub module: &'static str,
    pub family: AttentionStateFamily,
    pub family_label: &'static str,
    pub scirust_revision: &'static str,
    pub flat_revision: &'static str,
    pub data_seed: u64,
    pub shuffle_seed: u64,
    pub samples: usize,
    pub bins: usize,
    pub injection_amplitude: f64,
    pub label_block: usize,
    pub noise_stddev: f64,
    pub noise_site: &'static str,
    pub feature_rule: &'static str,
    pub flat_shape_summary: &'static str,
    pub kv_handling: &'static str,
    pub positive_mi_min_bits: f64,
    pub negative_mi_max_bits: f64,
    pub identity_abs_delta_max_bits: f64,
}

/// One family's positive/negative known-answer MI calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct AttentionStateFamilyResult {
    pub family: AttentionStateFamily,
    pub positive_mutual_information_bits: f64,
    pub negative_mutual_information_bits: f64,
    pub hidden_state_entropy_bits: f64,
    pub identity_audit: Option<FilteredMutualInformationAudit>,
    pub classification: AttentionStateKnownAnswerClass,
    pub provenance: AttentionStateProvenance,
}

/// Full Stage 0.7 panel over the preregistered attention observation families.
#[derive(Debug, Clone, PartialEq)]
pub struct AttentionStatePanelResult {
    pub families: Vec<AttentionStateFamilyResult>,
    pub all_passed: bool,
    pub identity_audit_enabled: bool,
}

/// Run the Stage 0.7 attention-state known-answer panel.
pub fn run_attention_state_information_panel(
    identity_audit: bool,
) -> Result<AttentionStatePanelResult, AttentionStateInformationError> {
    run_attention_state_information_panel_with_root(ATTENTION_STATE_SEED_ROOT, identity_audit)
}

/// Same panel with an explicit seed root (tests / deterministic overrides).
pub fn run_attention_state_information_panel_with_root(
    seed_root: u64,
    identity_audit: bool,
) -> Result<AttentionStatePanelResult, AttentionStateInformationError> {
    if ATTENTION_STATE_BINS < 2 {
        return Err(AttentionStateInformationError::ProtocolBins {
            bins: ATTENTION_STATE_BINS,
        });
    }

    let mut families = Vec::with_capacity(ATTENTION_STATE_PANEL_FAMILIES.len());
    for (index, family) in ATTENTION_STATE_PANEL_FAMILIES.iter().copied().enumerate() {
        let data_seed = derive_seed(seed_root, (index as u64).wrapping_add(1));
        let shuffle_seed = derive_seed(seed_root, (index as u64).wrapping_add(101));
        families.push(calibrate_family(
            family,
            data_seed,
            shuffle_seed,
            identity_audit,
        )?);
    }
    let all_passed = families
        .iter()
        .all(|result| result.classification == AttentionStateKnownAnswerClass::KnownAnswerPassed);
    Ok(AttentionStatePanelResult {
        families,
        all_passed,
        identity_audit_enabled: identity_audit,
    })
}

/// Calibrate a single attention observation family under the frozen protocol.
pub fn calibrate_attention_state_family(
    family: AttentionStateFamily,
    identity_audit: bool,
) -> Result<AttentionStateFamilyResult, AttentionStateInformationError> {
    let index = ATTENTION_STATE_PANEL_FAMILIES
        .iter()
        .position(|&candidate| candidate == family)
        .expect("panel families cover every AttentionStateFamily variant");
    let data_seed = derive_seed(ATTENTION_STATE_SEED_ROOT, (index as u64).wrapping_add(1));
    let shuffle_seed = derive_seed(ATTENTION_STATE_SEED_ROOT, (index as u64).wrapping_add(101));
    calibrate_family(family, data_seed, shuffle_seed, identity_audit)
}

fn calibrate_family(
    family: AttentionStateFamily,
    data_seed: u64,
    shuffle_seed: u64,
    identity_audit: bool,
) -> Result<AttentionStateFamilyResult, AttentionStateInformationError> {
    let base_residual = generate_attention_residual(family, data_seed)?;
    let unit_residual = unit_variance_residual(family, &base_residual)?;
    let labels = declared_balanced_labels(unit_residual.len());
    let positive_observation =
        inject_labelled_bit(&unit_residual, &labels, ATTENTION_STATE_INJECTION_AMPLITUDE);
    let mut negative_labels = labels.clone();
    fisher_yates_shuffle(&mut negative_labels, &mut SplitMix64::new(shuffle_seed));

    let positive_mi =
        histogram_mutual_information_bits(&positive_observation, &labels, ATTENTION_STATE_BINS)?;
    let negative_mi = histogram_mutual_information_bits(
        &positive_observation,
        &negative_labels,
        ATTENTION_STATE_BINS,
    )?;
    let hidden_state_entropy_bits = crate::discrete_entropy_bits(&labels)?;

    let identity = if identity_audit {
        Some(filtered_histogram_mutual_information_bits(
            &positive_observation,
            &labels,
            ATTENTION_STATE_BINS,
            &FrequencySelectiveFilterSpec::Identity,
        )?)
    } else {
        None
    };

    let identity_ok = match &identity {
        None => true,
        Some(audit) => {
            audit.information_delta_bits.abs() <= ATTENTION_STATE_IDENTITY_ABS_DELTA_MAX_BITS
        }
    };
    let classification = if positive_mi >= ATTENTION_STATE_POSITIVE_MI_MIN_BITS
        && negative_mi <= ATTENTION_STATE_NEGATIVE_MI_MAX_BITS
        && identity_ok
    {
        AttentionStateKnownAnswerClass::KnownAnswerPassed
    } else {
        AttentionStateKnownAnswerClass::KnownAnswerFailed
    };

    Ok(AttentionStateFamilyResult {
        family,
        positive_mutual_information_bits: positive_mi,
        negative_mutual_information_bits: negative_mi,
        hidden_state_entropy_bits,
        identity_audit: identity,
        classification,
        provenance: AttentionStateProvenance {
            module: ATTENTION_STATE_MODULE,
            family,
            family_label: family.label(),
            scirust_revision: U2_SCIRUST_REVISION,
            flat_revision: ATTENTION_STATE_FLAT_REVISION,
            data_seed,
            shuffle_seed,
            samples: ATTENTION_STATE_SAMPLES,
            bins: ATTENTION_STATE_BINS,
            injection_amplitude: ATTENTION_STATE_INJECTION_AMPLITUDE,
            label_block: ATTENTION_STATE_LABEL_BLOCK,
            noise_stddev: ATTENTION_STATE_NOISE_STDDEV,
            noise_site: match family.noise_site() {
                AttentionNoiseSite::Query => "query",
                AttentionNoiseSite::Key => "key",
                AttentionNoiseSite::Value => "value",
            },
            feature_rule: family.feature_rule(),
            flat_shape_summary: "batch=1,q_heads=2,kv_heads=1,seq_len=4,head_dim=4,causal=false,theta=10000,position_offset=7",
            kv_handling: "true_kv_cache_deferred; value_site_is_explicit_kv_analogue",
            positive_mi_min_bits: ATTENTION_STATE_POSITIVE_MI_MIN_BITS,
            negative_mi_max_bits: ATTENTION_STATE_NEGATIVE_MI_MAX_BITS,
            identity_abs_delta_max_bits: ATTENTION_STATE_IDENTITY_ABS_DELTA_MAX_BITS,
        },
    })
}

fn frozen_flat_shape() -> GroupedAttentionShape {
    GroupedAttentionShape {
        batch: FLAT_BATCH,
        q_heads: FLAT_Q_HEADS,
        kv_heads: FLAT_KV_HEADS,
        seq_len: FLAT_SEQ_LEN,
        head_dim: FLAT_HEAD_DIM,
    }
}

fn frozen_flat_config() -> FlatAttentionConfig {
    FlatAttentionConfig {
        causal: false,
        softmax_scale: None,
    }
}

fn frozen_rotary() -> RotaryEmbeddingConfig {
    RotaryEmbeddingConfig {
        theta: FLAT_ROTARY_THETA,
        position_offset: FLAT_POSITION_OFFSET,
    }
}

fn features_per_trial() -> usize {
    FLAT_BATCH * FLAT_Q_HEADS * FLAT_SEQ_LEN
}

fn trials_needed() -> usize {
    let per = features_per_trial();
    debug_assert!(ATTENTION_STATE_SAMPLES.is_multiple_of(per));
    ATTENTION_STATE_SAMPLES / per
}

fn generate_attention_residual(
    family: AttentionStateFamily,
    data_seed: u64,
) -> Result<Vec<f64>, AttentionStateInformationError> {
    let shape = frozen_flat_shape();
    let config = frozen_flat_config();
    let rotary = frozen_rotary();
    let mut residual = Vec::with_capacity(ATTENTION_STATE_SAMPLES);
    let trials = trials_needed();

    for trial in 0..trials {
        let trial_seed = derive_seed(data_seed, (trial as u64).wrapping_add(1));
        let (q, k, v) = synthesize_qkv(shape, trial_seed);
        let perturbation_seed = derive_seed(data_seed, (trial as u64).wrapping_add(10_001));
        let pair = evaluate_flat_rope_gaussian_perturbation_pair(
            &q,
            &k,
            &v,
            shape,
            config,
            rotary,
            AttentionPerturbationSpec {
                site: family.noise_site(),
                noise_stddev: ATTENTION_STATE_NOISE_STDDEV,
                seed: perturbation_seed,
            },
        )?;
        let chunk = match family {
            AttentionStateFamily::QuerySiteLseDelta => {
                per_query_abs_lse_delta(&pair.control.lse, &pair.perturbed.lse, shape)?
            }
            AttentionStateFamily::KvAnalogueValueSiteOutputDelta
            | AttentionStateFamily::KeySiteOutputDelta => per_query_mean_abs_output_delta(
                &pair.control.output,
                &pair.perturbed.output,
                shape,
            )?,
        };
        residual.extend(chunk);
    }

    finish_residual(family, residual)
}

fn synthesize_qkv(shape: GroupedAttentionShape, seed: u64) -> (Vec<f32>, Vec<f32>, Vec<f32>) {
    let q_len = shape.batch * shape.q_heads * shape.seq_len * shape.head_dim;
    let kv_len = shape.batch * shape.kv_heads * shape.seq_len * shape.head_dim;
    let mut rng = SplitMix64::new(seed);
    let mut q = Vec::with_capacity(q_len);
    let mut k = Vec::with_capacity(kv_len);
    let mut v = Vec::with_capacity(kv_len);
    for index in 0..q_len {
        let x = index as f32 + 1.0;
        let gaussian = rng.next_gaussian() as f32;
        q.push((0.17 * x).sin() + 0.03 * x + 0.05 * gaussian);
    }
    for index in 0..kv_len {
        let x = index as f32 + 1.0;
        let gaussian = rng.next_gaussian() as f32;
        k.push((0.23 * x).cos() - 0.02 * x + 0.05 * gaussian);
        let gaussian_v = rng.next_gaussian() as f32;
        v.push((0.11 * x).sin() + (0.07 * x).cos() + 0.05 * gaussian_v);
    }
    (q, k, v)
}

fn finish_residual(
    family: AttentionStateFamily,
    residual: Vec<f64>,
) -> Result<Vec<f64>, AttentionStateInformationError> {
    if residual.len() != ATTENTION_STATE_SAMPLES {
        return Err(AttentionStateInformationError::LengthMismatch {
            family,
            got: residual.len(),
            expected: ATTENTION_STATE_SAMPLES,
        });
    }
    for (index, value) in residual.iter().enumerate() {
        if !value.is_finite() {
            return Err(AttentionStateInformationError::NonFiniteResidual { family, index });
        }
    }
    let mean = residual.iter().sum::<f64>() / residual.len() as f64;
    let variance = residual
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / residual.len() as f64;
    if variance <= 1.0e-15 {
        return Err(AttentionStateInformationError::DegenerateResidual { family });
    }
    Ok(residual)
}

fn declared_balanced_labels(len: usize) -> Vec<usize> {
    (0..len)
        .map(|index| (index / ATTENTION_STATE_LABEL_BLOCK) % 2)
        .collect()
}

fn unit_variance_residual(
    family: AttentionStateFamily,
    residual: &[f64],
) -> Result<Vec<f64>, AttentionStateInformationError> {
    let n = residual.len() as f64;
    let mean = residual.iter().sum::<f64>() / n;
    let variance = residual
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / n;
    if !variance.is_finite() || variance <= 1.0e-15 {
        return Err(AttentionStateInformationError::DegenerateResidual { family });
    }
    let std = variance.sqrt();
    Ok(residual.iter().map(|value| (value - mean) / std).collect())
}

fn inject_labelled_bit(residual: &[f64], labels: &[usize], amplitude: f64) -> Vec<f64> {
    residual
        .iter()
        .zip(labels)
        .map(|(&value, &label)| {
            let signed = if label == 0 { -1.0 } else { 1.0 };
            value + amplitude * signed
        })
        .collect()
}

fn derive_seed(root: u64, stream: u64) -> u64 {
    root ^ stream.wrapping_mul(0x9e37_79b9_7f4a_7c15)
}

fn fisher_yates_shuffle(values: &mut [usize], rng: &mut SplitMix64) {
    for index in (1..values.len()).rev() {
        let swap_with = unbiased_index(rng, index + 1);
        values.swap(index, swap_with);
    }
}

fn unbiased_index(rng: &mut SplitMix64, upper_exclusive: usize) -> usize {
    debug_assert!(upper_exclusive > 0);
    let bound = upper_exclusive as u64;
    let threshold = bound.wrapping_neg() % bound;
    loop {
        let value = rng.next_u64();
        if value >= threshold {
            return (value % bound) as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::information::InformationError;

    #[test]
    fn panel_with_identity_audit_passes_all_families() {
        let panel = run_attention_state_information_panel(true).unwrap();
        assert!(panel.identity_audit_enabled);
        assert_eq!(panel.families.len(), 3);
        assert!(panel.all_passed);
        for result in &panel.families {
            assert_eq!(
                result.classification,
                AttentionStateKnownAnswerClass::KnownAnswerPassed
            );
            assert!(
                result.positive_mutual_information_bits >= ATTENTION_STATE_POSITIVE_MI_MIN_BITS,
                "{} positive MI={}",
                result.family.label(),
                result.positive_mutual_information_bits
            );
            assert!(
                result.negative_mutual_information_bits <= ATTENTION_STATE_NEGATIVE_MI_MAX_BITS,
                "{} negative MI={}",
                result.family.label(),
                result.negative_mutual_information_bits
            );
            let audit = result.identity_audit.as_ref().unwrap();
            assert!(
                audit.information_delta_bits.abs() <= ATTENTION_STATE_IDENTITY_ABS_DELTA_MAX_BITS
            );
            assert_eq!(result.provenance.module, ATTENTION_STATE_MODULE);
            assert_eq!(result.provenance.bins, ATTENTION_STATE_BINS);
            assert_eq!(result.provenance.scirust_revision, U2_SCIRUST_REVISION);
            assert_eq!(
                result.provenance.flat_revision,
                ATTENTION_STATE_FLAT_REVISION
            );
            assert!(result
                .provenance
                .kv_handling
                .contains("true_kv_cache_deferred"));
            assert!((result.hidden_state_entropy_bits - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn panel_is_seed_reproducible() {
        let a = run_attention_state_information_panel(true).unwrap();
        let b = run_attention_state_information_panel(true).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn distinct_roots_change_scores() {
        let a = run_attention_state_information_panel_with_root(1, false).unwrap();
        let b = run_attention_state_information_panel_with_root(2, false).unwrap();
        assert_ne!(
            a.families[0].provenance.data_seed,
            b.families[0].provenance.data_seed
        );
        assert_ne!(
            (
                a.families[0].positive_mutual_information_bits.to_bits(),
                a.families[0].negative_mutual_information_bits.to_bits()
            ),
            (
                b.families[0].positive_mutual_information_bits.to_bits(),
                b.families[0].negative_mutual_information_bits.to_bits()
            )
        );
    }

    #[test]
    fn positive_control_exceeds_negative_for_each_family() {
        for family in ATTENTION_STATE_PANEL_FAMILIES {
            let result = calibrate_attention_state_family(family, false).unwrap();
            assert!(
                result.positive_mutual_information_bits
                    > result.negative_mutual_information_bits + 0.1,
                "{}: pos={}, neg={}",
                family.label(),
                result.positive_mutual_information_bits,
                result.negative_mutual_information_bits
            );
        }
    }

    #[test]
    fn real_oracle_is_invoked_for_value_site_kv_analogue() {
        // Smoke path still exercises the real FLAT pair API for one trial.
        let shape = frozen_flat_shape();
        let (q, k, v) = synthesize_qkv(shape, 42);
        let pair = evaluate_flat_rope_gaussian_perturbation_pair(
            &q,
            &k,
            &v,
            shape,
            frozen_flat_config(),
            frozen_rotary(),
            AttentionPerturbationSpec {
                site: AttentionNoiseSite::Value,
                noise_stddev: ATTENTION_STATE_NOISE_STDDEV,
                seed: 7,
            },
        )
        .unwrap();
        assert!(pair.response.output_rms_delta > 0.0);
        assert_eq!(pair.response.lse_rms_delta, 0.0);
        let series =
            per_query_mean_abs_output_delta(&pair.control.output, &pair.perturbed.output, shape)
                .unwrap();
        assert_eq!(series.len(), features_per_trial());
        assert!(series.iter().any(|&value| value > 0.0));
    }

    #[test]
    fn injection_and_shuffle_helpers_are_length_preserving() {
        let residual = vec![0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8];
        let labels = declared_balanced_labels(residual.len());
        assert_eq!(labels.len(), residual.len());
        let injected = inject_labelled_bit(&residual, &labels, 2.0);
        assert_eq!(injected.len(), residual.len());
        let mut shuffled = labels.clone();
        fisher_yates_shuffle(&mut shuffled, &mut SplitMix64::new(99));
        assert_eq!(shuffled.len(), labels.len());
        let mut sorted_a = labels.clone();
        let mut sorted_b = shuffled.clone();
        sorted_a.sort_unstable();
        sorted_b.sort_unstable();
        assert_eq!(sorted_a, sorted_b);
    }

    #[test]
    fn malformed_estimator_inputs_fail_closed() {
        assert_eq!(
            histogram_mutual_information_bits(&[], &[], ATTENTION_STATE_BINS),
            Err(InformationError::EmptyInput)
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0, 1.0], &[0], ATTENTION_STATE_BINS),
            Err(InformationError::LengthMismatch {
                observed: 2,
                hidden_state: 1,
            })
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0, 1.0], &[0, 1], 1),
            Err(InformationError::TooFewBins { bins: 1 })
        );
    }

    #[test]
    fn length_mismatch_on_finish_is_protocol_error() {
        let err = finish_residual(
            AttentionStateFamily::KvAnalogueValueSiteOutputDelta,
            vec![1.0, 2.0],
        )
        .unwrap_err();
        assert!(matches!(
            err,
            AttentionStateInformationError::LengthMismatch {
                family: AttentionStateFamily::KvAnalogueValueSiteOutputDelta,
                got: 2,
                expected: ATTENTION_STATE_SAMPLES,
            }
        ));
    }

    #[test]
    fn degenerate_residual_fails_closed() {
        let err = finish_residual(
            AttentionStateFamily::QuerySiteLseDelta,
            vec![0.0; ATTENTION_STATE_SAMPLES],
        )
        .unwrap_err();
        assert!(matches!(
            err,
            AttentionStateInformationError::DegenerateResidual {
                family: AttentionStateFamily::QuerySiteLseDelta,
            }
        ));
    }

    #[test]
    fn control_kind_labels_are_stable() {
        assert_eq!(
            AttentionStateControlKind::PositiveInjectedBit.label(),
            "positive_injected_bit"
        );
        assert_eq!(
            AttentionStateControlKind::NegativeShuffledLabel.label(),
            "negative_shuffled_label"
        );
    }

    #[test]
    fn samples_divide_features_per_trial() {
        assert!(ATTENTION_STATE_SAMPLES.is_multiple_of(features_per_trial()));
        assert_eq!(
            trials_needed() * features_per_trial(),
            ATTENTION_STATE_SAMPLES
        );
    }
}
