//! Stage U1 synthetic/cross-mechanism panel for fluctuation-universality research.
//!
//! The panel is intentionally exploratory. It executes a fixed set of pairwise
//! comparisons without encoding a desired evidence outcome into the tests.
//! Results retain both fine-scale and terminal distances so an already-close
//! within-class pair is not misread as a failure merely because it has little
//! additional room to converge.

use crate::{
    gaussian_white_noise, ornstein_uhlenbeck_path, simulate_double_well,
    test_fluctuation_universality, DoubleWellLangevin, GaussianNoise, LangevinError,
    LangevinRun, NoiseInputError, UniversalityError, UniversalityEvidence,
    UniversalityTestConfig, UniversalityTestResult,
};
use std::error::Error;
use std::fmt::{Display, Formatter};

const MIN_PANEL_SAMPLES: usize = 256;
const MAX_PANEL_SAMPLES: usize = 1_000_000;

/// Fixed source identities used by the Stage U1 panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageU1Source {
    GaussianWhiteA,
    GaussianWhiteB,
    OrnsteinUhlenbeckFast,
    OrnsteinUhlenbeckSlow,
    DoubleWellLangevin,
}

/// Why a comparison is present in the preregistered panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageU1PairPurpose {
    /// Independent realizations from the same IID Gaussian family.
    WithinIidControl,
    /// Same OU mechanism with different correlation time.
    WithinOuFamily,
    /// IID versus correlated stochastic process: temporal-structure control.
    DistinctTemporalControl,
    /// Different stochastic/dynamical mechanisms: exploratory cross-mechanism pair.
    CrossMechanismExploratory,
}

/// Configuration frozen before inspecting Stage U1 outcomes.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU1PanelConfig {
    /// Samples retained per source after burn-in.
    pub samples: usize,
    /// Leading samples discarded from processes with state memory.
    pub burn_in_samples: usize,
    /// Multiscale/surrogate settings shared across pairwise tests.
    pub analysis: UniversalityTestConfig,
    /// Root seed from which deterministic source seeds are derived.
    pub data_seed: u64,
}

impl Default for StageU1PanelConfig {
    fn default() -> Self {
        Self {
            samples: 8_192,
            burn_in_samples: 1_024,
            analysis: UniversalityTestConfig::default(),
            data_seed: 0x5531_4e4f_4953_454c,
        }
    }
}

/// One retained pairwise outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU1Comparison {
    pub left: StageU1Source,
    pub right: StageU1Source,
    pub purpose: StageU1PairPurpose,
    /// Descriptor distance at the finest preregistered scale.
    pub fine_scale_distance: f64,
    /// Descriptor distance at the coarsest preregistered scale.
    pub terminal_scale_distance: f64,
    /// Fine distance minus terminal distance.
    pub convergence_score: f64,
    /// Absolute gap between log-variance scaling exponents.
    pub variance_scaling_exponent_gap: f64,
    pub empirical_p_value: f64,
    pub evidence: UniversalityEvidence,
}

/// Deterministic Stage U1 result. No aggregate "universal/not universal" bit is
/// emitted because this early panel is for calibration and hypothesis shaping.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU1PanelResult {
    pub comparisons: Vec<StageU1Comparison>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StageU1Error {
    InvalidSampleCount { samples: usize },
    BurnInTooLarge {
        burn_in_samples: usize,
        retained_samples: usize,
    },
    SampleCountOverflow,
    Noise(NoiseInputError),
    Langevin(LangevinError),
    Universality(UniversalityError),
}

impl Display for StageU1Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSampleCount { samples } => write!(
                f,
                "Stage U1 retained sample count must be in {MIN_PANEL_SAMPLES}..={MAX_PANEL_SAMPLES}, got {samples}"
            ),
            Self::BurnInTooLarge {
                burn_in_samples,
                retained_samples,
            } => write!(
                f,
                "burn-in {burn_in_samples} must not exceed retained sample count {retained_samples}"
            ),
            Self::SampleCountOverflow => f.write_str("Stage U1 total sample count overflowed usize"),
            Self::Noise(error) => write!(f, "noise source error: {error}"),
            Self::Langevin(error) => write!(f, "Langevin source error: {error}"),
            Self::Universality(error) => write!(f, "universality analysis error: {error}"),
        }
    }
}

impl Error for StageU1Error {}

impl From<NoiseInputError> for StageU1Error {
    fn from(value: NoiseInputError) -> Self {
        Self::Noise(value)
    }
}

impl From<LangevinError> for StageU1Error {
    fn from(value: LangevinError) -> Self {
        Self::Langevin(value)
    }
}

impl From<UniversalityError> for StageU1Error {
    fn from(value: UniversalityError) -> Self {
        Self::Universality(value)
    }
}

/// Execute the preregistered Stage U1 panel.
///
/// The source set contains two IID Gaussian realizations, two stationary OU
/// processes with distinct relaxation rates but matched stationary variance,
/// and one unforced bistable Langevin trajectory. Every pair is analyzed with
/// the same scale grid and surrogate count; only the surrogate RNG seed is
/// deterministically decorrelated per pair.
pub fn run_stage_u1_panel(
    config: &StageU1PanelConfig,
) -> Result<StageU1PanelResult, StageU1Error> {
    validate_panel_config(config)?;
    let total = config
        .samples
        .checked_add(config.burn_in_samples)
        .ok_or(StageU1Error::SampleCountOverflow)?;

    let gaussian_a = gaussian_white_noise(
        GaussianNoise::new(1.0, derive_seed(config.data_seed, 1))?,
        total,
    )?;
    let gaussian_b = gaussian_white_noise(
        GaussianNoise::new(1.0, derive_seed(config.data_seed, 2))?,
        total,
    )?;

    // For dX = theta(mu-X)dt + sigma dW, stationary variance is
    // sigma^2/(2*theta). Choosing sigma=sqrt(2*theta) matches it to one.
    let ou_fast_theta: f64 = 2.0;
    let ou_slow_theta: f64 = 0.25;
    let dt = 0.05;
    let ou_fast = ornstein_uhlenbeck_path(
        0.0,
        ou_fast_theta,
        0.0,
        (2.0 * ou_fast_theta).sqrt(),
        dt,
        total,
        derive_seed(config.data_seed, 3),
    )?;
    let ou_slow = ornstein_uhlenbeck_path(
        0.0,
        ou_slow_theta,
        0.0,
        (2.0 * ou_slow_theta).sqrt(),
        dt,
        total,
        derive_seed(config.data_seed, 4),
    )?;

    let double_well_model = DoubleWellLangevin::new(1.0, 1.0, 0.0, 0.1, 0.0)?;
    let run = LangevinRun::new(0.01, total, config.burn_in_samples.min(total - 1))?;
    let double_well = simulate_double_well(
        double_well_model,
        0.25,
        double_well_model.well_location(),
        run,
        derive_seed(config.data_seed, 5),
    )?
    .x;

    let gaussian_a = retain_tail(&gaussian_a, config.samples);
    let gaussian_b = retain_tail(&gaussian_b, config.samples);
    let ou_fast = retain_tail(&ou_fast, config.samples);
    let ou_slow = retain_tail(&ou_slow, config.samples);
    let double_well = retain_tail(&double_well, config.samples);

    let specifications = [
        (
            StageU1Source::GaussianWhiteA,
            StageU1Source::GaussianWhiteB,
            StageU1PairPurpose::WithinIidControl,
            gaussian_a.clone(),
            gaussian_b,
        ),
        (
            StageU1Source::OrnsteinUhlenbeckFast,
            StageU1Source::OrnsteinUhlenbeckSlow,
            StageU1PairPurpose::WithinOuFamily,
            ou_fast,
            ou_slow.clone(),
        ),
        (
            StageU1Source::GaussianWhiteA,
            StageU1Source::OrnsteinUhlenbeckSlow,
            StageU1PairPurpose::DistinctTemporalControl,
            gaussian_a,
            ou_slow.clone(),
        ),
        (
            StageU1Source::OrnsteinUhlenbeckSlow,
            StageU1Source::DoubleWellLangevin,
            StageU1PairPurpose::CrossMechanismExploratory,
            ou_slow,
            double_well,
        ),
    ];

    let mut comparisons = Vec::with_capacity(specifications.len());
    for (index, (left, right, purpose, left_series, right_series)) in
        specifications.into_iter().enumerate()
    {
        let mut analysis = config.analysis.clone();
        analysis.seed = derive_seed(config.analysis.seed, 100 + index as u64);
        let result = test_fluctuation_universality(&left_series, &right_series, &analysis)?;
        comparisons.push(compact_comparison(left, right, purpose, &result));
    }

    Ok(StageU1PanelResult { comparisons })
}

fn compact_comparison(
    left: StageU1Source,
    right: StageU1Source,
    purpose: StageU1PairPurpose,
    result: &UniversalityTestResult,
) -> StageU1Comparison {
    let fine_scale_distance = result
        .scale_distances
        .first()
        .map_or(0.0, |point| point.distance);
    let terminal_scale_distance = result
        .scale_distances
        .last()
        .map_or(0.0, |point| point.distance);
    StageU1Comparison {
        left,
        right,
        purpose,
        fine_scale_distance,
        terminal_scale_distance,
        convergence_score: result.convergence_score,
        variance_scaling_exponent_gap: result.variance_scaling_exponent_gap,
        empirical_p_value: result.empirical_p_value,
        evidence: result.evidence,
    }
}

fn validate_panel_config(config: &StageU1PanelConfig) -> Result<(), StageU1Error> {
    if config.samples < MIN_PANEL_SAMPLES || config.samples > MAX_PANEL_SAMPLES {
        return Err(StageU1Error::InvalidSampleCount {
            samples: config.samples,
        });
    }
    if config.burn_in_samples > config.samples {
        return Err(StageU1Error::BurnInTooLarge {
            burn_in_samples: config.burn_in_samples,
            retained_samples: config.samples,
        });
    }
    Ok(())
}

fn retain_tail(series: &[f64], retained: usize) -> Vec<f64> {
    series[series.len() - retained..].to_vec()
}

fn derive_seed(root: u64, stream: u64) -> u64 {
    root ^ stream.wrapping_mul(0x9e37_79b9_7f4a_7c15)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compact_test_config() -> StageU1PanelConfig {
        StageU1PanelConfig {
            samples: 2_048,
            burn_in_samples: 256,
            analysis: UniversalityTestConfig {
                scales: vec![1, 2, 4, 8, 16],
                surrogate_repetitions: 19,
                alpha: 0.05,
                seed: 123_456,
            },
            data_seed: 654_321,
        }
    }

    #[test]
    fn stage_u1_panel_is_complete_and_seed_reproducible() {
        let config = compact_test_config();
        let first = run_stage_u1_panel(&config).unwrap();
        let second = run_stage_u1_panel(&config).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.comparisons.len(), 4);
        assert_eq!(
            first.comparisons[0].purpose,
            StageU1PairPurpose::WithinIidControl
        );
        assert_eq!(
            first.comparisons[1].purpose,
            StageU1PairPurpose::WithinOuFamily
        );
        assert_eq!(
            first.comparisons[2].purpose,
            StageU1PairPurpose::DistinctTemporalControl
        );
        assert_eq!(
            first.comparisons[3].purpose,
            StageU1PairPurpose::CrossMechanismExploratory
        );
        assert!(first.comparisons.iter().all(|comparison| {
            comparison.fine_scale_distance.is_finite()
                && comparison.terminal_scale_distance.is_finite()
                && comparison.convergence_score.is_finite()
                && comparison.variance_scaling_exponent_gap.is_finite()
                && comparison.empirical_p_value.is_finite()
        }));
    }

    #[test]
    fn panel_rejects_undersampled_requests() {
        let mut config = compact_test_config();
        config.samples = 128;
        assert!(matches!(
            run_stage_u1_panel(&config),
            Err(StageU1Error::InvalidSampleCount { .. })
        ));
    }
}
