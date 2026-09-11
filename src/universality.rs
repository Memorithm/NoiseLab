//! Multiscale diagnostics for testing whether fluctuation series from different
//! mechanisms approach a common statistical representation under coarse graining.
//!
//! This module is deliberately conservative. A decreasing distance under coarse
//! graining is only a candidate universality signal. The stronger label emitted
//! here requires separation from independently shuffled surrogate controls, and
//! even that is not a proof of physical universality.

use scirust_sim::SplitMix64;
use std::error::Error;
use std::f64::consts::FRAC_PI_2;
use std::fmt::{Display, Formatter};

const MIN_COARSE_SAMPLES: usize = 8;
const MIN_SURROGATES: usize = 19;
const VARIANCE_FLOOR: f64 = 1.0e-15;

#[derive(Debug, Clone, PartialEq)]
pub enum UniversalityError {
    EmptyInput(&'static str),
    TooShort {
        name: &'static str,
        samples: usize,
        minimum: usize,
    },
    NonFinite {
        name: &'static str,
        index: usize,
        value: f64,
    },
    DegenerateSeries(&'static str),
    TooFewScales {
        provided: usize,
    },
    ZeroScale {
        index: usize,
    },
    NonIncreasingScales {
        previous: usize,
        current: usize,
    },
    ScaleTooLarge {
        scale: usize,
        samples: usize,
        minimum_coarse_samples: usize,
    },
    TooFewSurrogates {
        provided: usize,
        minimum: usize,
    },
    InvalidAlpha(f64),
}

impl Display for UniversalityError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput(name) => write!(f, "{name} must contain at least one sample"),
            Self::TooShort {
                name,
                samples,
                minimum,
            } => write!(
                f,
                "{name} has {samples} samples but at least {minimum} are required"
            ),
            Self::NonFinite {
                name,
                index,
                value,
            } => write!(f, "{name}[{index}] is not finite: {value}"),
            Self::DegenerateSeries(name) => {
                write!(f, "{name} has zero variance and cannot be standardized")
            }
            Self::TooFewScales { provided } => write!(
                f,
                "at least two coarse-graining scales are required, got {provided}"
            ),
            Self::ZeroScale { index } => write!(f, "scale at index {index} must be non-zero"),
            Self::NonIncreasingScales { previous, current } => write!(
                f,
                "coarse-graining scales must be strictly increasing: {previous} then {current}"
            ),
            Self::ScaleTooLarge {
                scale,
                samples,
                minimum_coarse_samples,
            } => write!(
                f,
                "scale {scale} leaves fewer than {minimum_coarse_samples} coarse samples from {samples} inputs"
            ),
            Self::TooFewSurrogates { provided, minimum } => write!(
                f,
                "surrogate test requires at least {minimum} repetitions, got {provided}"
            ),
            Self::InvalidAlpha(alpha) => {
                write!(f, "alpha must be finite and in (0, 1], got {alpha}")
            }
        }
    }
}

impl Error for UniversalityError {}

/// Statistical description of one coarse-grained representation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluctuationSignature {
    pub block_size: usize,
    pub coarse_samples: usize,
    pub dropped_samples: usize,
    /// Population variance after global unit-variance normalization of the raw series.
    pub variance: f64,
    /// Pearson lag-1 correlation of the coarse-grained series.
    pub lag1_autocorrelation: f64,
    /// Population excess kurtosis of the coarse-grained series.
    pub excess_kurtosis: f64,
    /// Fraction of adjacent, mean-centered samples that cross zero.
    pub sign_change_rate: f64,
}

/// Complete multiscale trace for one input series.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiscaleTrace {
    pub signatures: Vec<FluctuationSignature>,
    /// Least-squares slope of log(variance) versus log(block size).
    /// IID finite-variance noise approaches -1 for block means.
    pub variance_scaling_exponent: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScaleDistance {
    pub block_size: usize,
    pub distance: f64,
}

/// Interpretation emitted by the surrogate-gated comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UniversalityEvidence {
    /// Coarse graining did not reduce the declared descriptor distance.
    NoObservedConvergence,
    /// Distance decreased, but not more than expected from shuffled controls.
    CompatibleWithShuffledNull,
    /// Distance reduction exceeded the preregistered shuffled-surrogate gate.
    /// This is a candidate shared scaling structure, not a proof of universality.
    SurrogateSeparatedCandidate,
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniversalityTestConfig {
    pub scales: Vec<usize>,
    pub surrogate_repetitions: usize,
    pub alpha: f64,
    pub seed: u64,
}

impl Default for UniversalityTestConfig {
    fn default() -> Self {
        Self {
            scales: vec![1, 2, 4, 8, 16],
            surrogate_repetitions: 199,
            alpha: 0.05,
            seed: 0x4e4f_4953_454c_4142,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct UniversalityTestResult {
    pub trace_a: MultiscaleTrace,
    pub trace_b: MultiscaleTrace,
    pub scale_distances: Vec<ScaleDistance>,
    /// First-scale distance minus last-scale distance. Positive means convergence.
    pub convergence_score: f64,
    pub surrogate_mean_convergence: f64,
    pub surrogate_stddev_convergence: f64,
    /// One-sided randomization p-value with the standard +1 correction.
    pub empirical_p_value: f64,
    pub variance_scaling_exponent_gap: f64,
    pub evidence: UniversalityEvidence,
}

/// Compute a multiscale fluctuation trace after independently standardizing the
/// raw series to zero mean and unit population variance.
pub fn multiscale_trace(
    series: &[f64],
    scales: &[usize],
) -> Result<MultiscaleTrace, UniversalityError> {
    validate_scales(series.len(), scales)?;
    let standardized = standardize(series, "series")?;
    multiscale_trace_standardized(&standardized, scales)
}

/// Compare two fluctuation series across scales and test whether the observed
/// convergence exceeds independently shuffled surrogate controls.
///
/// Shuffling preserves each series' empirical one-point distribution while
/// destroying temporal organization. Passing this gate therefore means only
/// that the measured convergence is not explained by those marginals plus this
/// particular temporal-null model.
pub fn test_fluctuation_universality(
    series_a: &[f64],
    series_b: &[f64],
    config: &UniversalityTestConfig,
) -> Result<UniversalityTestResult, UniversalityError> {
    validate_config(series_a, series_b, config)?;

    let standardized_a = standardize(series_a, "series_a")?;
    let standardized_b = standardize(series_b, "series_b")?;
    let trace_a = multiscale_trace_standardized(&standardized_a, &config.scales)?;
    let trace_b = multiscale_trace_standardized(&standardized_b, &config.scales)?;
    let scale_distances = trace_distances(&trace_a, &trace_b);
    let convergence_score = convergence_score(&scale_distances);

    let mut rng = SplitMix64::new(config.seed);
    let mut surrogate_scores = Vec::with_capacity(config.surrogate_repetitions);

    for _ in 0..config.surrogate_repetitions {
        let mut surrogate_a = standardized_a.clone();
        let mut surrogate_b = standardized_b.clone();
        fisher_yates_shuffle(&mut surrogate_a, &mut rng);
        fisher_yates_shuffle(&mut surrogate_b, &mut rng);

        let surrogate_trace_a = multiscale_trace_standardized(&surrogate_a, &config.scales)?;
        let surrogate_trace_b = multiscale_trace_standardized(&surrogate_b, &config.scales)?;
        let distances = trace_distances(&surrogate_trace_a, &surrogate_trace_b);
        surrogate_scores.push(convergence_score(&distances));
    }

    let surrogate_mean_convergence = mean(&surrogate_scores);
    let surrogate_stddev_convergence = sample_stddev(&surrogate_scores);
    let at_least_as_extreme = surrogate_scores
        .iter()
        .filter(|&&score| score >= convergence_score)
        .count();
    let empirical_p_value = (at_least_as_extreme + 1) as f64
        / (config.surrogate_repetitions + 1) as f64;

    let evidence = if convergence_score <= 0.0 {
        UniversalityEvidence::NoObservedConvergence
    } else if empirical_p_value <= config.alpha {
        UniversalityEvidence::SurrogateSeparatedCandidate
    } else {
        UniversalityEvidence::CompatibleWithShuffledNull
    };

    Ok(UniversalityTestResult {
        variance_scaling_exponent_gap: (trace_a.variance_scaling_exponent
            - trace_b.variance_scaling_exponent)
            .abs(),
        trace_a,
        trace_b,
        scale_distances,
        convergence_score,
        surrogate_mean_convergence,
        surrogate_stddev_convergence,
        empirical_p_value,
        evidence,
    })
}

fn validate_config(
    series_a: &[f64],
    series_b: &[f64],
    config: &UniversalityTestConfig,
) -> Result<(), UniversalityError> {
    if config.surrogate_repetitions < MIN_SURROGATES {
        return Err(UniversalityError::TooFewSurrogates {
            provided: config.surrogate_repetitions,
            minimum: MIN_SURROGATES,
        });
    }
    if !config.alpha.is_finite() || config.alpha <= 0.0 || config.alpha > 1.0 {
        return Err(UniversalityError::InvalidAlpha(config.alpha));
    }
    validate_series(series_a, "series_a")?;
    validate_series(series_b, "series_b")?;
    validate_scales(series_a.len(), &config.scales)?;
    validate_scales(series_b.len(), &config.scales)?;
    Ok(())
}

fn validate_series(series: &[f64], name: &'static str) -> Result<(), UniversalityError> {
    if series.is_empty() {
        return Err(UniversalityError::EmptyInput(name));
    }
    if series.len() < MIN_COARSE_SAMPLES {
        return Err(UniversalityError::TooShort {
            name,
            samples: series.len(),
            minimum: MIN_COARSE_SAMPLES,
        });
    }
    for (index, &value) in series.iter().enumerate() {
        if !value.is_finite() {
            return Err(UniversalityError::NonFinite {
                name,
                index,
                value,
            });
        }
    }
    Ok(())
}

fn validate_scales(samples: usize, scales: &[usize]) -> Result<(), UniversalityError> {
    if scales.len() < 2 {
        return Err(UniversalityError::TooFewScales {
            provided: scales.len(),
        });
    }
    let mut previous = None;
    for (index, &scale) in scales.iter().enumerate() {
        if scale == 0 {
            return Err(UniversalityError::ZeroScale { index });
        }
        if let Some(prev) = previous {
            if scale <= prev {
                return Err(UniversalityError::NonIncreasingScales {
                    previous: prev,
                    current: scale,
                });
            }
        }
        if samples / scale < MIN_COARSE_SAMPLES {
            return Err(UniversalityError::ScaleTooLarge {
                scale,
                samples,
                minimum_coarse_samples: MIN_COARSE_SAMPLES,
            });
        }
        previous = Some(scale);
    }
    Ok(())
}

fn standardize(series: &[f64], name: &'static str) -> Result<Vec<f64>, UniversalityError> {
    validate_series(series, name)?;
    let mean_value = mean(series);
    let variance = series
        .iter()
        .map(|&value| {
            let centered = value - mean_value;
            centered * centered
        })
        .sum::<f64>()
        / series.len() as f64;
    if variance <= VARIANCE_FLOOR {
        return Err(UniversalityError::DegenerateSeries(name));
    }
    let scale = variance.sqrt();
    Ok(series
        .iter()
        .map(|&value| (value - mean_value) / scale)
        .collect())
}

fn multiscale_trace_standardized(
    standardized: &[f64],
    scales: &[usize],
) -> Result<MultiscaleTrace, UniversalityError> {
    validate_scales(standardized.len(), scales)?;
    let signatures = scales
        .iter()
        .map(|&block_size| signature_for_scale(standardized, block_size))
        .collect::<Vec<_>>();
    let variance_scaling_exponent = variance_scaling_exponent(&signatures);
    Ok(MultiscaleTrace {
        signatures,
        variance_scaling_exponent,
    })
}

fn signature_for_scale(series: &[f64], block_size: usize) -> FluctuationSignature {
    let coarse = block_mean(series, block_size);
    let dropped_samples = series.len() - coarse.len() * block_size;
    let mean_value = mean(&coarse);
    let variance = population_variance(&coarse, mean_value);
    let excess_kurtosis = if variance > VARIANCE_FLOOR {
        let fourth = coarse
            .iter()
            .map(|&value| {
                let centered = value - mean_value;
                centered.powi(4)
            })
            .sum::<f64>()
            / coarse.len() as f64;
        fourth / (variance * variance) - 3.0
    } else {
        0.0
    };

    FluctuationSignature {
        block_size,
        coarse_samples: coarse.len(),
        dropped_samples,
        variance,
        lag1_autocorrelation: lag1_correlation(&coarse),
        excess_kurtosis,
        sign_change_rate: sign_change_rate(&coarse, mean_value),
    }
}

fn block_mean(series: &[f64], block_size: usize) -> Vec<f64> {
    series
        .chunks_exact(block_size)
        .map(|chunk| chunk.iter().sum::<f64>() / block_size as f64)
        .collect()
}

fn population_variance(series: &[f64], mean_value: f64) -> f64 {
    series
        .iter()
        .map(|&value| {
            let centered = value - mean_value;
            centered * centered
        })
        .sum::<f64>()
        / series.len() as f64
}

fn lag1_correlation(series: &[f64]) -> f64 {
    if series.len() < 2 {
        return 0.0;
    }
    let left = &series[..series.len() - 1];
    let right = &series[1..];
    let left_mean = mean(left);
    let right_mean = mean(right);
    let mut numerator = 0.0;
    let mut left_energy = 0.0;
    let mut right_energy = 0.0;
    for (&x, &y) in left.iter().zip(right) {
        let dx = x - left_mean;
        let dy = y - right_mean;
        numerator += dx * dy;
        left_energy += dx * dx;
        right_energy += dy * dy;
    }
    let denominator = (left_energy * right_energy).sqrt();
    if denominator <= VARIANCE_FLOOR {
        0.0
    } else {
        (numerator / denominator).clamp(-1.0, 1.0)
    }
}

fn sign_change_rate(series: &[f64], mean_value: f64) -> f64 {
    if series.len() < 2 {
        return 0.0;
    }
    let changes = series
        .windows(2)
        .filter(|pair| {
            let left = pair[0] - mean_value;
            let right = pair[1] - mean_value;
            (left < 0.0 && right > 0.0) || (left > 0.0 && right < 0.0)
        })
        .count();
    changes as f64 / (series.len() - 1) as f64
}

fn signature_distance(a: FluctuationSignature, b: FluctuationSignature) -> f64 {
    let log_variance_gap = ((a.variance + VARIANCE_FLOOR) / (b.variance + VARIANCE_FLOOR))
        .ln();
    let autocorrelation_gap = a.lag1_autocorrelation - b.lag1_autocorrelation;
    let bounded_kurtosis_a = a.excess_kurtosis.atan() / FRAC_PI_2;
    let bounded_kurtosis_b = b.excess_kurtosis.atan() / FRAC_PI_2;
    let kurtosis_gap = bounded_kurtosis_a - bounded_kurtosis_b;
    let sign_change_gap = a.sign_change_rate - b.sign_change_rate;

    (log_variance_gap * log_variance_gap
        + autocorrelation_gap * autocorrelation_gap
        + kurtosis_gap * kurtosis_gap
        + sign_change_gap * sign_change_gap)
        .sqrt()
}

fn trace_distances(a: &MultiscaleTrace, b: &MultiscaleTrace) -> Vec<ScaleDistance> {
    a.signatures
        .iter()
        .zip(&b.signatures)
        .map(|(&left, &right)| ScaleDistance {
            block_size: left.block_size,
            distance: signature_distance(left, right),
        })
        .collect()
}

fn convergence_score(distances: &[ScaleDistance]) -> f64 {
    distances
        .first()
        .zip(distances.last())
        .map_or(0.0, |(first, last)| first.distance - last.distance)
}

fn variance_scaling_exponent(signatures: &[FluctuationSignature]) -> f64 {
    let xs = signatures
        .iter()
        .map(|signature| (signature.block_size as f64).ln())
        .collect::<Vec<_>>();
    let ys = signatures
        .iter()
        .map(|signature| (signature.variance + VARIANCE_FLOOR).ln())
        .collect::<Vec<_>>();
    let x_mean = mean(&xs);
    let y_mean = mean(&ys);
    let numerator = xs
        .iter()
        .zip(&ys)
        .map(|(&x, &y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>();
    let denominator = xs
        .iter()
        .map(|&x| {
            let dx = x - x_mean;
            dx * dx
        })
        .sum::<f64>();
    if denominator <= VARIANCE_FLOOR {
        0.0
    } else {
        numerator / denominator
    }
}

fn fisher_yates_shuffle(values: &mut [f64], rng: &mut SplitMix64) {
    for upper in (1..values.len()).rev() {
        let index = unbiased_index(rng, upper + 1);
        values.swap(upper, index);
    }
}

fn unbiased_index(rng: &mut SplitMix64, upper_exclusive: usize) -> usize {
    let bound = upper_exclusive as u64;
    let threshold = bound.wrapping_neg() % bound;
    loop {
        let value = rng.next_u64();
        if value >= threshold {
            return (value % bound) as usize;
        }
    }
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn sample_stddev(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let mean_value = mean(values);
    let variance = values
        .iter()
        .map(|&value| {
            let centered = value - mean_value;
            centered * centered
        })
        .sum::<f64>()
        / (values.len() - 1) as f64;
    variance.sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{gaussian_white_noise, ornstein_uhlenbeck_path, GaussianNoise};

    #[test]
    fn block_mean_scaling_for_white_noise_is_close_to_inverse_block_size() {
        let samples = gaussian_white_noise(GaussianNoise::new(1.0, 7).unwrap(), 65_536).unwrap();
        let trace = multiscale_trace(&samples, &[1, 2, 4, 8, 16, 32]).unwrap();
        let scale_four = trace
            .signatures
            .iter()
            .find(|signature| signature.block_size == 4)
            .unwrap();
        assert!((scale_four.variance - 0.25).abs() < 0.03);
        assert!((trace.variance_scaling_exponent + 1.0).abs() < 0.08);
    }

    #[test]
    fn identical_series_have_zero_descriptor_distance() {
        let samples = gaussian_white_noise(GaussianNoise::new(1.0, 11).unwrap(), 4096).unwrap();
        let config = UniversalityTestConfig {
            scales: vec![1, 2, 4, 8],
            surrogate_repetitions: 19,
            alpha: 0.05,
            seed: 22,
        };
        let result = test_fluctuation_universality(&samples, &samples, &config).unwrap();
        assert!(result
            .scale_distances
            .iter()
            .all(|point| point.distance.abs() < 1.0e-12));
        assert_eq!(result.evidence, UniversalityEvidence::NoObservedConvergence);
    }

    #[test]
    fn shuffled_surrogate_test_is_seed_reproducible() {
        let a = gaussian_white_noise(GaussianNoise::new(1.0, 101).unwrap(), 4096).unwrap();
        let b = ornstein_uhlenbeck_path(0.0, 0.7, 0.0, 1.0, 0.05, 4096, 202).unwrap();
        let config = UniversalityTestConfig {
            scales: vec![1, 2, 4, 8, 16],
            surrogate_repetitions: 39,
            alpha: 0.05,
            seed: 303,
        };
        let first = test_fluctuation_universality(&a, &b, &config).unwrap();
        let second = test_fluctuation_universality(&a, &b, &config).unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn invalid_scale_and_surrogate_count_fail_closed() {
        let samples = gaussian_white_noise(GaussianNoise::new(1.0, 1).unwrap(), 128).unwrap();
        let bad_scales = UniversalityTestConfig {
            scales: vec![1, 4, 4],
            surrogate_repetitions: 19,
            alpha: 0.05,
            seed: 1,
        };
        assert!(matches!(
            test_fluctuation_universality(&samples, &samples, &bad_scales),
            Err(UniversalityError::NonIncreasingScales { .. })
        ));

        let too_few_surrogates = UniversalityTestConfig {
            scales: vec![1, 2],
            surrogate_repetitions: 18,
            alpha: 0.05,
            seed: 1,
        };
        assert!(matches!(
            test_fluctuation_universality(&samples, &samples, &too_few_surrogates),
            Err(UniversalityError::TooFewSurrogates { .. })
        ));
    }
}
