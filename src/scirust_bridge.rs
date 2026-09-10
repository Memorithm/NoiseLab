//! Thin, provenance-preserving adapters over SciRust foundations.
//!
//! NoiseLab deliberately depends on pinned SciRust revisions instead of
//! copying numerical infrastructure into this repository.  This module is the
//! narrow boundary through which the first NoiseLab experiments consume those
//! upstream primitives.

use scirust_signal::{hanning, psd_centroid, psd_flatness, psd_spread, welch_psd};
use scirust_sim::{SplitMix64, stochastic::ou_path};
use std::error::Error;
use std::fmt::{Display, Formatter};

const MAX_SAMPLES: usize = 10_000_000;

/// Invalid request supplied to a NoiseLab/SciRust adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum NoiseInputError {
    /// A floating-point input was NaN or infinite.
    NonFinite(&'static str),
    /// A parameter that must be strictly positive was zero or negative.
    NonPositive(&'static str),
    /// The requested allocation exceeds the research-bench safety budget.
    TooManySamples { requested: usize, maximum: usize },
    /// A spectral request is structurally invalid.
    InvalidSpectrumRequest(&'static str),
    /// SciRust rejected the request.
    Upstream(String),
}

impl Display for NoiseInputError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(f, "{name} must be finite"),
            Self::NonPositive(name) => write!(f, "{name} must be strictly positive"),
            Self::TooManySamples { requested, maximum } => {
                write!(f, "requested {requested} samples exceeds safety maximum {maximum}")
            }
            Self::InvalidSpectrumRequest(message) => f.write_str(message),
            Self::Upstream(message) => write!(f, "SciRust rejected request: {message}"),
        }
    }
}

impl Error for NoiseInputError {}

/// Parameters for a deterministic zero-mean Gaussian white-noise realization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianNoise {
    /// Standard deviation of the generated process.
    pub sigma: f64,
    /// Explicit reproducibility seed.
    pub seed: u64,
}

impl GaussianNoise {
    /// Validate and construct Gaussian white-noise parameters.
    pub fn new(sigma: f64, seed: u64) -> Result<Self, NoiseInputError> {
        if !sigma.is_finite() {
            return Err(NoiseInputError::NonFinite("sigma"));
        }
        if sigma < 0.0 {
            return Err(NoiseInputError::NonPositive("sigma"));
        }
        Ok(Self { sigma, seed })
    }
}

/// Generate a deterministic Gaussian white-noise realization using SciRust's
/// validated `SplitMix64` + Box-Muller implementation.
pub fn gaussian_white_noise(
    params: GaussianNoise,
    samples: usize,
) -> Result<Vec<f64>, NoiseInputError> {
    validate_sample_count(samples)?;
    let mut rng = SplitMix64::new(params.seed);
    Ok((0..samples)
        .map(|_| params.sigma * rng.next_gaussian())
        .collect())
}

/// Generate an Ornstein-Uhlenbeck path through SciRust's exact transition law.
///
/// This is useful as the first correlated-noise/control process because the
/// upstream implementation has an analytic stationary variance oracle.
pub fn ornstein_uhlenbeck_path(
    x0: f64,
    theta: f64,
    mean: f64,
    sigma: f64,
    dt: f64,
    steps: usize,
    seed: u64,
) -> Result<Vec<f64>, NoiseInputError> {
    if steps > MAX_SAMPLES {
        return Err(NoiseInputError::TooManySamples {
            requested: steps,
            maximum: MAX_SAMPLES,
        });
    }
    ou_path(x0, theta, mean, sigma, dt, steps, seed)
        .map_err(|error| NoiseInputError::Upstream(error.to_string()))
}

/// Compact spectral characterization produced entirely by `scirust-signal`.
#[derive(Debug, Clone, PartialEq)]
pub struct ScirustSpectralSignature {
    /// Welch-averaged one-sided power spectral density.
    pub psd: Vec<f64>,
    /// Number of averaged Welch segments.
    pub segments: usize,
    /// Number of trailing samples not used by Welch segmentation.
    pub dropped_samples: usize,
    /// PSD-weighted center frequency in hertz.
    pub centroid_hz: f64,
    /// PSD-weighted spread around the center frequency in hertz.
    pub spread_hz: f64,
    /// Wiener spectral flatness over power, in `[0, 1]` for non-negative PSDs.
    pub flatness: f64,
}

/// Characterize a real-valued series with SciRust's Hann-windowed Welch PSD.
pub fn spectral_signature(
    signal: &[f64],
    sample_rate_hz: f64,
    segment_len: usize,
    overlap: usize,
) -> Result<ScirustSpectralSignature, NoiseInputError> {
    if !sample_rate_hz.is_finite() {
        return Err(NoiseInputError::NonFinite("sample_rate_hz"));
    }
    if sample_rate_hz <= 0.0 {
        return Err(NoiseInputError::NonPositive("sample_rate_hz"));
    }
    if segment_len < 2 || !segment_len.is_power_of_two() {
        return Err(NoiseInputError::InvalidSpectrumRequest(
            "segment_len must be a power of two greater than or equal to 2",
        ));
    }
    if overlap >= segment_len {
        return Err(NoiseInputError::InvalidSpectrumRequest(
            "overlap must be strictly less than segment_len",
        ));
    }

    let window = hanning(segment_len);
    let estimate = welch_psd(signal, segment_len, overlap, &window)
        .map_err(|error| NoiseInputError::Upstream(error.to_string()))?;

    let centroid_hz = psd_centroid(&estimate.psd, sample_rate_hz);
    let spread_hz = psd_spread(&estimate.psd, sample_rate_hz);
    let flatness = psd_flatness(&estimate.psd);

    Ok(ScirustSpectralSignature {
        psd: estimate.psd,
        segments: estimate.segments,
        dropped_samples: estimate.dropped_samples,
        centroid_hz,
        spread_hz,
        flatness,
    })
}

fn validate_sample_count(samples: usize) -> Result<(), NoiseInputError> {
    if samples > MAX_SAMPLES {
        return Err(NoiseInputError::TooManySamples {
            requested: samples,
            maximum: MAX_SAMPLES,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_generation_is_seed_reproducible() {
        let params = GaussianNoise::new(0.75, 42).unwrap();
        let a = gaussian_white_noise(params, 256).unwrap();
        let b = gaussian_white_noise(params, 256).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn gaussian_zero_sigma_is_exact_zero() {
        let params = GaussianNoise::new(0.0, 7).unwrap();
        let samples = gaussian_white_noise(params, 64).unwrap();
        assert!(samples.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn different_seeds_change_gaussian_realization() {
        let a = gaussian_white_noise(GaussianNoise::new(1.0, 1).unwrap(), 64).unwrap();
        let b = gaussian_white_noise(GaussianNoise::new(1.0, 2).unwrap(), 64).unwrap();
        assert_ne!(a, b);
    }

    #[test]
    fn ou_wrapper_keeps_upstream_reproducibility() {
        let a = ornstein_uhlenbeck_path(0.0, 0.8, 1.0, 0.3, 0.05, 128, 99).unwrap();
        let b = ornstein_uhlenbeck_path(0.0, 0.8, 1.0, 0.3, 0.05, 128, 99).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn spectral_signature_reports_welch_provenance() {
        let samples = gaussian_white_noise(GaussianNoise::new(1.0, 1234).unwrap(), 1024).unwrap();
        let signature = spectral_signature(&samples, 1024.0, 256, 128).unwrap();
        assert_eq!(signature.segments, 7);
        assert_eq!(signature.dropped_samples, 0);
        assert_eq!(signature.psd.len(), 129);
        assert!(signature.centroid_hz.is_finite());
        assert!(signature.spread_hz.is_finite());
        assert!(signature.flatness.is_finite());
    }
}
