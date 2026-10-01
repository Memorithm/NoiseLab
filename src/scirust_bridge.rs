//! Thin, provenance-preserving adapters over SciRust foundations.
//!
//! NoiseLab deliberately depends on pinned SciRust revisions instead of
//! copying numerical infrastructure into this repository.  This module is the
//! narrow boundary through which the first NoiseLab experiments consume those
//! upstream primitives.

use scirust_signal::{psd_centroid, psd_flatness, psd_spread, welch_psd};
use scirust_sim::{stochastic::ou_path, SplitMix64};
use std::error::Error;
use std::f64::consts::PI;
use std::fmt::{Display, Formatter};
use std::mem::size_of;

const MAX_SAMPLES: usize = 10_000_000;
const DEFAULT_MAX_SPECTRAL_SEGMENT_SAMPLES: usize = 1 << 20;
const DEFAULT_MAX_SPECTRAL_WORKING_BYTES: usize = 64 * 1024 * 1024;
const DEFAULT_MAX_SPECTRAL_WORK_UNITS: usize = 1_000_000_000;

/// Admission limits for a Welch spectral characterization.
///
/// `max_working_bytes` bounds the peak allocations implied by the pinned
/// SciRust implementation (window, accumulator, shaped segment, complex FFT
/// buffer and per-segment PSD). `max_work_units` bounds the number of samples
/// scanned plus the conservative radix-2 FFT work estimate over all segments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpectralLimits {
    /// Maximum number of input samples read by one request.
    pub max_signal_samples: usize,
    /// Maximum number of samples in one FFT segment.
    pub max_segment_samples: usize,
    /// Maximum admitted peak working-set bytes.
    pub max_working_bytes: usize,
    /// Maximum admitted conservative work units.
    pub max_work_units: usize,
}

impl Default for SpectralLimits {
    fn default() -> Self {
        Self {
            max_signal_samples: MAX_SAMPLES,
            max_segment_samples: DEFAULT_MAX_SPECTRAL_SEGMENT_SAMPLES,
            max_working_bytes: DEFAULT_MAX_SPECTRAL_WORKING_BYTES,
            max_work_units: DEFAULT_MAX_SPECTRAL_WORK_UNITS,
        }
    }
}

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
    /// A spectral request exceeds one of the configured admission limits.
    SpectralLimitExceeded {
        resource: &'static str,
        requested: usize,
        maximum: usize,
    },
    /// A spectral resource estimate overflowed and was rejected fail-closed.
    SpectralBudgetOverflow(&'static str),
    /// A bounded adapter-owned allocation failed.
    AllocationFailed {
        resource: &'static str,
        requested_bytes: usize,
    },
    /// A generated numeric output overflowed to NaN or infinity.
    NonFiniteOutput {
        operation: &'static str,
        sample_index: usize,
    },
    /// SciRust rejected the request.
    Upstream(String),
}

impl Display for NoiseInputError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(f, "{name} must be finite"),
            Self::NonPositive(name) => write!(f, "{name} must be strictly positive"),
            Self::TooManySamples { requested, maximum } => {
                write!(
                    f,
                    "requested {requested} samples exceeds safety maximum {maximum}"
                )
            }
            Self::InvalidSpectrumRequest(message) => f.write_str(message),
            Self::SpectralLimitExceeded {
                resource,
                requested,
                maximum,
            } => write!(
                f,
                "spectral {resource} request {requested} exceeds configured maximum {maximum}"
            ),
            Self::SpectralBudgetOverflow(resource) => {
                write!(f, "spectral {resource} estimate overflowed")
            }
            Self::AllocationFailed {
                resource,
                requested_bytes,
            } => write!(
                f,
                "failed to allocate {requested_bytes} bytes for {resource}"
            ),
            Self::NonFiniteOutput {
                operation,
                sample_index,
            } => write!(
                f,
                "{operation} produced a non-finite value at sample {sample_index}"
            ),
            Self::Upstream(message) => write!(f, "SciRust rejected request: {message}"),
        }
    }
}

impl Error for NoiseInputError {}

/// Parameters for a deterministic zero-mean Gaussian white-noise realization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaussianNoise {
    /// Standard deviation of the generated process.
    sigma: f64,
    /// Explicit reproducibility seed.
    seed: u64,
}

impl GaussianNoise {
    /// Validate and construct Gaussian white-noise parameters.
    pub fn new(sigma: f64, seed: u64) -> Result<Self, NoiseInputError> {
        let params = Self { sigma, seed };
        params.validate()?;
        Ok(params)
    }

    /// Return the validated standard deviation.
    pub fn sigma(&self) -> f64 {
        self.sigma
    }

    /// Return the explicit reproducibility seed.
    pub fn seed(&self) -> u64 {
        self.seed
    }

    fn validate(&self) -> Result<(), NoiseInputError> {
        if !self.sigma.is_finite() {
            return Err(NoiseInputError::NonFinite("sigma"));
        }
        if self.sigma < 0.0 {
            return Err(NoiseInputError::NonPositive("sigma"));
        }
        Ok(())
    }
}

/// Generate a deterministic Gaussian white-noise realization using SciRust's
/// validated `SplitMix64` + Box-Muller implementation.
pub fn gaussian_white_noise(
    params: GaussianNoise,
    samples: usize,
) -> Result<Vec<f64>, NoiseInputError> {
    params.validate()?;
    validate_sample_count(samples)?;

    let requested_bytes = samples
        .checked_mul(size_of::<f64>())
        .ok_or(NoiseInputError::AllocationFailed {
            resource: "Gaussian output buffer",
            requested_bytes: usize::MAX,
        })?;
    let mut output = Vec::new();
    output
        .try_reserve_exact(samples)
        .map_err(|_| NoiseInputError::AllocationFailed {
            resource: "Gaussian output buffer",
            requested_bytes,
        })?;

    let mut rng = SplitMix64::new(params.seed);
    for sample_index in 0..samples {
        let value = params.sigma * rng.next_gaussian();
        if !value.is_finite() {
            return Err(NoiseInputError::NonFiniteOutput {
                operation: "Gaussian white-noise scaling",
                sample_index,
            });
        }
        output.push(value);
    }
    Ok(output)
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
    spectral_signature_with_limits(
        signal,
        sample_rate_hz,
        segment_len,
        overlap,
        SpectralLimits::default(),
    )
}

/// Characterize a real-valued series under explicit allocation and work limits.
pub fn spectral_signature_with_limits(
    signal: &[f64],
    sample_rate_hz: f64,
    segment_len: usize,
    overlap: usize,
    limits: SpectralLimits,
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
    if segment_len > signal.len() {
        return Err(NoiseInputError::InvalidSpectrumRequest(
            "segment_len must not exceed signal length",
        ));
    }

    admit_spectral_request(signal.len(), segment_len, overlap, limits)?;

    let window = fallible_hanning(segment_len)?;
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

fn admit_spectral_request(
    signal_len: usize,
    segment_len: usize,
    overlap: usize,
    limits: SpectralLimits,
) -> Result<(), NoiseInputError> {
    require_spectral_limit("signal samples", signal_len, limits.max_signal_samples)?;
    require_spectral_limit("segment samples", segment_len, limits.max_segment_samples)?;

    let half_spectrum_len = segment_len
        .checked_div(2)
        .and_then(|value| value.checked_add(1))
        .ok_or(NoiseInputError::SpectralBudgetOverflow("working bytes"))?;
    let real_buffers = segment_len
        .checked_mul(2)
        .and_then(|value| value.checked_add(half_spectrum_len.checked_mul(2)?))
        .and_then(|value| value.checked_mul(size_of::<f64>()))
        .ok_or(NoiseInputError::SpectralBudgetOverflow("working bytes"))?;
    let fft_buffer = segment_len
        .checked_mul(size_of::<f64>().checked_mul(2).expect("two f64 values fit"))
        .ok_or(NoiseInputError::SpectralBudgetOverflow("working bytes"))?;
    let working_bytes = real_buffers
        .checked_add(fft_buffer)
        .ok_or(NoiseInputError::SpectralBudgetOverflow("working bytes"))?;
    require_spectral_limit("working bytes", working_bytes, limits.max_working_bytes)?;

    let hop = segment_len - overlap;
    let segments = 1 + (signal_len - segment_len) / hop;
    let fft_stages = segment_len.ilog2() as usize;
    let per_segment_work = segment_len
        .checked_mul(
            fft_stages
                .checked_add(3)
                .ok_or(NoiseInputError::SpectralBudgetOverflow("work units"))?,
        )
        .ok_or(NoiseInputError::SpectralBudgetOverflow("work units"))?;
    let work_units = segments
        .checked_mul(per_segment_work)
        .and_then(|value| value.checked_add(signal_len))
        .and_then(|value| value.checked_add(segment_len))
        .ok_or(NoiseInputError::SpectralBudgetOverflow("work units"))?;
    require_spectral_limit("work units", work_units, limits.max_work_units)
}

fn require_spectral_limit(
    resource: &'static str,
    requested: usize,
    maximum: usize,
) -> Result<(), NoiseInputError> {
    if requested > maximum {
        return Err(NoiseInputError::SpectralLimitExceeded {
            resource,
            requested,
            maximum,
        });
    }
    Ok(())
}

fn fallible_hanning(segment_len: usize) -> Result<Vec<f64>, NoiseInputError> {
    let requested_bytes = segment_len
        .checked_mul(size_of::<f64>())
        .ok_or(NoiseInputError::SpectralBudgetOverflow("window bytes"))?;
    let mut window = Vec::new();
    window
        .try_reserve_exact(segment_len)
        .map_err(|_| NoiseInputError::AllocationFailed {
            resource: "Hann window",
            requested_bytes,
        })?;

    let denominator = (segment_len - 1) as f64;
    window.extend(
        (0..segment_len).map(|index| 0.5 * (1.0 - (2.0 * PI * index as f64 / denominator).cos())),
    );
    Ok(window)
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
        assert_eq!(params.sigma(), 0.0);
        assert_eq!(params.seed(), 7);
        let samples = gaussian_white_noise(params, 64).unwrap();
        assert!(samples.iter().all(|&x| x == 0.0));
    }

    #[test]
    fn gaussian_boundary_revalidates_literal_parameters() {
        let nan = GaussianNoise {
            sigma: f64::NAN,
            seed: 1,
        };
        assert_eq!(
            gaussian_white_noise(nan, 1),
            Err(NoiseInputError::NonFinite("sigma"))
        );

        let negative = GaussianNoise {
            sigma: -1.0,
            seed: 1,
        };
        assert_eq!(
            gaussian_white_noise(negative, 1),
            Err(NoiseInputError::NonPositive("sigma"))
        );

        let infinite = GaussianNoise {
            sigma: f64::INFINITY,
            seed: 1,
        };
        assert_eq!(
            gaussian_white_noise(infinite, 1),
            Err(NoiseInputError::NonFinite("sigma"))
        );
    }

    #[test]
    fn gaussian_boundary_rejects_non_finite_generated_output() {
        let params = GaussianNoise {
            sigma: f64::MAX,
            seed: 42,
        };
        assert!(matches!(
            gaussian_white_noise(params, 256),
            Err(NoiseInputError::NonFiniteOutput {
                operation: "Gaussian white-noise scaling",
                ..
            })
        ));
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

    #[test]
    fn spectral_signature_rejects_segment_longer_than_signal() {
        let signal = [0.0; 8];
        assert_eq!(
            spectral_signature(&signal, 1.0, 16, 0),
            Err(NoiseInputError::InvalidSpectrumRequest(
                "segment_len must not exceed signal length"
            ))
        );
    }

    #[test]
    fn spectral_limits_reject_segment_before_window_allocation() {
        let signal = [0.0; 256];
        let limits = SpectralLimits {
            max_segment_samples: 128,
            ..SpectralLimits::default()
        };
        assert_eq!(
            spectral_signature_with_limits(&signal, 1.0, 256, 0, limits),
            Err(NoiseInputError::SpectralLimitExceeded {
                resource: "segment samples",
                requested: 256,
                maximum: 128,
            })
        );
    }

    #[test]
    fn spectral_limits_reject_signal_read_budget() {
        let signal = [0.0; 256];
        let limits = SpectralLimits {
            max_signal_samples: 128,
            ..SpectralLimits::default()
        };
        assert_eq!(
            spectral_signature_with_limits(&signal, 1.0, 128, 0, limits),
            Err(NoiseInputError::SpectralLimitExceeded {
                resource: "signal samples",
                requested: 256,
                maximum: 128,
            })
        );
    }

    #[test]
    fn spectral_limits_reject_work_amplification_from_overlap() {
        let signal = [0.0; 1024];
        let limits = SpectralLimits {
            max_work_units: 1_000,
            ..SpectralLimits::default()
        };
        assert!(matches!(
            spectral_signature_with_limits(&signal, 1.0, 256, 255, limits),
            Err(NoiseInputError::SpectralLimitExceeded {
                resource: "work units",
                ..
            })
        ));
    }

    #[test]
    fn spectral_limits_reject_working_set_before_window_allocation() {
        let signal = [0.0; 256];
        let limits = SpectralLimits {
            max_working_bytes: 1_024,
            ..SpectralLimits::default()
        };
        assert!(matches!(
            spectral_signature_with_limits(&signal, 1.0, 256, 0, limits),
            Err(NoiseInputError::SpectralLimitExceeded {
                resource: "working bytes",
                ..
            })
        ));
    }

    #[test]
    fn fallible_hann_window_matches_pinned_scirust_method() {
        let expected = scirust_signal::hanning(256);
        let actual = fallible_hanning(256).unwrap();
        assert_eq!(actual, expected);
    }
}
