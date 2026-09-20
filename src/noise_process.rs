//! Reproducible scalar noise-process specifications.
//!
//! This module keeps process selection and parameter semantics explicit while
//! delegating random-number generation and the exact OU transition law to the
//! pinned SciRust adapters.  The OU scale exposed here is the stationary
//! standard deviation, not the diffusion coefficient used by the SDE.

use crate::scirust_bridge::{
    gaussian_white_noise, ornstein_uhlenbeck_path, GaussianNoise, NoiseInputError,
};
use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

const INITIAL_STATE_SEED_SALT: u64 = 0x6a09_e667_f3bc_c909;
const PATH_SEED_SALT: u64 = 0xbb67_ae85_84ca_a73b;

/// A stationary Ornstein-Uhlenbeck process specification.
///
/// The generated process has zero mean and stationary standard deviation
/// approximately equal to stationary_stddev.  The first state is sampled from
/// that stationary distribution, so the finite path does not begin at a
/// deterministic zero-state transient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrnsteinUhlenbeckNoise {
    /// Desired stationary standard deviation.
    pub stationary_stddev: f64,
    /// Mean-reversion rate, strictly positive.
    pub theta: f64,
    /// Sample interval, strictly positive.
    pub dt: f64,
    /// Explicit reproducibility seed.
    pub seed: u64,
}

impl OrnsteinUhlenbeckNoise {
    /// Validate and construct a stationary OU specification.
    pub fn new(
        stationary_stddev: f64,
        theta: f64,
        dt: f64,
        seed: u64,
    ) -> Result<Self, NoiseProcessError> {
        require_finite("stationary_stddev", stationary_stddev)?;
        require_finite("theta", theta)?;
        require_finite("dt", dt)?;
        if stationary_stddev < 0.0 {
            return Err(NoiseProcessError::Negative("stationary_stddev"));
        }
        if theta <= 0.0 {
            return Err(NoiseProcessError::NonPositive("theta"));
        }
        if dt <= 0.0 {
            return Err(NoiseProcessError::NonPositive("dt"));
        }
        Ok(Self {
            stationary_stddev,
            theta,
            dt,
            seed,
        })
    }

    /// Generate a deterministic stationary OU realization.
    pub fn generate(self, samples: usize) -> Result<Vec<f64>, NoiseProcessError> {
        if samples == 0 {
            return Ok(Vec::new());
        }
        if self.stationary_stddev == 0.0 {
            return Ok(vec![0.0; samples]);
        }

        let mut initial_rng = SplitMix64::new(self.seed ^ INITIAL_STATE_SEED_SALT);
        let initial_state = self.stationary_stddev * initial_rng.next_gaussian();

        // For dx = theta * (mean - x) dt + sigma dW, the stationary variance
        // is sigma^2 / (2 * theta).  Convert the public stationary scale to
        // the diffusion coefficient expected by SciRust.
        let diffusion_sigma = self.stationary_stddev * (2.0 * self.theta).sqrt();
        ornstein_uhlenbeck_path(
            initial_state,
            self.theta,
            0.0,
            diffusion_sigma,
            self.dt,
            samples,
            self.seed ^ PATH_SEED_SALT,
        )
        .map_err(NoiseProcessError::Upstream)
    }
}

/// A typed process family used by NoiseLab experiments.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NoiseProcessSpec {
    /// Independent zero-mean Gaussian samples.
    GaussianWhite(GaussianNoise),
    /// Zero-mean stationary Ornstein-Uhlenbeck samples.
    OrnsteinUhlenbeck(OrnsteinUhlenbeckNoise),
}

impl NoiseProcessSpec {
    /// Construct a white-Gaussian specification.
    pub fn gaussian_white(sigma: f64, seed: u64) -> Result<Self, NoiseProcessError> {
        GaussianNoise::new(sigma, seed)
            .map(Self::GaussianWhite)
            .map_err(NoiseProcessError::Upstream)
    }

    /// Construct a stationary OU specification.
    pub fn ornstein_uhlenbeck(
        stationary_stddev: f64,
        theta: f64,
        dt: f64,
        seed: u64,
    ) -> Result<Self, NoiseProcessError> {
        OrnsteinUhlenbeckNoise::new(stationary_stddev, theta, dt, seed)
            .map(Self::OrnsteinUhlenbeck)
    }

    /// Generate the selected reproducible process.
    pub fn generate(self, samples: usize) -> Result<Vec<f64>, NoiseProcessError> {
        match self {
            Self::GaussianWhite(params) => {
                gaussian_white_noise(params, samples).map_err(NoiseProcessError::Upstream)
            }
            Self::OrnsteinUhlenbeck(params) => params.generate(samples),
        }
    }
}

/// Validation or upstream failure while generating a typed noise process.
#[derive(Debug, Clone, PartialEq)]
pub enum NoiseProcessError {
    /// A floating-point parameter was NaN or infinite.
    NonFinite(&'static str),
    /// A parameter that may be zero was negative.
    Negative(&'static str),
    /// A strictly positive parameter was zero or negative.
    NonPositive(&'static str),
    /// The pinned SciRust adapter rejected the request.
    Upstream(NoiseInputError),
}

impl Display for NoiseProcessError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(formatter, "{name} must be finite"),
            Self::Negative(name) => write!(formatter, "{name} must be non-negative"),
            Self::NonPositive(name) => write!(formatter, "{name} must be strictly positive"),
            Self::Upstream(error) => write!(formatter, "noise source rejected request: {error}"),
        }
    }
}

impl Error for NoiseProcessError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Upstream(error) => Some(error),
            _ => None,
        }
    }
}

fn require_finite(name: &'static str, value: f64) -> Result<(), NoiseProcessError> {
    if !value.is_finite() {
        return Err(NoiseProcessError::NonFinite(name));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn white_process_is_reproducible_through_typed_spec() {
        let spec = NoiseProcessSpec::gaussian_white(0.75, 42).unwrap();
        assert_eq!(spec.generate(256).unwrap(), spec.generate(256).unwrap());
    }

    #[test]
    fn zero_scale_is_exact_zero_for_both_process_families() {
        let white = NoiseProcessSpec::gaussian_white(0.0, 42).unwrap();
        let ou = NoiseProcessSpec::ornstein_uhlenbeck(0.0, 0.8, 0.05, 99).unwrap();
        let white_values = white.generate(64).unwrap();
        assert!(white_values.iter().all(|&value| value == 0.0));
        assert!(ou.generate(64).unwrap().iter().all(|&value| value == 0.0));
    }

    #[test]
    fn stationary_ou_process_is_reproducible_and_seeded() {
        let first = NoiseProcessSpec::ornstein_uhlenbeck(1.0, 0.8, 0.05, 99)
            .unwrap()
            .generate(128)
            .unwrap();
        let second = NoiseProcessSpec::ornstein_uhlenbeck(1.0, 0.8, 0.05, 99)
            .unwrap()
            .generate(128)
            .unwrap();
        let different = NoiseProcessSpec::ornstein_uhlenbeck(1.0, 0.8, 0.05, 100)
            .unwrap()
            .generate(128)
            .unwrap();
        assert_eq!(first.len(), 128);
        assert_eq!(first, second);
        assert_ne!(first, different);
    }

    #[test]
    fn invalid_ou_parameters_fail_closed() {
        assert!(matches!(
            OrnsteinUhlenbeckNoise::new(-1.0, 0.8, 0.05, 1),
            Err(NoiseProcessError::Negative("stationary_stddev"))
        ));
        assert!(matches!(
            OrnsteinUhlenbeckNoise::new(1.0, 0.0, 0.05, 1),
            Err(NoiseProcessError::NonPositive("theta"))
        ));
        assert!(matches!(
            OrnsteinUhlenbeckNoise::new(1.0, 0.8, f64::NAN, 1),
            Err(NoiseProcessError::NonFinite("dt"))
        ));
    }

    #[test]
    fn empty_generation_is_allocation_free() {
        let spec = NoiseProcessSpec::ornstein_uhlenbeck(1.0, 0.8, 0.05, 99).unwrap();
        assert!(spec.generate(0).unwrap().is_empty());
    }
}
