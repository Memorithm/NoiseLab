//! FitzHugh-Nagumo coherence-resonance calibration.
//!
//! This module studies the canonical excitable system
//!
//! `epsilon * dx/dt = x - x^3/3 - y`
//! `dy/dt = x + a + sigma * xi(t)`
//!
//! with Gaussian white noise on the slow variable. In the deterministic
//! excitable regime `a > 1`, the stable fixed point does not spike. Noise can
//! trigger large excursions; coherence resonance is the non-monotone regime in
//! which the inter-spike intervals become most regular at an intermediate noise
//! amplitude.
//!
//! NoiseLab borrows [`scirust_sim::SplitMix64`] for all stochastic draws. The
//! explicit stochastic step used here is deliberately local to this calibration
//! experiment; it is not presented as a general-purpose SDE solver.

use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

const MAX_STEPS: usize = 10_000_000;

/// Parameters for the canonical excitable FitzHugh-Nagumo model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FitzHughNagumo {
    /// Positive fast/slow time-scale ratio.
    pub epsilon: f64,
    /// Excitability parameter. This calibration requires `a > 1`.
    pub a: f64,
}

impl FitzHughNagumo {
    /// Validate a model in the positive excitable regime.
    pub fn new(epsilon: f64, a: f64) -> Result<Self, FhnError> {
        if !epsilon.is_finite() {
            return Err(FhnError::NonFinite("epsilon"));
        }
        if epsilon <= 0.0 {
            return Err(FhnError::NonPositive("epsilon"));
        }
        if !a.is_finite() {
            return Err(FhnError::NonFinite("a"));
        }
        if a <= 1.0 {
            return Err(FhnError::NotExcitable { a });
        }
        Ok(Self { epsilon, a })
    }

    /// Deterministic stable fixed point `[x*, y*]` for the declared model.
    #[must_use]
    pub fn equilibrium(self) -> [f64; 2] {
        let x = -self.a;
        [x, self.a.powi(3) / 3.0 - self.a]
    }
}

/// Fixed-step and spike-measurement settings for one FHN realization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FhnRun {
    /// Euler step in the model's time unit.
    pub dt: f64,
    /// Number of stochastic steps.
    pub steps: usize,
    /// Leading steps ignored before spike measurement.
    pub burn_in_steps: usize,
    /// Upward activator crossing used to mark one spike.
    pub spike_threshold: f64,
    /// Minimum number of measured spikes required for an inter-spike CV.
    pub min_spikes: usize,
}

impl FhnRun {
    /// Validate run and measurement settings.
    pub fn new(
        dt: f64,
        steps: usize,
        burn_in_steps: usize,
        spike_threshold: f64,
        min_spikes: usize,
    ) -> Result<Self, FhnError> {
        if !dt.is_finite() {
            return Err(FhnError::NonFinite("dt"));
        }
        if dt <= 0.0 {
            return Err(FhnError::NonPositive("dt"));
        }
        if steps == 0 {
            return Err(FhnError::NonPositive("steps"));
        }
        if steps > MAX_STEPS {
            return Err(FhnError::TooManySteps {
                requested: steps,
                maximum: MAX_STEPS,
            });
        }
        if burn_in_steps >= steps {
            return Err(FhnError::InvalidBurnIn {
                burn_in_steps,
                steps,
            });
        }
        if !spike_threshold.is_finite() {
            return Err(FhnError::NonFinite("spike_threshold"));
        }
        if min_spikes < 3 {
            return Err(FhnError::TooFewRequiredSpikes { min_spikes });
        }
        Ok(Self {
            dt,
            steps,
            burn_in_steps,
            spike_threshold,
            min_spikes,
        })
    }
}

/// Spike train extracted directly from one stochastic realization.
#[derive(Debug, Clone, PartialEq)]
pub struct FhnSpikeTrain {
    /// Linearly interpolated upward-crossing times after burn-in.
    pub spike_times: Vec<f64>,
    /// Final `[x, y]` state after all requested steps.
    pub final_state: [f64; 2],
}

/// Inter-spike regularity for one realization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterSpikeStats {
    /// Number of detected spikes.
    pub spikes: usize,
    /// Mean interval between consecutive spikes.
    pub mean_interval: f64,
    /// Sample standard deviation of inter-spike intervals.
    pub sample_stddev: f64,
    /// Coefficient of variation `sample_stddev / mean_interval`.
    pub coefficient_of_variation: f64,
}

/// Aggregated coherence response at one noise amplitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoherenceResponse {
    /// Gaussian white-noise amplitude `sigma` in `dy = ... dt + sigma dW`.
    pub noise_amplitude: f64,
    /// Mean inter-spike coefficient of variation over deterministic seeds.
    pub mean_cv: f64,
    /// Sample standard deviation of replicate CV values.
    pub sample_stddev_cv: f64,
    /// Mean inter-spike interval averaged over replicates.
    pub mean_interval: f64,
    /// Number of deterministic replicate seeds.
    pub replicates: usize,
}

/// Conservative evidence for an interior minimum of inter-spike CV.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoherenceMinimum {
    /// Noise amplitude of the best sampled interior mean CV.
    pub noise_amplitude: f64,
    /// Mean CV at the selected sample.
    pub mean_cv: f64,
    /// Standard error of the selected mean CV.
    pub standard_error: f64,
    /// Upper penalized CV `mean + w*SE` at the selected interior minimum.
    pub conservative_minimum: f64,
    /// Smaller lower penalized CV among the two scan edges.
    pub conservative_edge: f64,
    /// `conservative_edge - conservative_minimum`; positive for returned evidence.
    pub conservative_separation: f64,
    /// Uncertainty multiplier `w` used by the heuristic.
    pub uncertainty_weight: f64,
}

/// Result of a complete coherence-resonance sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct CoherenceCalibration {
    /// Responses in the same order as the requested noise amplitudes.
    pub responses: Vec<CoherenceResponse>,
    /// Strongest uncertainty-separated sampled interior minimum, if present.
    pub interior_minimum: Option<CoherenceMinimum>,
}

/// Validation or simulation failure in the FHN calibration.
#[derive(Debug, Clone, PartialEq)]
pub enum FhnError {
    /// A named scalar is NaN or infinite.
    NonFinite(&'static str),
    /// A named positive scalar is zero or negative.
    NonPositive(&'static str),
    /// The calibration requires the deterministic excitable regime `a > 1`.
    NotExcitable { a: f64 },
    /// Requested stochastic step budget exceeds the hard safety limit.
    TooManySteps { requested: usize, maximum: usize },
    /// Burn-in consumes the whole requested run.
    InvalidBurnIn { burn_in_steps: usize, steps: usize },
    /// A CV requires at least three spikes.
    TooFewRequiredSpikes { min_spikes: usize },
    /// At least one deterministic replicate seed is required.
    NoSeeds,
    /// A noise sweep requires at least three coordinates.
    SweepTooShort,
    /// Noise amplitudes must be finite, non-negative and strictly increasing.
    InvalidNoiseAmplitude { index: usize },
    /// State integration left the finite floating-point domain.
    NonFiniteState { step: usize },
    /// A standalone spike train was malformed or too short for CV measurement.
    InvalidSpikeTrain {
        /// Number of spikes actually supplied.
        observed: usize,
        /// Minimum number of spikes requested.
        required: usize,
    },
    /// A realization in a calibrated sweep did not produce enough spikes.
    InsufficientSpikes {
        noise_amplitude: f64,
        seed: u64,
        observed: usize,
        required: usize,
    },
    /// The uncertainty multiplier must be finite and non-negative.
    InvalidUncertaintyWeight,
}

impl Display for FhnError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(formatter, "{name} must be finite"),
            Self::NonPositive(name) => write!(formatter, "{name} must be strictly positive"),
            Self::NotExcitable { a } => write!(
                formatter,
                "FitzHugh-Nagumo coherence calibration requires a > 1, got {a}"
            ),
            Self::TooManySteps { requested, maximum } => write!(
                formatter,
                "requested {requested} FHN steps exceeds safety maximum {maximum}"
            ),
            Self::InvalidBurnIn {
                burn_in_steps,
                steps,
            } => write!(
                formatter,
                "burn_in_steps {burn_in_steps} must be smaller than steps {steps}"
            ),
            Self::TooFewRequiredSpikes { min_spikes } => write!(
                formatter,
                "min_spikes must be at least 3 to define interval variability, got {min_spikes}"
            ),
            Self::NoSeeds => formatter.write_str("at least one deterministic seed is required"),
            Self::SweepTooShort => formatter.write_str("noise sweep needs at least three points"),
            Self::InvalidNoiseAmplitude { index } => write!(
                formatter,
                "noise amplitude at index {index} must be finite, non-negative and strictly increasing"
            ),
            Self::NonFiniteState { step } => {
                write!(formatter, "FHN state became non-finite at step {step}")
            }
            Self::InvalidSpikeTrain { observed, required } => write!(
                formatter,
                "spike train has {observed} usable spikes, fewer than required {required}, or contains non-increasing/non-finite times"
            ),
            Self::InsufficientSpikes {
                noise_amplitude,
                seed,
                observed,
                required,
            } => write!(
                formatter,
                "noise amplitude {noise_amplitude} with seed {seed} produced {observed} spikes, fewer than required {required}"
            ),
            Self::InvalidUncertaintyWeight => {
                formatter.write_str("uncertainty weight must be finite and non-negative")
            }
        }
    }
}

impl Error for FhnError {}

/// Simulate one FHN realization and extract post-burn-in spikes.
///
/// The discretization is explicit Euler for the fast deterministic equation and
/// Euler-Maruyama for additive noise on the slow variable. Gaussian draws come
/// from SciRust's deterministic [`SplitMix64`].
pub fn simulate_fhn_spikes(
    model: FitzHughNagumo,
    noise_amplitude: f64,
    run: FhnRun,
    seed: u64,
) -> Result<FhnSpikeTrain, FhnError> {
    if !noise_amplitude.is_finite() {
        return Err(FhnError::NonFinite("noise_amplitude"));
    }
    if noise_amplitude < 0.0 {
        return Err(FhnError::InvalidNoiseAmplitude { index: 0 });
    }

    let mut rng = SplitMix64::new(seed);
    let [mut x, mut y] = model.equilibrium();
    let stochastic_scale = noise_amplitude * run.dt.sqrt();
    let mut spike_times = Vec::new();

    for step in 0..run.steps {
        let old_x = x;
        let old_y = y;
        let dx = (old_x - old_x * old_x * old_x / 3.0 - old_y) / model.epsilon;
        let dy = old_x + model.a;
        x = old_x + run.dt * dx;
        y = old_y + run.dt * dy + stochastic_scale * rng.next_gaussian();
        if !x.is_finite() || !y.is_finite() {
            return Err(FhnError::NonFiniteState { step });
        }

        if step >= run.burn_in_steps && old_x < run.spike_threshold && x >= run.spike_threshold {
            let fraction = (run.spike_threshold - old_x) / (x - old_x);
            let crossing_time = (step as f64 + fraction) * run.dt;
            spike_times.push(crossing_time);
        }
    }

    Ok(FhnSpikeTrain {
        spike_times,
        final_state: [x, y],
    })
}

/// Compute inter-spike statistics for one measured spike train.
pub fn inter_spike_stats(
    train: &FhnSpikeTrain,
    min_spikes: usize,
) -> Result<InterSpikeStats, FhnError> {
    if min_spikes < 3 {
        return Err(FhnError::TooFewRequiredSpikes { min_spikes });
    }
    let valid_times = train.spike_times.iter().all(|time| time.is_finite())
        && train.spike_times.windows(2).all(|pair| pair[1] > pair[0]);
    if train.spike_times.len() < min_spikes || !valid_times {
        return Err(FhnError::InvalidSpikeTrain {
            observed: train.spike_times.len(),
            required: min_spikes,
        });
    }

    let intervals: Vec<f64> = train
        .spike_times
        .windows(2)
        .map(|pair| pair[1] - pair[0])
        .collect();
    let mean_interval = intervals.iter().sum::<f64>() / intervals.len() as f64;
    let sum_sq = intervals
        .iter()
        .map(|interval| {
            let delta = *interval - mean_interval;
            delta * delta
        })
        .sum::<f64>();
    let sample_stddev = (sum_sq / (intervals.len() - 1) as f64).sqrt();
    Ok(InterSpikeStats {
        spikes: train.spike_times.len(),
        mean_interval,
        sample_stddev,
        coefficient_of_variation: sample_stddev / mean_interval,
    })
}

/// Sweep noise amplitude and seek an uncertainty-separated interior CV minimum.
///
/// Lower inter-spike CV means more regular noise-induced excursions. The
/// returned minimum is evidence for a sampled non-monotone regularity optimum;
/// it does not by itself prove a universal coherence-resonance law.
pub fn calibrate_coherence_resonance(
    model: FitzHughNagumo,
    run: FhnRun,
    noise_amplitudes: &[f64],
    seeds: &[u64],
    uncertainty_weight: f64,
) -> Result<CoherenceCalibration, FhnError> {
    if noise_amplitudes.len() < 3 {
        return Err(FhnError::SweepTooShort);
    }
    if seeds.is_empty() {
        return Err(FhnError::NoSeeds);
    }
    if !uncertainty_weight.is_finite() || uncertainty_weight < 0.0 {
        return Err(FhnError::InvalidUncertaintyWeight);
    }
    for (index, &noise) in noise_amplitudes.iter().enumerate() {
        if !noise.is_finite() || noise < 0.0 || (index > 0 && noise <= noise_amplitudes[index - 1])
        {
            return Err(FhnError::InvalidNoiseAmplitude { index });
        }
    }

    let mut responses = Vec::with_capacity(noise_amplitudes.len());
    for &noise_amplitude in noise_amplitudes {
        let mut cvs = Vec::with_capacity(seeds.len());
        let mut intervals = Vec::with_capacity(seeds.len());
        for &seed in seeds {
            let train = simulate_fhn_spikes(model, noise_amplitude, run, seed)?;
            if train.spike_times.len() < run.min_spikes {
                return Err(FhnError::InsufficientSpikes {
                    noise_amplitude,
                    seed,
                    observed: train.spike_times.len(),
                    required: run.min_spikes,
                });
            }
            let stats = inter_spike_stats(&train, run.min_spikes)?;
            cvs.push(stats.coefficient_of_variation);
            intervals.push(stats.mean_interval);
        }
        let mean_cv = mean(&cvs);
        responses.push(CoherenceResponse {
            noise_amplitude,
            mean_cv,
            sample_stddev_cv: sample_stddev(&cvs, mean_cv),
            mean_interval: mean(&intervals),
            replicates: seeds.len(),
        });
    }

    let interior_minimum = detect_coherence_minimum(&responses, uncertainty_weight);
    Ok(CoherenceCalibration {
        responses,
        interior_minimum,
    })
}

fn detect_coherence_minimum(
    responses: &[CoherenceResponse],
    uncertainty_weight: f64,
) -> Option<CoherenceMinimum> {
    let (index, minimum) = responses[1..responses.len() - 1]
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| left.mean_cv.total_cmp(&right.mean_cv))
        .map(|(offset, response)| (offset + 1, response))?;
    let minimum_se = standard_error(minimum);
    let conservative_minimum = minimum.mean_cv + uncertainty_weight * minimum_se;
    let left = &responses[0];
    let right = &responses[responses.len() - 1];
    let conservative_left = left.mean_cv - uncertainty_weight * standard_error(left);
    let conservative_right = right.mean_cv - uncertainty_weight * standard_error(right);
    let conservative_edge = conservative_left.min(conservative_right);
    let conservative_separation = conservative_edge - conservative_minimum;

    let strict_local_minimum = minimum.mean_cv < responses[index - 1].mean_cv
        && minimum.mean_cv < responses[index + 1].mean_cv;
    if !strict_local_minimum || conservative_separation <= 0.0 {
        return None;
    }

    Some(CoherenceMinimum {
        noise_amplitude: minimum.noise_amplitude,
        mean_cv: minimum.mean_cv,
        standard_error: minimum_se,
        conservative_minimum,
        conservative_edge,
        conservative_separation,
        uncertainty_weight,
    })
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn sample_stddev(values: &[f64], sample_mean: f64) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let sum_sq = values
        .iter()
        .map(|value| {
            let delta = *value - sample_mean;
            delta * delta
        })
        .sum::<f64>();
    (sum_sq / (values.len() - 1) as f64).sqrt()
}

fn standard_error(response: &CoherenceResponse) -> f64 {
    response.sample_stddev_cv / (response.replicates as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model() -> FitzHughNagumo {
        FitzHughNagumo::new(0.01, 1.05).unwrap()
    }

    fn run() -> FhnRun {
        FhnRun::new(0.001, 80_000, 10_000, 0.0, 5).unwrap()
    }

    #[test]
    fn deterministic_equilibrium_matches_closed_form() {
        let model = model();
        let [x, y] = model.equilibrium();
        assert!((x + 1.05).abs() < 1e-15);
        assert!((x - x.powi(3) / 3.0 - y).abs() < 1e-15);
        assert!((x + model.a).abs() < 1e-15);
    }

    #[test]
    fn same_seed_reproduces_spike_train_exactly() {
        let a = simulate_fhn_spikes(model(), 0.075, run(), 42).unwrap();
        let b = simulate_fhn_spikes(model(), 0.075, run(), 42).unwrap();
        assert_eq!(a, b);
        assert!(a.spike_times.len() >= run().min_spikes);
    }

    #[test]
    fn zero_noise_stays_at_the_excitable_fixed_point() {
        let train = simulate_fhn_spikes(model(), 0.0, run(), 7).unwrap();
        assert!(train.spike_times.is_empty());
        let equilibrium = model().equilibrium();
        assert!((train.final_state[0] - equilibrium[0]).abs() < 1e-12);
        assert!((train.final_state[1] - equilibrium[1]).abs() < 1e-12);
    }

    #[test]
    fn noise_sweep_has_an_interior_regular_spiking_optimum() {
        let noise = [0.02, 0.03, 0.05, 0.075, 0.10, 0.20, 0.40, 0.70];
        let seeds = [11, 23, 37, 41, 53, 67, 79, 97];
        let calibration =
            calibrate_coherence_resonance(model(), run(), &noise, &seeds, 2.0).unwrap();
        let minimum = calibration
            .interior_minimum
            .expect("expected an uncertainty-separated intermediate CV minimum");

        assert!(
            (0.03..=0.20).contains(&minimum.noise_amplitude),
            "unexpected coherence optimum at sigma={}",
            minimum.noise_amplitude
        );
        assert!(minimum.conservative_separation > 0.0);
        assert!(
            minimum.mean_cv < calibration.responses[0].mean_cv,
            "intermediate noise must be more regular than low-noise activation"
        );
        assert!(
            minimum.mean_cv < calibration.responses.last().unwrap().mean_cv,
            "intermediate noise must be more regular than strong-noise spiking"
        );
    }

    #[test]
    fn malformed_calibrations_fail_closed() {
        assert!(FitzHughNagumo::new(0.01, 1.0).is_err());
        assert!(FhnRun::new(0.0, 100, 10, 0.0, 5).is_err());
        assert!(calibrate_coherence_resonance(model(), run(), &[0.1, 0.2], &[1], 1.0).is_err());
        assert!(
            calibrate_coherence_resonance(model(), run(), &[0.1, 0.1, 0.2], &[1], 1.0).is_err()
        );
        let malformed = FhnSpikeTrain {
            spike_times: vec![1.0, 0.5, 2.0],
            final_state: [0.0, 0.0],
        };
        assert!(inter_spike_stats(&malformed, 3).is_err());
    }
}
