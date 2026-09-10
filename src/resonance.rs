//! Resonance and system-dependent operating-point search.
//!
//! NoiseLab keeps an operating-point optimum, an interior response peak, and
//! a model-derived resonance condition as distinct evidence classes. None is
//! promoted to a universal "perfect moment" law.

use std::f64::consts::FRAC_1_SQRT_2;

/// A candidate perturbation operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OperatingPoint {
    /// Perturbation onset time in the target system's time unit.
    pub onset_time: f64,
    /// Non-negative perturbation amplitude or intensity parameter.
    pub amplitude: f64,
    /// Non-negative carrier/forcing frequency in hertz.
    pub frequency_hz: f64,
    /// Perturbation phase in radians.
    pub phase_rad: f64,
}

impl OperatingPoint {
    /// Construct a finite operating point.
    pub fn new(
        onset_time: f64,
        amplitude: f64,
        frequency_hz: f64,
        phase_rad: f64,
    ) -> Result<Self, ResonanceInputError> {
        for (name, value) in [
            ("onset_time", onset_time),
            ("amplitude", amplitude),
            ("frequency_hz", frequency_hz),
            ("phase_rad", phase_rad),
        ] {
            if !value.is_finite() {
                return Err(ResonanceInputError::NonFinite(name));
            }
        }
        if amplitude < 0.0 {
            return Err(ResonanceInputError::Negative("amplitude"));
        }
        if frequency_hz < 0.0 {
            return Err(ResonanceInputError::Negative("frequency_hz"));
        }
        Ok(Self {
            onset_time,
            amplitude,
            frequency_hz,
            phase_rad,
        })
    }
}

/// Validation failures for resonance formulas and searches.
#[derive(Debug, Clone, PartialEq)]
pub enum ResonanceInputError {
    /// A named scalar was NaN or infinite.
    NonFinite(&'static str),
    /// A named scalar that must be non-negative was negative.
    Negative(&'static str),
    /// A named scalar that must be strictly positive was zero or negative.
    NonPositive(&'static str),
    /// A search did not contain enough paired replicates.
    TooFewReplicates { supplied: usize, required: usize },
    /// A sweep contained fewer than three samples.
    SweepTooShort,
    /// Two sweep coordinates were identical.
    DuplicateSweepCoordinate,
}

/// Configuration for paired operating-point evaluation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SearchConfig {
    /// Minimum number of paired seeds required before ranking candidates.
    pub min_replicates: usize,
    /// Penalty multiplier applied to the standard error when ranking.
    ///
    /// The score is `mean_uplift - uncertainty_weight * standard_error`.
    /// This is a ranking heuristic, not automatically a confidence bound.
    pub uncertainty_weight: f64,
    /// If true, every paired replicate must improve on its matched control.
    pub require_positive_worst_case: bool,
}

impl Default for SearchConfig {
    fn default() -> Self {
        Self {
            min_replicates: 3,
            uncertainty_weight: 1.0,
            require_positive_worst_case: false,
        }
    }
}

impl SearchConfig {
    fn validate(self) -> Result<Self, ResonanceInputError> {
        if self.min_replicates == 0 {
            return Err(ResonanceInputError::TooFewReplicates {
                supplied: 0,
                required: 1,
            });
        }
        require_nonnegative("uncertainty_weight", self.uncertainty_weight)?;
        Ok(self)
    }
}

/// Paired estimate for one candidate operating point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CandidateEstimate {
    /// Evaluated operating point.
    pub point: OperatingPoint,
    /// Mean candidate-minus-control utility over paired seeds.
    pub mean_uplift: f64,
    /// Sample standard deviation of paired uplifts.
    pub sample_stddev: f64,
    /// Standard error of the paired-uplift mean.
    pub standard_error: f64,
    /// `mean_uplift - uncertainty_weight * standard_error`.
    pub conservative_score: f64,
    /// Lowest paired uplift observed across the supplied seeds.
    pub worst_uplift: f64,
    /// Number of paired replicates.
    pub replicates: usize,
}

/// Result of a paired operating-point search.
#[derive(Debug, Clone, PartialEq)]
pub struct SearchResult {
    /// Mean utility of the unperturbed/control evaluations.
    pub control_mean: f64,
    /// Estimates in the same order as the candidate input.
    pub candidates: Vec<CandidateEstimate>,
    /// Best candidate that passed the configured positive-uplift gate.
    pub best: Option<CandidateEstimate>,
}

/// Failures that can occur while evaluating a paired operating-point search.
#[derive(Debug)]
pub enum SearchError<E> {
    /// Invalid search configuration.
    Input(ResonanceInputError),
    /// Target evaluator returned its own error.
    Evaluation(E),
    /// Target evaluator returned NaN or infinity.
    NonFiniteUtility {
        /// `None` identifies the control; `Some(i)` identifies candidate `i`.
        candidate_index: Option<usize>,
        /// Seed of the failed paired replicate.
        seed: u64,
    },
}

/// Search candidate operating points using matched control/candidate seeds.
///
/// The evaluator receives `None` for the control and `Some(point)` for an
/// intervention. Reusing each seed for both sides is a common-random-numbers
/// design that pairs stochastic background variation.
pub fn search_operating_points<F, E>(
    candidates: &[OperatingPoint],
    seeds: &[u64],
    config: SearchConfig,
    mut evaluate: F,
) -> Result<SearchResult, SearchError<E>>
where
    F: FnMut(Option<OperatingPoint>, u64) -> Result<f64, E>,
{
    let config = config.validate().map_err(SearchError::Input)?;
    if seeds.len() < config.min_replicates {
        return Err(SearchError::Input(ResonanceInputError::TooFewReplicates {
            supplied: seeds.len(),
            required: config.min_replicates,
        }));
    }

    let mut controls = Vec::with_capacity(seeds.len());
    for &seed in seeds {
        let score = evaluate(None, seed).map_err(SearchError::Evaluation)?;
        if !score.is_finite() {
            return Err(SearchError::NonFiniteUtility {
                candidate_index: None,
                seed,
            });
        }
        controls.push(score);
    }
    let control_mean = mean(&controls);

    let mut estimates = Vec::with_capacity(candidates.len());
    for (candidate_index, &point) in candidates.iter().enumerate() {
        let mut uplifts = Vec::with_capacity(seeds.len());
        for (&seed, &control) in seeds.iter().zip(&controls) {
            let score = evaluate(Some(point), seed).map_err(SearchError::Evaluation)?;
            if !score.is_finite() {
                return Err(SearchError::NonFiniteUtility {
                    candidate_index: Some(candidate_index),
                    seed,
                });
            }
            uplifts.push(score - control);
        }

        let mean_uplift = mean(&uplifts);
        let sample_stddev = sample_stddev(&uplifts, mean_uplift);
        let standard_error = sample_stddev / (uplifts.len() as f64).sqrt();
        let conservative_score = mean_uplift - config.uncertainty_weight * standard_error;
        let worst_uplift = uplifts.iter().copied().fold(f64::INFINITY, f64::min);
        estimates.push(CandidateEstimate {
            point,
            mean_uplift,
            sample_stddev,
            standard_error,
            conservative_score,
            worst_uplift,
            replicates: uplifts.len(),
        });
    }

    let best = estimates
        .iter()
        .copied()
        .filter(|estimate| estimate.conservative_score > 0.0)
        .filter(|estimate| !config.require_positive_worst_case || estimate.worst_uplift > 0.0)
        .max_by(|a, b| a.conservative_score.total_cmp(&b.conservative_score));

    Ok(SearchResult {
        control_mean,
        candidates: estimates,
        best,
    })
}

/// One sample from a scalar response sweep.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweepSample {
    /// Swept parameter coordinate, such as noise intensity or frequency.
    pub coordinate: f64,
    /// Declared response metric at this coordinate.
    pub response: f64,
}

/// A detected interior local response peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InteriorPeak {
    /// Coordinate of the peak sample.
    pub coordinate: f64,
    /// Response at the peak sample.
    pub response: f64,
    /// Peak height above the larger immediate neighbor.
    pub local_prominence: f64,
}

/// Detect the strongest sampled interior local peak of a one-dimensional sweep.
///
/// This supports only the statement that a sampled response is non-monotone
/// with an interior local maximum; it does not identify the mechanism.
pub fn detect_interior_response_peak(
    samples: &[SweepSample],
    minimum_prominence: f64,
) -> Result<Option<InteriorPeak>, ResonanceInputError> {
    if samples.len() < 3 {
        return Err(ResonanceInputError::SweepTooShort);
    }
    require_nonnegative("minimum_prominence", minimum_prominence)?;

    let mut ordered = samples.to_vec();
    for sample in &ordered {
        if !sample.coordinate.is_finite() {
            return Err(ResonanceInputError::NonFinite("sweep coordinate"));
        }
        if !sample.response.is_finite() {
            return Err(ResonanceInputError::NonFinite("sweep response"));
        }
    }
    ordered.sort_by(|a, b| a.coordinate.total_cmp(&b.coordinate));
    if ordered
        .windows(2)
        .any(|pair| pair[0].coordinate == pair[1].coordinate)
    {
        return Err(ResonanceInputError::DuplicateSweepCoordinate);
    }

    let mut best: Option<InteriorPeak> = None;
    for window in ordered.windows(3) {
        let [left, center, right] = [window[0], window[1], window[2]];
        let prominence = center.response - left.response.max(right.response);
        if prominence >= minimum_prominence
            && center.response > left.response
            && center.response > right.response
        {
            let peak = InteriorPeak {
                coordinate: center.coordinate,
                response: center.response,
                local_prominence: prominence,
            };
            if best.is_none_or(|current| peak.response > current.response) {
                best = Some(peak);
            }
        }
    }
    Ok(best)
}

/// Displacement-amplitude resonance frequency for a linearly damped oscillator.
///
/// For `x'' + 2*zeta*omega_n*x' + omega_n^2*x = F*cos(omega*t)/m`, the
/// steady-state displacement amplitude peaks at
/// `omega_r = omega_n * sqrt(1 - 2*zeta^2)` when `zeta < 1/sqrt(2)`.
pub fn linear_displacement_resonance_angular_frequency(
    omega_n: f64,
    damping_ratio: f64,
) -> Result<Option<f64>, ResonanceInputError> {
    require_positive("omega_n", omega_n)?;
    require_nonnegative("damping_ratio", damping_ratio)?;
    if damping_ratio >= FRAC_1_SQRT_2 {
        return Ok(None);
    }
    Ok(Some(
        omega_n * (1.0 - 2.0 * damping_ratio * damping_ratio).sqrt(),
    ))
}

/// Simplified weak-noise Kramers rate `r(D) = r0 * exp(-barrier / D)`.
///
/// `barrier` and `noise_intensity` must use the same energy-like scale. This
/// is a model formula and is not assumed outside its validity regime.
pub fn kramers_escape_rate(
    barrier: f64,
    noise_intensity: f64,
    prefactor_rate: f64,
) -> Result<f64, ResonanceInputError> {
    require_positive("barrier", barrier)?;
    require_positive("noise_intensity", noise_intensity)?;
    require_positive("prefactor_rate", prefactor_rate)?;
    Ok(prefactor_rate * (-barrier / noise_intensity).exp())
}

/// Classical symmetric-bistable stochastic-resonance rate-matching target.
///
/// The adiabatic heuristic asks for approximately one switch per half forcing
/// period, hence `r_target ~= 2 * f_signal`.
pub fn stochastic_resonance_target_rate(
    signal_frequency_hz: f64,
) -> Result<f64, ResonanceInputError> {
    require_positive("signal_frequency_hz", signal_frequency_hz)?;
    Ok(2.0 * signal_frequency_hz)
}

/// Solve the simplified Kramers model for the classical half-period target.
///
/// Returns `None` when the target rate is at or above the model prefactor, for
/// which no finite positive intensity solves the declared equation.
pub fn matched_kramers_noise_intensity(
    barrier: f64,
    prefactor_rate: f64,
    signal_frequency_hz: f64,
) -> Result<Option<f64>, ResonanceInputError> {
    require_positive("barrier", barrier)?;
    require_positive("prefactor_rate", prefactor_rate)?;
    let target_rate = stochastic_resonance_target_rate(signal_frequency_hz)?;
    if target_rate >= prefactor_rate {
        return Ok(None);
    }
    Ok(Some(barrier / (prefactor_rate / target_rate).ln()))
}

fn require_positive(name: &'static str, value: f64) -> Result<(), ResonanceInputError> {
    if !value.is_finite() {
        return Err(ResonanceInputError::NonFinite(name));
    }
    if value <= 0.0 {
        return Err(ResonanceInputError::NonPositive(name));
    }
    Ok(())
}

fn require_nonnegative(name: &'static str, value: f64) -> Result<(), ResonanceInputError> {
    if !value.is_finite() {
        return Err(ResonanceInputError::NonFinite(name));
    }
    if value < 0.0 {
        return Err(ResonanceInputError::Negative(name));
    }
    Ok(())
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
            let delta = value - sample_mean;
            delta * delta
        })
        .sum::<f64>();
    (sum_sq / (values.len() - 1) as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use scirust_sim::apd::{Apd, ApdParams};

    #[test]
    fn linear_resonance_matches_closed_form() {
        let omega_n = 10.0;
        let zeta = 0.2;
        let actual = linear_displacement_resonance_angular_frequency(omega_n, zeta)
            .unwrap()
            .unwrap();
        let expected = omega_n * (1.0_f64 - 2.0 * zeta * zeta).sqrt();
        assert!((actual - expected).abs() < 1e-14);
        assert_eq!(
            linear_displacement_resonance_angular_frequency(omega_n, 0.8).unwrap(),
            None
        );
    }

    #[test]
    fn kramers_matching_inverts_the_declared_rate_model() {
        let barrier = 3.0;
        let prefactor = 100.0;
        let signal_hz = 2.0;
        let intensity = matched_kramers_noise_intensity(barrier, prefactor, signal_hz)
            .unwrap()
            .unwrap();
        let rate = kramers_escape_rate(barrier, intensity, prefactor).unwrap();
        assert!((rate - 2.0 * signal_hz).abs() < 1e-12);
    }

    #[test]
    fn interior_peak_requires_non_monotone_local_maximum() {
        let samples = [
            SweepSample {
                coordinate: 0.1,
                response: 1.0,
            },
            SweepSample {
                coordinate: 0.5,
                response: 3.0,
            },
            SweepSample {
                coordinate: 1.0,
                response: 2.0,
            },
            SweepSample {
                coordinate: 2.0,
                response: 1.5,
            },
        ];
        let peak = detect_interior_response_peak(&samples, 0.5)
            .unwrap()
            .unwrap();
        assert_eq!(peak.coordinate, 0.5);
        assert_eq!(peak.local_prominence, 1.0);
    }

    #[test]
    fn paired_search_cancels_seed_dependent_background() {
        let candidates =
            [1.0, 2.0, 3.0].map(|amplitude| OperatingPoint::new(0.0, amplitude, 0.0, 0.0).unwrap());
        let seeds = [11, 12, 13, 14, 15];
        let result = search_operating_points(
            &candidates,
            &seeds,
            SearchConfig::default(),
            |point, seed| -> Result<f64, ()> {
                let background = (seed % 7) as f64 * 0.125;
                let uplift = point
                    .map(|candidate| 4.0 - (candidate.amplitude - 2.0).powi(2))
                    .unwrap_or(0.0);
                Ok(background + uplift)
            },
        )
        .unwrap();
        assert_eq!(result.best.unwrap().point.amplitude, 2.0);
        assert!(
            result
                .candidates
                .iter()
                .all(|candidate| candidate.sample_stddev == 0.0)
        );
    }

    #[test]
    fn search_finds_scirust_apd_intermediate_optimal_gain() {
        fn apd_snr(gain: f64) -> f64 {
            Apd::new(ApdParams {
                responsivity: 0.8,
                gain,
                ionization_ratio: 0.02,
                dark_current: 1.0e-10,
                r_load: 1.0e4,
                temperature: 300.0,
                bandwidth: 1.0e9,
                optical_power: 1.0e-9,
            })
            .unwrap()
            .snr()
        }

        let gains = [20.0, 50.0, 100.0, 200.0, 1000.0];
        let candidates = gains.map(|gain| OperatingPoint::new(0.0, gain, 0.0, 0.0).unwrap());
        let result = search_operating_points(
            &candidates,
            &[1, 2, 3],
            SearchConfig::default(),
            |point, _seed| -> Result<f64, ()> {
                Ok(apd_snr(
                    point.map(|candidate| candidate.amplitude).unwrap_or(1.0),
                ))
            },
        )
        .unwrap();
        assert_eq!(result.best.unwrap().point.amplitude, 100.0);
    }
}
