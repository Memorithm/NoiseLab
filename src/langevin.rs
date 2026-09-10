//! Executable bistable Langevin controls for stochastic-resonance research.
//!
//! NoiseLab uses the overdamped quartic double well
//!
//! `dx = (a*x - b*x^3 + A*cos(2*pi*f*t + phi))*dt + sqrt(2*D)*dW`
//!
//! as a deliberately simple calibration target.  Its unforced potential,
//! barrier height and weak-noise Kramers prefactor are analytic, while the
//! periodically forced response must be measured.  The simulator borrows
//! SciRust's deterministic `SplitMix64` Gaussian source; no second RNG is
//! introduced here.

use crate::resonance::{
    detect_interior_response_peak, matched_kramers_noise_intensity, InteriorPeak,
    ResonanceInputError, SweepSample,
};
use scirust_sim::SplitMix64;
use std::error::Error;
use std::f64::consts::{SQRT_2, TAU};
use std::fmt::{Display, Formatter};

const MAX_STEPS: usize = 10_000_000;

/// Symmetric quartic bistable potential with a subthreshold sinusoidal tilt.
///
/// At zero forcing the potential is `U(x) = b*x^4/4 - a*x^2/2`, with minima
/// at `±sqrt(a/b)` and a barrier at zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DoubleWellLangevin {
    /// Positive quadratic-instability coefficient.
    pub a: f64,
    /// Positive quartic-restoring coefficient.
    pub b: f64,
    /// Non-negative sinusoidal tilt amplitude.
    pub forcing_amplitude: f64,
    /// Positive forcing frequency in hertz.
    pub forcing_frequency_hz: f64,
    /// Forcing phase in radians.
    pub phase_rad: f64,
}

impl DoubleWellLangevin {
    /// Validate and construct a double-well experiment model.
    pub fn new(
        a: f64,
        b: f64,
        forcing_amplitude: f64,
        forcing_frequency_hz: f64,
        phase_rad: f64,
    ) -> Result<Self, LangevinError> {
        for (name, value) in [
            ("a", a),
            ("b", b),
            ("forcing_amplitude", forcing_amplitude),
            ("forcing_frequency_hz", forcing_frequency_hz),
            ("phase_rad", phase_rad),
        ] {
            if !value.is_finite() {
                return Err(LangevinError::NonFinite(name));
            }
        }
        if a <= 0.0 {
            return Err(LangevinError::NonPositive("a"));
        }
        if b <= 0.0 {
            return Err(LangevinError::NonPositive("b"));
        }
        if forcing_amplitude < 0.0 {
            return Err(LangevinError::Negative("forcing_amplitude"));
        }
        if forcing_frequency_hz <= 0.0 {
            return Err(LangevinError::NonPositive("forcing_frequency_hz"));
        }
        Ok(Self {
            a,
            b,
            forcing_amplitude,
            forcing_frequency_hz,
            phase_rad,
        })
    }

    /// Positive location of either unforced potential minimum.
    #[must_use]
    pub fn well_location(self) -> f64 {
        (self.a / self.b).sqrt()
    }

    /// Unforced barrier height `Delta U = a^2/(4b)`.
    #[must_use]
    pub fn barrier_height(self) -> f64 {
        self.a * self.a / (4.0 * self.b)
    }

    /// Weak-noise overdamped Kramers prefactor for unit mobility.
    ///
    /// For this potential, `|U''(0)| = a` and
    /// `U''(±sqrt(a/b)) = 2a`, hence
    /// `r0 = sqrt(|U''(0)| U''(xmin))/(2*pi) = a*sqrt(2)/(2*pi)`.
    #[must_use]
    pub fn kramers_prefactor_rate(self) -> f64 {
        self.a * SQRT_2 / TAU
    }

    /// Static tilt magnitude at which a metastable well disappears.
    ///
    /// A forcing amplitude below this value is subthreshold in the static
    /// quartic-potential sense.  Stochastic resonance experiments normally use
    /// a weak/subthreshold coherent input so that noise, rather than the
    /// deterministic drive alone, enables switching.
    #[must_use]
    pub fn critical_static_tilt(self) -> f64 {
        2.0 * self.a.powf(1.5) / (3.0 * (3.0 * self.b).sqrt())
    }

    /// Whether this experiment's forcing is below the static barrier-removal
    /// threshold.
    #[must_use]
    pub fn forcing_is_subthreshold(self) -> bool {
        self.forcing_amplitude < self.critical_static_tilt()
    }

    fn drift(self, t: f64, x: f64) -> f64 {
        self.a * x - self.b * x * x * x
            + self.forcing_amplitude
                * (TAU * self.forcing_frequency_hz * t + self.phase_rad).cos()
    }
}

/// Fixed-step settings for one Euler-Maruyama realization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LangevinRun {
    /// Integration step in seconds.
    pub dt: f64,
    /// Total number of integration steps.
    pub steps: usize,
    /// Number of leading samples excluded from response measurement.
    pub burn_in_steps: usize,
}

impl LangevinRun {
    /// Validate run settings.
    pub fn new(dt: f64, steps: usize, burn_in_steps: usize) -> Result<Self, LangevinError> {
        if !dt.is_finite() {
            return Err(LangevinError::NonFinite("dt"));
        }
        if dt <= 0.0 {
            return Err(LangevinError::NonPositive("dt"));
        }
        if steps == 0 {
            return Err(LangevinError::NonPositive("steps"));
        }
        if steps > MAX_STEPS {
            return Err(LangevinError::TooManySteps {
                requested: steps,
                maximum: MAX_STEPS,
            });
        }
        if burn_in_steps >= steps {
            return Err(LangevinError::InvalidBurnIn {
                burn_in_steps,
                steps,
            });
        }
        Ok(Self {
            dt,
            steps,
            burn_in_steps,
        })
    }
}

/// One reproducible scalar Langevin trajectory.
#[derive(Debug, Clone, PartialEq)]
pub struct LangevinTrajectory {
    /// State after every accepted Euler-Maruyama step.
    pub x: Vec<f64>,
    /// Fixed sample interval.
    pub dt: f64,
}

/// Phase-locked response at the declared forcing frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoherentResponse {
    /// Cosine projection of the measured signal.
    pub in_phase: f64,
    /// Sine projection of the measured signal.
    pub quadrature: f64,
    /// Magnitude `sqrt(in_phase^2 + quadrature^2)`.
    pub amplitude: f64,
}

/// Averaged response for one noise intensity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoiseResponse {
    /// Langevin noise intensity `D` in `sqrt(2D) dW`.
    pub noise_intensity: f64,
    /// Mean phase-locked switching amplitude over seeds.
    pub mean_coherent_amplitude: f64,
    /// Sample standard deviation over seeds.
    pub sample_stddev: f64,
    /// Number of deterministic replicates.
    pub replicates: usize,
}

/// Result of a stochastic-resonance calibration sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct KramersCalibration {
    /// Classical weak-noise half-period matching prediction, when finite.
    pub predicted_noise_intensity: Option<f64>,
    /// Measured responses in the input intensity order.
    pub responses: Vec<NoiseResponse>,
    /// Strongest sampled interior response peak, if one exists.
    pub interior_peak: Option<InteriorPeak>,
}

/// Input or simulation failure in the Langevin calibration layer.
#[derive(Debug, Clone, PartialEq)]
pub enum LangevinError {
    /// A named scalar was NaN or infinite.
    NonFinite(&'static str),
    /// A named scalar that must be positive was zero or negative.
    NonPositive(&'static str),
    /// A named scalar that must be non-negative was negative.
    Negative(&'static str),
    /// Step budget exceeded.
    TooManySteps { requested: usize, maximum: usize },
    /// Burn-in consumed the whole trajectory.
    InvalidBurnIn { burn_in_steps: usize, steps: usize },
    /// No deterministic replicate seeds were supplied.
    NoSeeds,
    /// A state left the finite floating-point domain.
    NonFiniteState { step: usize },
    /// Existing resonance helper rejected a derived/input value.
    Resonance(ResonanceInputError),
}

impl Display for LangevinError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite(name) => write!(formatter, "{name} must be finite"),
            Self::NonPositive(name) => write!(formatter, "{name} must be strictly positive"),
            Self::Negative(name) => write!(formatter, "{name} must be non-negative"),
            Self::TooManySteps { requested, maximum } => write!(
                formatter,
                "requested {requested} Langevin steps exceeds safety maximum {maximum}"
            ),
            Self::InvalidBurnIn {
                burn_in_steps,
                steps,
            } => write!(
                formatter,
                "burn_in_steps {burn_in_steps} must be smaller than steps {steps}"
            ),
            Self::NoSeeds => formatter.write_str("at least one deterministic seed is required"),
            Self::NonFiniteState { step } => {
                write!(formatter, "Langevin state became non-finite at step {step}")
            }
            Self::Resonance(error) => write!(formatter, "resonance input error: {error:?}"),
        }
    }
}

impl Error for LangevinError {}

/// Simulate one deterministic realization with Euler-Maruyama.
///
/// The numerical scheme is intentionally local to this concrete experimental
/// model for now; SciRust supplies the RNG.  NoiseLab does not present this as
/// a general-purpose SDE solver.
pub fn simulate_double_well(
    model: DoubleWellLangevin,
    noise_intensity: f64,
    initial_state: f64,
    run: LangevinRun,
    seed: u64,
) -> Result<LangevinTrajectory, LangevinError> {
    if !noise_intensity.is_finite() {
        return Err(LangevinError::NonFinite("noise_intensity"));
    }
    if noise_intensity < 0.0 {
        return Err(LangevinError::Negative("noise_intensity"));
    }
    if !initial_state.is_finite() {
        return Err(LangevinError::NonFinite("initial_state"));
    }

    let mut rng = SplitMix64::new(seed);
    let stochastic_scale = (2.0 * noise_intensity * run.dt).sqrt();
    let mut state = initial_state;
    let mut trajectory = Vec::with_capacity(run.steps);
    for step in 0..run.steps {
        let t = step as f64 * run.dt;
        state += model.drift(t, state) * run.dt + stochastic_scale * rng.next_gaussian();
        if !state.is_finite() {
            return Err(LangevinError::NonFiniteState { step });
        }
        trajectory.push(state);
    }
    Ok(LangevinTrajectory {
        x: trajectory,
        dt: run.dt,
    })
}

/// Measure phase locking of *well occupancy* to the coherent drive.
///
/// Each state is mapped to `-1/+1` before projection.  This intentionally
/// measures noise-assisted switching rather than small intra-well oscillations,
/// which can otherwise look like useful coherent response even when no barrier
/// crossing occurs.
pub fn coherent_switching_response(
    model: DoubleWellLangevin,
    trajectory: &LangevinTrajectory,
    burn_in_steps: usize,
) -> Result<CoherentResponse, LangevinError> {
    if burn_in_steps >= trajectory.x.len() {
        return Err(LangevinError::InvalidBurnIn {
            burn_in_steps,
            steps: trajectory.x.len(),
        });
    }
    if !trajectory.dt.is_finite() || trajectory.dt <= 0.0 {
        return Err(LangevinError::NonPositive("trajectory.dt"));
    }

    let mut cosine = 0.0;
    let mut sine = 0.0;
    let mut count = 0usize;
    for (step, &state) in trajectory.x.iter().enumerate().skip(burn_in_steps) {
        let t = (step + 1) as f64 * trajectory.dt;
        let phase = TAU * model.forcing_frequency_hz * t + model.phase_rad;
        let occupancy = if state >= 0.0 { 1.0 } else { -1.0 };
        cosine += occupancy * phase.cos();
        sine += occupancy * phase.sin();
        count += 1;
    }
    let scale = 2.0 / count as f64;
    let in_phase = scale * cosine;
    let quadrature = scale * sine;
    Ok(CoherentResponse {
        in_phase,
        quadrature,
        amplitude: in_phase.hypot(quadrature),
    })
}

/// Sweep noise intensity and test whether the sampled switching response has an
/// interior non-monotone maximum.
///
/// The returned Kramers prediction is an analytic *control*, not an instruction
/// to force the empirical maximum to that location.  The measured sweep is
/// computed independently using identical seed sets at every intensity.
pub fn calibrate_kramers_matching(
    model: DoubleWellLangevin,
    run: LangevinRun,
    noise_intensities: &[f64],
    seeds: &[u64],
    minimum_prominence: f64,
) -> Result<KramersCalibration, LangevinError> {
    if seeds.is_empty() {
        return Err(LangevinError::NoSeeds);
    }
    let predicted_noise_intensity = matched_kramers_noise_intensity(
        model.barrier_height(),
        model.kramers_prefactor_rate(),
        model.forcing_frequency_hz,
    )
    .map_err(LangevinError::Resonance)?;

    let mut responses = Vec::with_capacity(noise_intensities.len());
    for &noise_intensity in noise_intensities {
        if !noise_intensity.is_finite() {
            return Err(LangevinError::NonFinite("noise_intensity"));
        }
        if noise_intensity < 0.0 {
            return Err(LangevinError::Negative("noise_intensity"));
        }
        let mut amplitudes = Vec::with_capacity(seeds.len());
        for &seed in seeds {
            let trajectory = simulate_double_well(
                model,
                noise_intensity,
                -model.well_location(),
                run,
                seed,
            )?;
            amplitudes.push(
                coherent_switching_response(model, &trajectory, run.burn_in_steps)?.amplitude,
            );
        }
        let mean = amplitudes.iter().sum::<f64>() / amplitudes.len() as f64;
        let sample_stddev = if amplitudes.len() < 2 {
            0.0
        } else {
            let sum_sq = amplitudes
                .iter()
                .map(|value| {
                    let delta = *value - mean;
                    delta * delta
                })
                .sum::<f64>();
            (sum_sq / (amplitudes.len() - 1) as f64).sqrt()
        };
        responses.push(NoiseResponse {
            noise_intensity,
            mean_coherent_amplitude: mean,
            sample_stddev,
            replicates: amplitudes.len(),
        });
    }

    let sweep: Vec<SweepSample> = responses
        .iter()
        .map(|response| SweepSample {
            coordinate: response.noise_intensity,
            response: response.mean_coherent_amplitude,
        })
        .collect();
    let interior_peak = detect_interior_response_peak(&sweep, minimum_prominence)
        .map_err(LangevinError::Resonance)?;

    Ok(KramersCalibration {
        predicted_noise_intensity,
        responses,
        interior_peak,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn model() -> DoubleWellLangevin {
        DoubleWellLangevin::new(1.0, 1.0, 0.20, 0.05, 0.0).unwrap()
    }

    #[test]
    fn quartic_oracles_match_closed_forms() {
        let model = model();
        assert!((model.well_location() - 1.0).abs() < 1e-15);
        assert!((model.barrier_height() - 0.25).abs() < 1e-15);
        assert!((model.kramers_prefactor_rate() - SQRT_2 / TAU).abs() < 1e-15);
        assert!(model.forcing_is_subthreshold());
    }

    #[test]
    fn deterministic_seed_reproduces_the_entire_path() {
        let run = LangevinRun::new(0.02, 10_000, 1_000).unwrap();
        let a = simulate_double_well(model(), 0.3, -1.0, run, 42).unwrap();
        let b = simulate_double_well(model(), 0.3, -1.0, run, 42).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn zero_noise_path_is_seed_independent() {
        let run = LangevinRun::new(0.02, 2_000, 100).unwrap();
        let a = simulate_double_well(model(), 0.0, -1.0, run, 1).unwrap();
        let b = simulate_double_well(model(), 0.0, -1.0, run, 999).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn kramers_control_predicts_finite_noise_for_calibration_model() {
        let model = model();
        let predicted = matched_kramers_noise_intensity(
            model.barrier_height(),
            model.kramers_prefactor_rate(),
            model.forcing_frequency_hz,
        )
        .unwrap()
        .unwrap();
        assert!(predicted > 0.2 && predicted < 0.4, "predicted D={predicted}");
    }

    #[test]
    fn noise_sweep_has_an_interior_phase_locked_switching_peak() {
        // Fifty forcing periods, discarding the first ten.  Fixed seeds make
        // this an executable calibration rather than a flaky Monte Carlo test.
        let model = model();
        let period = 1.0 / model.forcing_frequency_hz;
        let dt = 0.02;
        let steps_per_period = (period / dt) as usize;
        let run = LangevinRun::new(dt, 50 * steps_per_period, 10 * steps_per_period).unwrap();
        let intensities = [0.05, 0.10, 0.20, 0.30, 0.45, 0.70, 1.00];
        let seeds = [11, 23, 37, 41, 53, 67, 79, 97];
        let calibration =
            calibrate_kramers_matching(model, run, &intensities, &seeds, 0.005).unwrap();
        let peak = calibration
            .interior_peak
            .expect("bistable switching response should be non-monotone");
        let predicted = calibration.predicted_noise_intensity.unwrap();

        assert!(
            (0.20..=0.70).contains(&peak.coordinate),
            "unexpected empirical peak at D={} (prediction D={predicted})",
            peak.coordinate
        );
        assert!(
            peak.coordinate / predicted < 2.5 && predicted / peak.coordinate < 2.5,
            "empirical D={} is not in the same broad weak-noise regime as prediction D={predicted}",
            peak.coordinate
        );
    }

    #[test]
    fn invalid_inputs_fail_closed() {
        assert!(DoubleWellLangevin::new(0.0, 1.0, 0.1, 0.05, 0.0).is_err());
        assert!(LangevinRun::new(0.02, 100, 100).is_err());
        let run = LangevinRun::new(0.02, 100, 10).unwrap();
        assert!(simulate_double_well(model(), -0.1, -1.0, run, 1).is_err());
    }
}
