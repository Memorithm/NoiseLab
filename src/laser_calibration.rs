//! Blind calibration of NoiseLab's period measurement on the SciRust
//! semiconductor-laser rate equations.
//!
//! The measurement path never receives the analytic relaxation frequency. It
//! perturbs the above-threshold operating point, integrates the actual SciRust
//! model, and estimates a period from the resulting photon-density trajectory.
//! Only after that measurement is complete is the independent linearized
//! prediction computed for calibration.

use scirust_sim::laser::{LaserParams, SemiconductorLaser};
use scirust_sim::simulate;
use std::error::Error;
use std::f64::consts::TAU;
use std::fmt::{Display, Formatter};

/// Numerical settings used to excite and observe relaxation ringing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaserRingMeasurement {
    /// Total simulation duration in seconds.
    pub duration_s: f64,
    /// Fixed RK4 step in seconds.
    pub step_s: f64,
    /// Fractional downward perturbation of the steady-state photon density.
    /// Must lie strictly between zero and one.
    pub photon_kick_fraction: f64,
}

impl Default for LaserRingMeasurement {
    fn default() -> Self {
        Self {
            duration_s: 8.0,
            step_s: 5.0e-4,
            photon_kick_fraction: 0.30,
        }
    }
}

impl LaserRingMeasurement {
    fn validate(self) -> Result<Self, LaserCalibrationError> {
        for (name, value) in [
            ("duration_s", self.duration_s),
            ("step_s", self.step_s),
            ("photon_kick_fraction", self.photon_kick_fraction),
        ] {
            if !value.is_finite() {
                return Err(LaserCalibrationError::InvalidSetting(name));
            }
        }
        if self.duration_s <= 0.0 {
            return Err(LaserCalibrationError::InvalidSetting("duration_s"));
        }
        if self.step_s <= 0.0 || self.step_s >= self.duration_s {
            return Err(LaserCalibrationError::InvalidSetting("step_s"));
        }
        if !(0.0 < self.photon_kick_fraction && self.photon_kick_fraction < 1.0) {
            return Err(LaserCalibrationError::InvalidSetting(
                "photon_kick_fraction",
            ));
        }
        Ok(self)
    }
}

/// Evidence returned by the blind laser relaxation calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaserRelaxationCalibration {
    /// Frequency inferred only from the simulated photon-density trajectory.
    pub measured_frequency_hz: f64,
    /// Period inferred from that same trajectory.
    pub measured_period_s: f64,
    /// SciRust model's undamped natural relaxation frequency.
    pub undamped_frequency_hz: f64,
    /// Damped small-signal frequency obtained independently from the
    /// linearized Jacobian about the operating point.
    pub predicted_damped_frequency_hz: f64,
    /// Exponential amplitude-decay rate of the linearized ringing.
    pub decay_rate_per_s: f64,
    /// Relative disagreement between measured and predicted damped frequency.
    pub relative_error: f64,
    /// Steady-state photon density about which the trajectory rings.
    pub steady_state_photon_density: f64,
}

/// Failures of the laser calibration protocol.
#[derive(Debug, Clone, PartialEq)]
pub enum LaserCalibrationError {
    /// A measurement setting is non-finite or outside its declared domain.
    InvalidSetting(&'static str),
    /// SciRust rejected the laser parameters.
    InvalidModel(String),
    /// The current analytic operating-point oracle is the beta=0 limit only.
    AnalyticOracleRequiresZeroBeta,
    /// The configured pump is not above the beta=0 threshold.
    NoAboveThresholdOperatingPoint,
    /// The linearized eigenvalues are real rather than oscillatory.
    NonOscillatoryLinearization,
    /// SciRust failed to integrate the requested trajectory.
    Simulation(String),
    /// The expected photon-density column was absent from the trajectory.
    MissingPhotonSeries,
    /// The late trajectory did not contain enough same-direction crossings to
    /// estimate more than one full period.
    InsufficientLateCrossings,
}

impl Display for LaserCalibrationError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSetting(name) => {
                write!(formatter, "invalid laser calibration setting: {name}")
            }
            Self::InvalidModel(message) => {
                write!(formatter, "invalid SciRust laser model: {message}")
            }
            Self::AnalyticOracleRequiresZeroBeta => formatter.write_str(
                "the closed-form laser calibration oracle applies only when beta = 0",
            ),
            Self::NoAboveThresholdOperatingPoint => formatter.write_str(
                "laser calibration requires an above-threshold steady state with positive photon density",
            ),
            Self::NonOscillatoryLinearization => formatter.write_str(
                "the configured laser linearization has no oscillatory relaxation mode",
            ),
            Self::Simulation(message) => write!(formatter, "SciRust simulation failed: {message}"),
            Self::MissingPhotonSeries => {
                formatter.write_str("SciRust trajectory has no photon-density series")
            }
            Self::InsufficientLateCrossings => formatter.write_str(
                "late photon-density trajectory has too few upward crossings for a period estimate",
            ),
        }
    }
}

impl Error for LaserCalibrationError {}

/// Run a blind relaxation-frequency calibration against SciRust's laser model.
///
/// The measured frequency is obtained from simulation before the analytic
/// frequency is evaluated. This prevents the expected answer from influencing
/// the period estimator while still allowing the result to be scored against a
/// known physical control.
pub fn calibrate_laser_relaxation(
    params: LaserParams,
    settings: LaserRingMeasurement,
) -> Result<LaserRelaxationCalibration, LaserCalibrationError> {
    let settings = settings.validate()?;
    if params.beta != 0.0 {
        return Err(LaserCalibrationError::AnalyticOracleRequiresZeroBeta);
    }

    let g0 = params.g0;
    let tau_n = params.tau_n;
    let model = SemiconductorLaser::new(params)
        .map_err(|error| LaserCalibrationError::InvalidModel(error.to_string()))?;
    let steady_photon = model.steady_state_photon_density();
    if !steady_photon.is_finite() || steady_photon <= 0.0 {
        return Err(LaserCalibrationError::NoAboveThresholdOperatingPoint);
    }

    let measured_period_s = measure_relaxation_period(&model, steady_photon, settings)?;
    let measured_frequency_hz = 1.0 / measured_period_s;

    // This analytic path is intentionally evaluated only after the trajectory
    // has been measured. It follows the independent Jacobian calculation used
    // by SciRust Studio: omega_n^2 = g0*s_ss/tau_p and
    // gamma = (1/tau_n + g0*s_ss)/2.
    let undamped_frequency_hz = model.relaxation_frequency();
    let omega_n = TAU * undamped_frequency_hz;
    let decay_rate_per_s = 0.5 * (1.0 / tau_n + g0 * steady_photon);
    let omega_d_squared = omega_n * omega_n - decay_rate_per_s * decay_rate_per_s;
    if omega_d_squared <= 0.0 {
        return Err(LaserCalibrationError::NonOscillatoryLinearization);
    }
    let predicted_damped_frequency_hz = omega_d_squared.sqrt() / TAU;
    let relative_error = (measured_frequency_hz - predicted_damped_frequency_hz).abs()
        / predicted_damped_frequency_hz;

    Ok(LaserRelaxationCalibration {
        measured_frequency_hz,
        measured_period_s,
        undamped_frequency_hz,
        predicted_damped_frequency_hz,
        decay_rate_per_s,
        relative_error,
        steady_state_photon_density: steady_photon,
    })
}

fn measure_relaxation_period(
    model: &SemiconductorLaser,
    steady_photon: f64,
    settings: LaserRingMeasurement,
) -> Result<f64, LaserCalibrationError> {
    let steady_carrier = model.steady_state_carrier_density();
    let initial_photon = steady_photon * (1.0 - settings.photon_kick_fraction);
    let initial_state = model.initial_state(steady_carrier, initial_photon);
    let trajectory = simulate(
        model,
        &initial_state,
        0.0,
        settings.duration_s,
        settings.step_s,
    )
    .map_err(|error| LaserCalibrationError::Simulation(error.to_string()))?;
    let photon = trajectory
        .column(1)
        .ok_or(LaserCalibrationError::MissingPhotonSeries)?;
    period_from_second_half(&trajectory.t, &photon, steady_photon)
        .ok_or(LaserCalibrationError::InsufficientLateCrossings)
}

/// Private NoiseLab copy of SciRust Studio's trajectory-period estimator.
///
/// Provenance: `scirust-studio-runtime/src/measure.rs` at SciRust revision
/// `f57d598bf03e5dfb16ec6423e4e43105a77540d3`. The upstream function is
/// `pub(crate)`, so NoiseLab keeps this control local rather than pinning an
/// unmerged SciRust branch. The algorithm uses linearly interpolated upward
/// crossings from the second half of the trajectory, exactly matching the
/// qualified SciRust Studio control used by the laser adapter.
fn period_from_second_half(t: &[f64], values: &[f64], level: f64) -> Option<f64> {
    let (first, last) = (*t.first()?, *t.last()?);
    let midpoint = first + (last - first) / 2.0;
    let late: Vec<f64> = upward_crossing_times(t, values, level)
        .into_iter()
        .filter(|crossing| *crossing >= midpoint)
        .collect();
    if late.len() < 3 {
        return None;
    }
    let periods = (late.len() - 1) as f64;
    Some((late[late.len() - 1] - late[0]) / periods)
}

fn upward_crossing_times(t: &[f64], values: &[f64], level: f64) -> Vec<f64> {
    let mut crossings = Vec::new();
    for i in 1..values.len().min(t.len()) {
        let (a, b) = (values[i - 1] - level, values[i] - level);
        if a < 0.0 && b >= 0.0 {
            let span = b - a;
            let fraction = if span != 0.0 { -a / span } else { 0.0 };
            crossings.push(t[i - 1] + fraction * (t[i] - t[i - 1]));
        }
    }
    crossings
}

#[cfg(test)]
mod tests {
    use super::*;

    fn normalized_laser(beta: f64, pump: f64) -> LaserParams {
        LaserParams {
            g0: 1.0,
            n_t: 1.0,
            tau_n: 1.0,
            tau_p: 0.01,
            gamma: 1.0,
            beta,
            pump,
        }
    }

    #[test]
    fn blind_trajectory_measurement_recovers_damped_relaxation_frequency() {
        let calibration = calibrate_laser_relaxation(
            normalized_laser(0.0, 150.0),
            LaserRingMeasurement::default(),
        )
        .unwrap();

        assert!(
            calibration.relative_error < 2.0e-3,
            "measured {} Hz, predicted damped {} Hz, relative error {}",
            calibration.measured_frequency_hz,
            calibration.predicted_damped_frequency_hz,
            calibration.relative_error
        );
        assert!(
            (calibration.measured_frequency_hz - calibration.predicted_damped_frequency_hz).abs()
                < (calibration.measured_frequency_hz - calibration.undamped_frequency_hz).abs(),
            "measurement should resolve the damping correction"
        );
    }

    #[test]
    fn period_estimator_preserves_scirust_known_sine_control() {
        let period = 0.7;
        let omega = TAU / period;
        let step = period / 37.0;
        let t: Vec<f64> = (0..400).map(|i| i as f64 * step).collect();
        let values: Vec<f64> = t.iter().map(|time| (omega * time).sin()).collect();
        let measured = period_from_second_half(&t, &values, 0.0).expect("many crossings");
        assert!((measured - period).abs() / period < 1.0e-4);
    }

    #[test]
    fn below_threshold_configuration_is_rejected_as_a_calibration_target() {
        let error = calibrate_laser_relaxation(
            normalized_laser(0.0, 50.0),
            LaserRingMeasurement::default(),
        )
        .unwrap_err();
        assert_eq!(error, LaserCalibrationError::NoAboveThresholdOperatingPoint);
    }

    #[test]
    fn beta_positive_does_not_reuse_the_beta_zero_oracle() {
        let error = calibrate_laser_relaxation(
            normalized_laser(1.0e-4, 150.0),
            LaserRingMeasurement::default(),
        )
        .unwrap_err();
        assert_eq!(error, LaserCalibrationError::AnalyticOracleRequiresZeroBeta);
    }

    #[test]
    fn invalid_measurement_settings_fail_closed() {
        let settings = LaserRingMeasurement {
            duration_s: 8.0,
            step_s: 0.0,
            photon_kick_fraction: 0.30,
        };
        assert!(matches!(
            calibrate_laser_relaxation(normalized_laser(0.0, 150.0), settings),
            Err(LaserCalibrationError::InvalidSetting("step_s"))
        ));
    }
}
