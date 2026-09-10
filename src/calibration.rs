//! Executable calibration targets for resonance-search machinery.
//!
//! These are deliberately simple systems with known answers. NoiseLab should
//! recover their optima from simulation before the same search code is trusted
//! on systems whose best operating point is unknown.

use crate::resonance::{
    InteriorPeak, ResonanceInputError, SweepSample, detect_interior_response_peak,
};
use scirust_sim::{System, simulate};
use std::f64::consts::TAU;

/// A harmonically driven, linearly damped oscillator in normalized form:
///
/// `x'' + 2*zeta*omega_n*x' + omega_n^2*x = a*cos(omega*t)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DrivenLinearOscillator {
    omega_n: f64,
    damping_ratio: f64,
    forcing_acceleration: f64,
    forcing_omega: f64,
}

impl DrivenLinearOscillator {
    /// Validate and construct a driven linear oscillator.
    pub fn new(
        omega_n: f64,
        damping_ratio: f64,
        forcing_acceleration: f64,
        forcing_omega: f64,
    ) -> Result<Self, CalibrationError> {
        require_positive("omega_n", omega_n)?;
        require_nonnegative("damping_ratio", damping_ratio)?;
        require_nonnegative("forcing_acceleration", forcing_acceleration)?;
        require_positive("forcing_omega", forcing_omega)?;
        Ok(Self {
            omega_n,
            damping_ratio,
            forcing_acceleration,
            forcing_omega,
        })
    }

    /// Natural angular frequency in radians per unit time.
    pub fn omega_n(self) -> f64 {
        self.omega_n
    }

    /// Damping ratio `zeta`.
    pub fn damping_ratio(self) -> f64 {
        self.damping_ratio
    }

    /// Drive angular frequency in radians per unit time.
    pub fn forcing_omega(self) -> f64 {
        self.forcing_omega
    }

    /// Closed-form steady-state displacement amplitude for this linear model.
    pub fn analytic_steady_state_amplitude(self) -> f64 {
        let detuning = self.omega_n * self.omega_n - self.forcing_omega * self.forcing_omega;
        let damping =
            2.0 * self.damping_ratio * self.omega_n * self.forcing_omega;
        self.forcing_acceleration / (detuning * detuning + damping * damping).sqrt()
    }
}

impl System for DrivenLinearOscillator {
    fn dim(&self) -> usize {
        2
    }

    fn derivatives(&self, t: f64, y: &[f64], dydt: &mut [f64]) {
        dydt[0] = y[1];
        dydt[1] = self.forcing_acceleration * (self.forcing_omega * t).cos()
            - 2.0 * self.damping_ratio * self.omega_n * y[1]
            - self.omega_n * self.omega_n * y[0];
    }
}

/// Numerical measurement settings for one driven-oscillator frequency.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResponseMeasurement {
    /// Number of exponential damping time constants discarded as transient.
    pub settling_time_constants: f64,
    /// Number of drive periods retained for amplitude measurement.
    pub measurement_periods: usize,
    /// Integration samples per shortest relevant period.
    pub samples_per_period: usize,
}

impl Default for ResponseMeasurement {
    fn default() -> Self {
        Self {
            settling_time_constants: 10.0,
            measurement_periods: 8,
            samples_per_period: 160,
        }
    }
}

impl ResponseMeasurement {
    fn validate(self) -> Result<Self, CalibrationError> {
        require_positive("settling_time_constants", self.settling_time_constants)?;
        if self.measurement_periods < 2 {
            return Err(CalibrationError::TooFewMeasurementPeriods(
                self.measurement_periods,
            ));
        }
        if self.samples_per_period < 32 {
            return Err(CalibrationError::TooFewSamplesPerPeriod(
                self.samples_per_period,
            ));
        }
        Ok(self)
    }
}

/// Failures produced by executable resonance calibrations.
#[derive(Debug)]
pub enum CalibrationError {
    /// Invalid named scalar.
    InvalidScalar {
        /// Stable field name.
        name: &'static str,
        /// Rejected value.
        value: f64,
        /// Whether zero would have been accepted.
        zero_allowed: bool,
    },
    /// Fewer than two periods were requested for the response measurement.
    TooFewMeasurementPeriods(usize),
    /// Too few integration samples per period were requested.
    TooFewSamplesPerPeriod(usize),
    /// A zero-damping model has no finite exponential settling time.
    UndampedSettlingUndefined,
    /// SciRust's integration engine rejected or failed the run.
    Simulation(String),
    /// The simulated trajectory did not contain a usable measurement window.
    EmptyMeasurementWindow,
    /// The generic interior-peak detector rejected the sweep.
    PeakDetection(ResonanceInputError),
}

/// Simulate one oscillator and estimate its steady-state displacement amplitude
/// from half the peak-to-peak response after an explicit transient discard.
pub fn measure_steady_state_displacement_amplitude(
    oscillator: DrivenLinearOscillator,
    settings: ResponseMeasurement,
) -> Result<f64, CalibrationError> {
    let settings = settings.validate()?;
    if oscillator.damping_ratio == 0.0 {
        return Err(CalibrationError::UndampedSettlingUndefined);
    }

    let decay_rate = oscillator.damping_ratio * oscillator.omega_n;
    let settle_time = settings.settling_time_constants / decay_rate;
    let drive_period = TAU / oscillator.forcing_omega;
    let measurement_time = settings.measurement_periods as f64 * drive_period;
    let t_end = settle_time + measurement_time;

    let fastest_omega = oscillator.omega_n.max(oscillator.forcing_omega);
    let fastest_period = TAU / fastest_omega;
    let step = fastest_period / settings.samples_per_period as f64;

    let trajectory = simulate(&oscillator, &[0.0, 0.0], 0.0, t_end, step)
        .map_err(|error| CalibrationError::Simulation(error.to_string()))?;

    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    let mut retained = 0usize;
    for (&time, state) in trajectory.t.iter().zip(&trajectory.y) {
        if time >= settle_time {
            let displacement = state[0];
            minimum = minimum.min(displacement);
            maximum = maximum.max(displacement);
            retained += 1;
        }
    }
    if retained < 2 {
        return Err(CalibrationError::EmptyMeasurementWindow);
    }
    Ok(0.5 * (maximum - minimum))
}

/// Result of a blind frequency sweep over the driven linear oscillator.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearResonanceCalibration {
    /// Simulated response samples keyed by angular drive frequency.
    pub sweep: Vec<SweepSample>,
    /// Strongest sampled interior peak, if one was found.
    pub peak: Option<InteriorPeak>,
}

/// Sweep angular drive frequencies and recover the strongest interior response
/// peak entirely from SciRust-integrated trajectories.
pub fn sweep_linear_resonance(
    omega_n: f64,
    damping_ratio: f64,
    forcing_acceleration: f64,
    drive_omegas: &[f64],
    settings: ResponseMeasurement,
    minimum_prominence: f64,
) -> Result<LinearResonanceCalibration, CalibrationError> {
    require_positive("omega_n", omega_n)?;
    require_nonnegative("damping_ratio", damping_ratio)?;
    require_nonnegative("forcing_acceleration", forcing_acceleration)?;
    if drive_omegas.len() < 3 {
        return Err(CalibrationError::PeakDetection(
            ResonanceInputError::SweepTooShort,
        ));
    }

    let mut sweep = Vec::with_capacity(drive_omegas.len());
    for &drive_omega in drive_omegas {
        let oscillator = DrivenLinearOscillator::new(
            omega_n,
            damping_ratio,
            forcing_acceleration,
            drive_omega,
        )?;
        let response = measure_steady_state_displacement_amplitude(oscillator, settings)?;
        sweep.push(SweepSample {
            coordinate: drive_omega,
            response,
        });
    }
    let peak = detect_interior_response_peak(&sweep, minimum_prominence)
        .map_err(CalibrationError::PeakDetection)?;
    Ok(LinearResonanceCalibration { sweep, peak })
}

fn require_positive(name: &'static str, value: f64) -> Result<(), CalibrationError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(CalibrationError::InvalidScalar {
            name,
            value,
            zero_allowed: false,
        });
    }
    Ok(())
}

fn require_nonnegative(name: &'static str, value: f64) -> Result<(), CalibrationError> {
    if !value.is_finite() || value < 0.0 {
        return Err(CalibrationError::InvalidScalar {
            name,
            value,
            zero_allowed: true,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_displacement_resonance_angular_frequency;

    #[test]
    fn numerical_amplitude_matches_linear_closed_form_away_from_transient() {
        let oscillator = DrivenLinearOscillator::new(10.0, 0.2, 1.0, 8.0).unwrap();
        let measured = measure_steady_state_displacement_amplitude(
            oscillator,
            ResponseMeasurement::default(),
        )
        .unwrap();
        let expected = oscillator.analytic_steady_state_amplitude();
        let relative_error = (measured - expected).abs() / expected;
        assert!(
            relative_error < 0.01,
            "measured={measured}, expected={expected}, relative_error={relative_error}"
        );
    }

    #[test]
    fn blind_sweep_recovers_the_known_displacement_resonance() {
        let omega_n = 10.0;
        let zeta = 0.2;
        let expected = linear_displacement_resonance_angular_frequency(omega_n, zeta)
            .unwrap()
            .unwrap();
        let drive_omegas: Vec<f64> = (0..=40).map(|index| 7.5 + 0.1 * index as f64).collect();
        let calibration = sweep_linear_resonance(
            omega_n,
            zeta,
            1.0,
            &drive_omegas,
            ResponseMeasurement::default(),
            1.0e-5,
        )
        .unwrap();
        let measured = calibration.peak.expect("interior peak").coordinate;
        assert!(
            (measured - expected).abs() <= 0.11,
            "measured omega={measured}, analytic omega={expected}"
        );
    }

    #[test]
    fn high_damping_sweep_has_no_false_interior_resonance_peak() {
        let drive_omegas: Vec<f64> = (0..=30).map(|index| 0.5 + 0.25 * index as f64).collect();
        let calibration = sweep_linear_resonance(
            10.0,
            0.8,
            1.0,
            &drive_omegas,
            ResponseMeasurement::default(),
            1.0e-6,
        )
        .unwrap();
        assert!(calibration.peak.is_none());
    }
}
