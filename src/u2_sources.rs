//! Frozen residual-series generators for the Stage U2 cross-mechanism panel.
//!
//! Parameters, seeds and extraction rules are fixed **outcome-blind**: they are
//! not chosen or revised after inspecting U2 scores or decisions. Each family
//! records provenance fields so a later science run can cite exact settings.
//!
//! ## Frozen residual definitions
//!
//! 1. **Driven damped oscillator** — Euler–Maruyama integration of the driven
//!    linear oscillator with additive velocity noise (`SplitMix64`). Residual =
//!    post-burn-in displacement minus the retained-window sample mean.
//! 2. **Semiconductor laser** — SciRust rate-equation photon-density trajectory
//!    after a fractional photon kick from the above-threshold steady state.
//!    Residual = post-burn-in `s(t) - s_ss` (analytic steady photon density).
//! 3. **Bistable Langevin** — unforced double-well trajectory via
//!    [`simulate_double_well`]. Residual = post-burn-in state minus the
//!    retained-window sample mean.
//! 4. **FitzHugh–Nagumo** — same Euler / Euler–Maruyama discretization as
//!    [`simulate_fhn_spikes`], recording the activator. Residual = post-burn-in
//!    `x(t) - x*` where `x*` is the deterministic equilibrium activator.

use crate::u2_plan::U2SourceFamily;
use crate::u2_readiness::U2_SCIRUST_REVISION;
use crate::{
    simulate_double_well, DoubleWellLangevin, DrivenLinearOscillator, LangevinError, LangevinRun,
};
use scirust_sim::laser::{LaserParams, SemiconductorLaser};
use scirust_sim::{simulate, SplitMix64};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Science-facing retained length (matches Stage U1 panel samples).
pub const U2_SOURCE_SAMPLES: usize = 8_192;
/// Leading samples discarded before residual extraction.
pub const U2_SOURCE_BURN_IN: usize = 1_024;
/// Root from which deterministic per-family data seeds are derived.
pub const U2_DATA_SEED_ROOT: u64 = 0x5532_4e4f_4953_454c;

const OSC_OMEGA_N: f64 = 1.0;
const OSC_DAMPING_RATIO: f64 = 0.05;
const OSC_FORCING_ACCEL: f64 = 1.0;
const OSC_FORCING_OMEGA: f64 = 0.95;
const OSC_NOISE_SIGMA: f64 = 0.05;
const OSC_DT: f64 = 0.05;

const LASER_G0: f64 = 1.0;
const LASER_N_T: f64 = 1.0;
const LASER_TAU_N: f64 = 1.0;
const LASER_TAU_P: f64 = 0.01;
const LASER_GAMMA: f64 = 1.0;
const LASER_BETA: f64 = 0.0;
const LASER_PUMP: f64 = 150.0;
const LASER_STEP_S: f64 = 5.0e-4;
const LASER_PHOTON_KICK: f64 = 0.30;

const BISTABLE_A: f64 = 1.0;
const BISTABLE_B: f64 = 1.0;
const BISTABLE_NOISE_D: f64 = 0.10;
const BISTABLE_DT: f64 = 0.01;

const FHN_EPSILON: f64 = 0.01;
const FHN_A: f64 = 1.05;
const FHN_NOISE_SIGMA: f64 = 0.05;
const FHN_DT: f64 = 0.001;

/// Provenance attached to one frozen residual series.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct U2ResidualProvenance {
    pub family: U2SourceFamily,
    pub scirust_revision: &'static str,
    pub data_seed: u64,
    pub samples: usize,
    pub burn_in_samples: usize,
    pub extraction_rule: &'static str,
    pub parameter_summary: &'static str,
}

/// One outcome-blind residual series with provenance.
#[derive(Debug, Clone, PartialEq)]
pub struct U2ResidualSeries {
    pub family: U2SourceFamily,
    pub residual: Vec<f64>,
    pub provenance: U2ResidualProvenance,
}

#[derive(Debug, Clone, PartialEq)]
pub enum U2SourceError {
    Langevin(LangevinError),
    Calibration(String),
    Laser(String),
    Simulation(String),
    NonFiniteState {
        family: &'static str,
        step: usize,
    },
    DegenerateResidual {
        family: U2SourceFamily,
    },
    LengthMismatch {
        family: U2SourceFamily,
        got: usize,
        expected: usize,
    },
}

impl Display for U2SourceError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Langevin(error) => write!(f, "bistable Langevin residual error: {error}"),
            Self::Calibration(message) => write!(f, "oscillator residual error: {message}"),
            Self::Laser(message) => write!(f, "laser residual error: {message}"),
            Self::Simulation(message) => write!(f, "simulation error: {message}"),
            Self::NonFiniteState { family, step } => {
                write!(f, "{family} residual became non-finite at step {step}")
            }
            Self::DegenerateResidual { family } => {
                write!(f, "{family:?} residual has zero empirical variance")
            }
            Self::LengthMismatch {
                family,
                got,
                expected,
            } => write!(
                f,
                "{family:?} residual length {got} != expected retained length {expected}"
            ),
        }
    }
}

impl Error for U2SourceError {}

impl From<LangevinError> for U2SourceError {
    fn from(value: LangevinError) -> Self {
        Self::Langevin(value)
    }
}

/// Generate all four frozen residual series from a shared data-seed root.
pub fn generate_u2_residuals(data_seed_root: u64) -> Result<[U2ResidualSeries; 4], U2SourceError> {
    Ok([
        generate_driven_oscillator_residual(derive_seed(data_seed_root, 1))?,
        generate_laser_residual(derive_seed(data_seed_root, 2))?,
        generate_bistable_residual(derive_seed(data_seed_root, 3))?,
        generate_fhn_residual(derive_seed(data_seed_root, 4))?,
    ])
}

/// Look up one residual by frozen family identity.
pub fn residual_for_family(
    residuals: &[U2ResidualSeries; 4],
    family: U2SourceFamily,
) -> &U2ResidualSeries {
    residuals
        .iter()
        .find(|series| series.family == family)
        .expect("all four U2 families are generated together")
}

fn generate_driven_oscillator_residual(seed: u64) -> Result<U2ResidualSeries, U2SourceError> {
    // Validate the linear plant through the existing constructor without using
    // its deterministic amplitude helper for the residual itself.
    DrivenLinearOscillator::new(
        OSC_OMEGA_N,
        OSC_DAMPING_RATIO,
        OSC_FORCING_ACCEL,
        OSC_FORCING_OMEGA,
    )
    .map_err(|error| U2SourceError::Calibration(format!("{error:?}")))?;

    let total = U2_SOURCE_SAMPLES
        .checked_add(U2_SOURCE_BURN_IN)
        .expect("U2 oscillator sample budget fits in usize");
    let mut rng = SplitMix64::new(seed);
    let mut x = 0.0;
    let mut v = 0.0;
    let mut trajectory = Vec::with_capacity(total);
    let damping = 2.0 * OSC_DAMPING_RATIO * OSC_OMEGA_N;
    let restoring = OSC_OMEGA_N * OSC_OMEGA_N;
    let noise_scale = OSC_NOISE_SIGMA * OSC_DT.sqrt();

    for step in 0..total {
        let t = step as f64 * OSC_DT;
        let force = OSC_FORCING_ACCEL * (OSC_FORCING_OMEGA * t).cos();
        let dv = (force - damping * v - restoring * x) * OSC_DT + noise_scale * rng.next_gaussian();
        let dx = v * OSC_DT;
        x += dx;
        v += dv;
        if !x.is_finite() || !v.is_finite() {
            return Err(U2SourceError::NonFiniteState {
                family: "driven_damped_oscillator",
                step,
            });
        }
        trajectory.push(x);
    }

    let residual = window_mean_residual(retain_tail(&trajectory, U2_SOURCE_SAMPLES));
    finish_series(
        U2SourceFamily::DrivenDampedOscillator,
        seed,
        residual,
        "post-burn-in displacement minus retained-window sample mean",
        "omega_n=1, zeta=0.05, a=1, omega=0.95, sigma=0.05, dt=0.05, samples=8192, burn_in=1024",
    )
}

fn generate_laser_residual(seed: u64) -> Result<U2ResidualSeries, U2SourceError> {
    // Seed is recorded for provenance; the SciRust laser plant used here is
    // deterministic given parameters and initial kick.
    let _ = seed;
    let params = LaserParams {
        g0: LASER_G0,
        n_t: LASER_N_T,
        tau_n: LASER_TAU_N,
        tau_p: LASER_TAU_P,
        gamma: LASER_GAMMA,
        beta: LASER_BETA,
        pump: LASER_PUMP,
    };
    let model =
        SemiconductorLaser::new(params).map_err(|error| U2SourceError::Laser(error.to_string()))?;
    let steady_photon = model.steady_state_photon_density();
    if !steady_photon.is_finite() || steady_photon <= 0.0 {
        return Err(U2SourceError::Laser(
            "laser residual requires above-threshold steady photon density".into(),
        ));
    }
    let steady_carrier = model.steady_state_carrier_density();
    let initial_photon = steady_photon * (1.0 - LASER_PHOTON_KICK);
    let initial_state = model.initial_state(steady_carrier, initial_photon);
    let total = U2_SOURCE_SAMPLES
        .checked_add(U2_SOURCE_BURN_IN)
        .expect("U2 laser sample budget fits in usize");
    let t_end = (total as f64) * LASER_STEP_S;
    let trajectory = simulate(&model, &initial_state, 0.0, t_end, LASER_STEP_S)
        .map_err(|error| U2SourceError::Simulation(error.to_string()))?;
    let photon = trajectory
        .column(1)
        .ok_or_else(|| U2SourceError::Laser("missing photon-density column".into()))?;
    if photon.len() < total {
        return Err(U2SourceError::Laser(format!(
            "laser trajectory too short: got {} samples, need {total}",
            photon.len()
        )));
    }
    let retained = retain_tail(&photon[..total], U2_SOURCE_SAMPLES);
    let residual: Vec<f64> = retained
        .into_iter()
        .map(|value| value - steady_photon)
        .collect();
    finish_series(
        U2SourceFamily::SemiconductorLaser,
        seed,
        residual,
        "post-burn-in photon density minus analytic steady-state photon density",
        "LaserParams{g0=1,n_t=1,tau_n=1,tau_p=0.01,gamma=1,beta=0,pump=150}, kick=0.30, step=5e-4, samples=8192, burn_in=1024",
    )
}

fn generate_bistable_residual(seed: u64) -> Result<U2ResidualSeries, U2SourceError> {
    let model = DoubleWellLangevin::new(BISTABLE_A, BISTABLE_B, 0.0, 0.1, 0.0)?;
    let total = U2_SOURCE_SAMPLES
        .checked_add(U2_SOURCE_BURN_IN)
        .expect("U2 bistable sample budget fits in usize");
    let run = LangevinRun::new(BISTABLE_DT, total, U2_SOURCE_BURN_IN)?;
    let trajectory =
        simulate_double_well(model, BISTABLE_NOISE_D, model.well_location(), run, seed)?;
    let residual = window_mean_residual(retain_tail(&trajectory.x, U2_SOURCE_SAMPLES));
    finish_series(
        U2SourceFamily::BistableLangevin,
        seed,
        residual,
        "post-burn-in double-well state minus retained-window sample mean",
        "a=1,b=1,A=0,f=0.1,D=0.10,dt=0.01,x0=well,samples=8192,burn_in=1024",
    )
}

fn generate_fhn_residual(seed: u64) -> Result<U2ResidualSeries, U2SourceError> {
    let total = U2_SOURCE_SAMPLES
        .checked_add(U2_SOURCE_BURN_IN)
        .expect("U2 FHN sample budget fits in usize");
    let mut rng = SplitMix64::new(seed);
    let x_eq = -FHN_A;
    let y_eq = FHN_A.powi(3) / 3.0 - FHN_A;
    let mut x = x_eq;
    let mut y = y_eq;
    let stochastic_scale = FHN_NOISE_SIGMA * FHN_DT.sqrt();
    let mut trajectory = Vec::with_capacity(total);

    for step in 0..total {
        let old_x = x;
        let old_y = y;
        let dx = (old_x - old_x * old_x * old_x / 3.0 - old_y) / FHN_EPSILON;
        let dy = old_x + FHN_A;
        x = old_x + FHN_DT * dx;
        y = old_y + FHN_DT * dy + stochastic_scale * rng.next_gaussian();
        if !x.is_finite() || !y.is_finite() {
            return Err(U2SourceError::NonFiniteState {
                family: "fitzhugh_nagumo",
                step,
            });
        }
        trajectory.push(x);
    }

    let retained = retain_tail(&trajectory, U2_SOURCE_SAMPLES);
    let residual: Vec<f64> = retained.into_iter().map(|value| value - x_eq).collect();
    finish_series(
        U2SourceFamily::FitzHughNagumo,
        seed,
        residual,
        "post-burn-in activator minus deterministic equilibrium activator x*=-a",
        "epsilon=0.01,a=1.05,sigma=0.05,dt=0.001,samples=8192,burn_in=1024",
    )
}

fn finish_series(
    family: U2SourceFamily,
    data_seed: u64,
    residual: Vec<f64>,
    extraction_rule: &'static str,
    parameter_summary: &'static str,
) -> Result<U2ResidualSeries, U2SourceError> {
    if residual.len() != U2_SOURCE_SAMPLES {
        return Err(U2SourceError::LengthMismatch {
            family,
            got: residual.len(),
            expected: U2_SOURCE_SAMPLES,
        });
    }
    if !residual.iter().all(|value| value.is_finite()) {
        return Err(U2SourceError::NonFiniteState {
            family: "residual",
            step: 0,
        });
    }
    let mean = residual.iter().sum::<f64>() / residual.len() as f64;
    let variance = residual
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / residual.len() as f64;
    if variance <= 1.0e-15 {
        return Err(U2SourceError::DegenerateResidual { family });
    }
    Ok(U2ResidualSeries {
        family,
        residual,
        provenance: U2ResidualProvenance {
            family,
            scirust_revision: U2_SCIRUST_REVISION,
            data_seed,
            samples: U2_SOURCE_SAMPLES,
            burn_in_samples: U2_SOURCE_BURN_IN,
            extraction_rule,
            parameter_summary,
        },
    })
}

fn window_mean_residual(values: Vec<f64>) -> Vec<f64> {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    values.into_iter().map(|value| value - mean).collect()
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

    #[test]
    fn residuals_are_seed_reproducible_and_correct_length() {
        let first = generate_u2_residuals(U2_DATA_SEED_ROOT).unwrap();
        let second = generate_u2_residuals(U2_DATA_SEED_ROOT).unwrap();
        assert_eq!(first, second);
        for series in &first {
            assert_eq!(series.residual.len(), U2_SOURCE_SAMPLES);
            assert_eq!(series.provenance.scirust_revision, U2_SCIRUST_REVISION);
            assert!(series.residual.iter().all(|value| value.is_finite()));
        }
    }

    #[test]
    fn distinct_roots_change_at_least_one_stochastic_family() {
        let a = generate_u2_residuals(1).unwrap();
        let b = generate_u2_residuals(2).unwrap();
        // Laser is deterministic given params; stochastic families must move.
        assert_ne!(
            a[0].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            b[0].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        assert_ne!(
            a[2].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            b[2].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        assert_ne!(
            a[3].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            b[3].residual
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
    }
}
