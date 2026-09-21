//! Stage 0.6 calibration: residual-series ↔ discrete-state mutual information.
//!
//! This module wires NoiseLab residual-family plants (driven oscillator,
//! bistable Langevin, FitzHugh–Nagumo, semiconductor laser) to the Stage 0
//! equal-width histogram mutual-information estimator, with optional Stage 0.5
//! identity-filter audit. It is an exploratory **known-answer** calibration
//! panel: it does **not** claim that natural residuals always encode hidden
//! state, establish causality, or discover system structure.
//!
//! ## Protocol
//!
//! For each preregistered residual family the panel:
//!
//! 1. Generates a residual series `R` under frozen parameters, seeds and
//!    extraction rules (shorter than Stage U2 for cheap calibration).
//! 2. Declares a balanced discrete label `Z` (block-alternating bit).
//! 3. Unit-normalizes `R` to empirical standard deviation one, then builds a
//!    **positive** residual `N = R_unit + A · (2Z − 1)` that injects one
//!    labelled bit by design.
//! 4. Builds a **negative** control that keeps the same `N` but Fisher–Yates
//!    shuffles `Z` under a frozen seed, destroying label dependence.
//! 5. Optionally audits Stage 0.5 identity filtering on the positive pair
//!    (expect `delta_I ≈ 0` under frozen raw bin edges).
//!
//! Classification is fail-closed against frozen MI thresholds. Thresholds,
//! bins, seeds and injection amplitude must not be retuned after inspecting
//! outcomes.

use crate::filtered_information::{
    filtered_histogram_mutual_information_bits, FilteredInformationError,
    FilteredMutualInformationAudit, FrequencySelectiveFilterSpec,
};
use crate::information::{histogram_mutual_information_bits, InformationError};
use crate::u2_readiness::U2_SCIRUST_REVISION;
use crate::{
    simulate_double_well, DoubleWellLangevin, DrivenLinearOscillator, LangevinError, LangevinRun,
};
use scirust_sim::laser::{LaserParams, SemiconductorLaser};
use scirust_sim::{simulate, SplitMix64};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Module identity retained in every panel provenance record.
pub const RESIDUAL_STATE_MODULE: &str = "residual_state_information";

/// Retained residual length for Stage 0.6 (shorter than Stage U2 by design).
pub const RESIDUAL_STATE_SAMPLES: usize = 2_048;
/// Leading samples discarded before residual extraction.
pub const RESIDUAL_STATE_BURN_IN: usize = 256;
/// Frozen equal-width bin count for the Stage 0 histogram estimator.
pub const RESIDUAL_STATE_BINS: usize = 8;
/// Root from which per-family data and shuffle seeds are derived.
pub const RESIDUAL_STATE_SEED_ROOT: u64 = 0x5230_365f_5354_4154;
/// Additive amplitude of the injected labelled bit after unit-variance
/// normalization of the base residual `R`.
pub const RESIDUAL_STATE_INJECTION_AMPLITUDE: f64 = 2.0;
/// Block length for the declared balanced binary label `Z`.
pub const RESIDUAL_STATE_LABEL_BLOCK: usize = 64;
/// Minimum positive-control MI (bits) required to pass known-answer.
pub const RESIDUAL_STATE_POSITIVE_MI_MIN_BITS: f64 = 0.25;
/// Maximum negative-control MI (bits) allowed to pass known-answer.
pub const RESIDUAL_STATE_NEGATIVE_MI_MAX_BITS: f64 = 0.08;
/// Maximum absolute identity-filter delta (bits) allowed when audited.
pub const RESIDUAL_STATE_IDENTITY_ABS_DELTA_MAX_BITS: f64 = 1.0e-12;

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

/// Preregistered residual families for the Stage 0.6 known-answer panel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ResidualStateFamily {
    DrivenDampedOscillator,
    BistableLangevin,
    FitzHughNagumo,
    SemiconductorLaser,
}

impl ResidualStateFamily {
    /// Stable protocol label for provenance.
    pub fn label(self) -> &'static str {
        match self {
            Self::DrivenDampedOscillator => "driven_damped_oscillator",
            Self::BistableLangevin => "bistable_langevin",
            Self::FitzHughNagumo => "fitzhugh_nagumo",
            Self::SemiconductorLaser => "semiconductor_laser",
        }
    }
}

/// Declaration order of the Stage 0.6 panel (frozen).
pub const RESIDUAL_STATE_PANEL_FAMILIES: [ResidualStateFamily; 4] = [
    ResidualStateFamily::DrivenDampedOscillator,
    ResidualStateFamily::BistableLangevin,
    ResidualStateFamily::FitzHughNagumo,
    ResidualStateFamily::SemiconductorLaser,
];

/// Known-answer classification for one family under frozen thresholds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidualStateKnownAnswerClass {
    /// Positive MI above threshold and negative MI below threshold (and
    /// identity delta within bound when audited).
    KnownAnswerPassed,
    /// Estimators ran but controls missed frozen thresholds.
    KnownAnswerFailed,
}

/// Control arm identity retained in provenance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidualStateControlKind {
    /// `N = unit(R) + A·(2Z−1)` with declared balanced `Z`.
    PositiveInjectedBit,
    /// Same `N`, Fisher–Yates-shuffled `Z` under a frozen seed.
    NegativeShuffledLabel,
}

impl ResidualStateControlKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::PositiveInjectedBit => "positive_injected_bit",
            Self::NegativeShuffledLabel => "negative_shuffled_label",
        }
    }
}

/// Fail-closed errors for the Stage 0.6 panel.
#[derive(Debug, Clone, PartialEq)]
pub enum ResidualStateInformationError {
    Information(InformationError),
    Filtered(FilteredInformationError),
    Langevin(LangevinError),
    Calibration(String),
    Laser(String),
    Simulation(String),
    NonFiniteState {
        family: &'static str,
        step: usize,
    },
    DegenerateResidual {
        family: ResidualStateFamily,
    },
    LengthMismatch {
        family: ResidualStateFamily,
        got: usize,
        expected: usize,
    },
    ProtocolBins {
        bins: usize,
    },
}

impl Display for ResidualStateInformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Information(error) => write!(f, "residual-state information error: {error}"),
            Self::Filtered(error) => write!(f, "residual-state filtered audit error: {error}"),
            Self::Langevin(error) => write!(f, "residual-state Langevin error: {error}"),
            Self::Calibration(message) => write!(f, "residual-state oscillator error: {message}"),
            Self::Laser(message) => write!(f, "residual-state laser error: {message}"),
            Self::Simulation(message) => write!(f, "residual-state simulation error: {message}"),
            Self::NonFiniteState { family, step } => {
                write!(f, "{family} residual became non-finite at step {step}")
            }
            Self::DegenerateResidual { family } => {
                write!(f, "{} residual has zero empirical variance", family.label())
            }
            Self::LengthMismatch {
                family,
                got,
                expected,
            } => write!(
                f,
                "{} residual length {got} != expected {expected}",
                family.label()
            ),
            Self::ProtocolBins { bins } => write!(
                f,
                "Stage 0.6 requires at least 2 histogram bins, got {bins}"
            ),
        }
    }
}

impl Error for ResidualStateInformationError {}

impl From<InformationError> for ResidualStateInformationError {
    fn from(value: InformationError) -> Self {
        Self::Information(value)
    }
}

impl From<FilteredInformationError> for ResidualStateInformationError {
    fn from(value: FilteredInformationError) -> Self {
        Self::Filtered(value)
    }
}

impl From<LangevinError> for ResidualStateInformationError {
    fn from(value: LangevinError) -> Self {
        Self::Langevin(value)
    }
}

/// Provenance attached to one family calibration result.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualStateProvenance {
    pub module: &'static str,
    pub family: ResidualStateFamily,
    pub family_label: &'static str,
    pub scirust_revision: &'static str,
    pub data_seed: u64,
    pub shuffle_seed: u64,
    pub samples: usize,
    pub burn_in_samples: usize,
    pub bins: usize,
    pub injection_amplitude: f64,
    pub label_block: usize,
    pub extraction_rule: &'static str,
    pub parameter_summary: &'static str,
    pub positive_mi_min_bits: f64,
    pub negative_mi_max_bits: f64,
    pub identity_abs_delta_max_bits: f64,
}

/// One family's positive/negative known-answer MI calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualStateFamilyResult {
    pub family: ResidualStateFamily,
    pub positive_mutual_information_bits: f64,
    pub negative_mutual_information_bits: f64,
    pub hidden_state_entropy_bits: f64,
    pub identity_audit: Option<FilteredMutualInformationAudit>,
    pub classification: ResidualStateKnownAnswerClass,
    pub provenance: ResidualStateProvenance,
}

/// Full Stage 0.6 panel over the preregistered residual families.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidualStatePanelResult {
    pub families: Vec<ResidualStateFamilyResult>,
    pub all_passed: bool,
    pub identity_audit_enabled: bool,
}

/// Run the Stage 0.6 residual-state known-answer panel.
///
/// When `identity_audit` is true, each family also runs a Stage 0.5 identity
/// filter audit on the positive control (expect near-zero delta). Protocol
/// inputs are the frozen module constants; callers must not retune them after
/// inspecting outcomes.
pub fn run_residual_state_information_panel(
    identity_audit: bool,
) -> Result<ResidualStatePanelResult, ResidualStateInformationError> {
    run_residual_state_information_panel_with_root(RESIDUAL_STATE_SEED_ROOT, identity_audit)
}

/// Same panel with an explicit seed root (tests / deterministic overrides).
pub fn run_residual_state_information_panel_with_root(
    seed_root: u64,
    identity_audit: bool,
) -> Result<ResidualStatePanelResult, ResidualStateInformationError> {
    if RESIDUAL_STATE_BINS < 2 {
        return Err(ResidualStateInformationError::ProtocolBins {
            bins: RESIDUAL_STATE_BINS,
        });
    }

    let mut families = Vec::with_capacity(RESIDUAL_STATE_PANEL_FAMILIES.len());
    for (index, family) in RESIDUAL_STATE_PANEL_FAMILIES.iter().copied().enumerate() {
        let data_seed = derive_seed(seed_root, (index as u64).wrapping_add(1));
        let shuffle_seed = derive_seed(seed_root, (index as u64).wrapping_add(101));
        families.push(calibrate_family(
            family,
            data_seed,
            shuffle_seed,
            identity_audit,
        )?);
    }
    let all_passed = families
        .iter()
        .all(|result| result.classification == ResidualStateKnownAnswerClass::KnownAnswerPassed);
    Ok(ResidualStatePanelResult {
        families,
        all_passed,
        identity_audit_enabled: identity_audit,
    })
}

/// Calibrate a single residual family under the frozen known-answer protocol.
pub fn calibrate_residual_state_family(
    family: ResidualStateFamily,
    identity_audit: bool,
) -> Result<ResidualStateFamilyResult, ResidualStateInformationError> {
    let index = RESIDUAL_STATE_PANEL_FAMILIES
        .iter()
        .position(|&candidate| candidate == family)
        .expect("panel families cover every ResidualStateFamily variant");
    let data_seed = derive_seed(RESIDUAL_STATE_SEED_ROOT, (index as u64).wrapping_add(1));
    let shuffle_seed = derive_seed(RESIDUAL_STATE_SEED_ROOT, (index as u64).wrapping_add(101));
    calibrate_family(family, data_seed, shuffle_seed, identity_audit)
}

fn calibrate_family(
    family: ResidualStateFamily,
    data_seed: u64,
    shuffle_seed: u64,
    identity_audit: bool,
) -> Result<ResidualStateFamilyResult, ResidualStateInformationError> {
    let (base_residual, extraction_rule, parameter_summary) =
        generate_family_residual(family, data_seed)?;
    let unit_residual = unit_variance_residual(family, &base_residual)?;
    let labels = declared_balanced_labels(unit_residual.len());
    let positive_observation =
        inject_labelled_bit(&unit_residual, &labels, RESIDUAL_STATE_INJECTION_AMPLITUDE);
    let mut negative_labels = labels.clone();
    fisher_yates_shuffle(&mut negative_labels, &mut SplitMix64::new(shuffle_seed));

    let positive_mi =
        histogram_mutual_information_bits(&positive_observation, &labels, RESIDUAL_STATE_BINS)?;
    let negative_mi = histogram_mutual_information_bits(
        &positive_observation,
        &negative_labels,
        RESIDUAL_STATE_BINS,
    )?;
    let hidden_state_entropy_bits = crate::discrete_entropy_bits(&labels)?;

    let identity = if identity_audit {
        Some(filtered_histogram_mutual_information_bits(
            &positive_observation,
            &labels,
            RESIDUAL_STATE_BINS,
            &FrequencySelectiveFilterSpec::Identity,
        )?)
    } else {
        None
    };

    let identity_ok = match &identity {
        None => true,
        Some(audit) => {
            audit.information_delta_bits.abs() <= RESIDUAL_STATE_IDENTITY_ABS_DELTA_MAX_BITS
        }
    };
    let classification = if positive_mi >= RESIDUAL_STATE_POSITIVE_MI_MIN_BITS
        && negative_mi <= RESIDUAL_STATE_NEGATIVE_MI_MAX_BITS
        && identity_ok
    {
        ResidualStateKnownAnswerClass::KnownAnswerPassed
    } else {
        ResidualStateKnownAnswerClass::KnownAnswerFailed
    };

    Ok(ResidualStateFamilyResult {
        family,
        positive_mutual_information_bits: positive_mi,
        negative_mutual_information_bits: negative_mi,
        hidden_state_entropy_bits,
        identity_audit: identity,
        classification,
        provenance: ResidualStateProvenance {
            module: RESIDUAL_STATE_MODULE,
            family,
            family_label: family.label(),
            scirust_revision: U2_SCIRUST_REVISION,
            data_seed,
            shuffle_seed,
            samples: RESIDUAL_STATE_SAMPLES,
            burn_in_samples: RESIDUAL_STATE_BURN_IN,
            bins: RESIDUAL_STATE_BINS,
            injection_amplitude: RESIDUAL_STATE_INJECTION_AMPLITUDE,
            label_block: RESIDUAL_STATE_LABEL_BLOCK,
            extraction_rule,
            parameter_summary,
            positive_mi_min_bits: RESIDUAL_STATE_POSITIVE_MI_MIN_BITS,
            negative_mi_max_bits: RESIDUAL_STATE_NEGATIVE_MI_MAX_BITS,
            identity_abs_delta_max_bits: RESIDUAL_STATE_IDENTITY_ABS_DELTA_MAX_BITS,
        },
    })
}

fn generate_family_residual(
    family: ResidualStateFamily,
    seed: u64,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
    match family {
        ResidualStateFamily::DrivenDampedOscillator => generate_oscillator_residual(seed),
        ResidualStateFamily::BistableLangevin => generate_bistable_residual(seed),
        ResidualStateFamily::FitzHughNagumo => generate_fhn_residual(seed),
        ResidualStateFamily::SemiconductorLaser => generate_laser_residual(seed),
    }
}

fn generate_oscillator_residual(
    seed: u64,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
    DrivenLinearOscillator::new(
        OSC_OMEGA_N,
        OSC_DAMPING_RATIO,
        OSC_FORCING_ACCEL,
        OSC_FORCING_OMEGA,
    )
    .map_err(|error| ResidualStateInformationError::Calibration(format!("{error:?}")))?;

    let total = RESIDUAL_STATE_SAMPLES
        .checked_add(RESIDUAL_STATE_BURN_IN)
        .expect("Stage 0.6 oscillator sample budget fits in usize");
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
            return Err(ResidualStateInformationError::NonFiniteState {
                family: "driven_damped_oscillator",
                step,
            });
        }
        trajectory.push(x);
    }

    let residual = window_mean_residual(retain_tail(&trajectory, RESIDUAL_STATE_SAMPLES));
    finish_residual(
        ResidualStateFamily::DrivenDampedOscillator,
        residual,
        "post-burn-in displacement minus retained-window sample mean",
        "omega_n=1, zeta=0.05, a=1, omega=0.95, sigma=0.05, dt=0.05, samples=2048, burn_in=256",
    )
}

fn generate_laser_residual(
    seed: u64,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
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
    let model = SemiconductorLaser::new(params)
        .map_err(|error| ResidualStateInformationError::Laser(error.to_string()))?;
    let steady_photon = model.steady_state_photon_density();
    if !steady_photon.is_finite() || steady_photon <= 0.0 {
        return Err(ResidualStateInformationError::Laser(
            "laser residual requires above-threshold steady photon density".into(),
        ));
    }
    let steady_carrier = model.steady_state_carrier_density();
    let initial_photon = steady_photon * (1.0 - LASER_PHOTON_KICK);
    let initial_state = model.initial_state(steady_carrier, initial_photon);
    let total = RESIDUAL_STATE_SAMPLES
        .checked_add(RESIDUAL_STATE_BURN_IN)
        .expect("Stage 0.6 laser sample budget fits in usize");
    let t_end = (total as f64) * LASER_STEP_S;
    let trajectory = simulate(&model, &initial_state, 0.0, t_end, LASER_STEP_S)
        .map_err(|error| ResidualStateInformationError::Simulation(error.to_string()))?;
    let photon = trajectory.column(1).ok_or_else(|| {
        ResidualStateInformationError::Laser("missing photon-density column".into())
    })?;
    if photon.len() < total {
        return Err(ResidualStateInformationError::Laser(format!(
            "laser trajectory too short: got {} samples, need {total}",
            photon.len()
        )));
    }
    let retained = retain_tail(&photon[..total], RESIDUAL_STATE_SAMPLES);
    let residual: Vec<f64> = retained
        .into_iter()
        .map(|value| value - steady_photon)
        .collect();
    finish_residual(
        ResidualStateFamily::SemiconductorLaser,
        residual,
        "post-burn-in photon density minus analytic steady-state photon density",
        "LaserParams{g0=1,n_t=1,tau_n=1,tau_p=0.01,gamma=1,beta=0,pump=150}, kick=0.30, step=5e-4, samples=2048, burn_in=256",
    )
}

fn generate_bistable_residual(
    seed: u64,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
    let model = DoubleWellLangevin::new(BISTABLE_A, BISTABLE_B, 0.0, 0.1, 0.0)?;
    let total = RESIDUAL_STATE_SAMPLES
        .checked_add(RESIDUAL_STATE_BURN_IN)
        .expect("Stage 0.6 bistable sample budget fits in usize");
    let run = LangevinRun::new(BISTABLE_DT, total, RESIDUAL_STATE_BURN_IN)?;
    let trajectory =
        simulate_double_well(model, BISTABLE_NOISE_D, model.well_location(), run, seed)?;
    let residual = window_mean_residual(retain_tail(&trajectory.x, RESIDUAL_STATE_SAMPLES));
    finish_residual(
        ResidualStateFamily::BistableLangevin,
        residual,
        "post-burn-in double-well state minus retained-window sample mean",
        "a=1,b=1,A=0,f=0.1,D=0.10,dt=0.01,x0=well,samples=2048,burn_in=256",
    )
}

fn generate_fhn_residual(
    seed: u64,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
    let total = RESIDUAL_STATE_SAMPLES
        .checked_add(RESIDUAL_STATE_BURN_IN)
        .expect("Stage 0.6 FHN sample budget fits in usize");
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
            return Err(ResidualStateInformationError::NonFiniteState {
                family: "fitzhugh_nagumo",
                step,
            });
        }
        trajectory.push(x);
    }

    let retained = retain_tail(&trajectory, RESIDUAL_STATE_SAMPLES);
    let residual: Vec<f64> = retained.into_iter().map(|value| value - x_eq).collect();
    finish_residual(
        ResidualStateFamily::FitzHughNagumo,
        residual,
        "post-burn-in activator minus deterministic equilibrium activator x*=-a",
        "epsilon=0.01,a=1.05,sigma=0.05,dt=0.001,samples=2048,burn_in=256",
    )
}

fn finish_residual(
    family: ResidualStateFamily,
    residual: Vec<f64>,
    extraction_rule: &'static str,
    parameter_summary: &'static str,
) -> Result<(Vec<f64>, &'static str, &'static str), ResidualStateInformationError> {
    if residual.len() != RESIDUAL_STATE_SAMPLES {
        return Err(ResidualStateInformationError::LengthMismatch {
            family,
            got: residual.len(),
            expected: RESIDUAL_STATE_SAMPLES,
        });
    }
    if !residual.iter().all(|value| value.is_finite()) {
        return Err(ResidualStateInformationError::NonFiniteState {
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
        return Err(ResidualStateInformationError::DegenerateResidual { family });
    }
    Ok((residual, extraction_rule, parameter_summary))
}

fn declared_balanced_labels(len: usize) -> Vec<usize> {
    (0..len)
        .map(|index| (index / RESIDUAL_STATE_LABEL_BLOCK) % 2)
        .collect()
}

fn unit_variance_residual(
    family: ResidualStateFamily,
    residual: &[f64],
) -> Result<Vec<f64>, ResidualStateInformationError> {
    let n = residual.len() as f64;
    let mean = residual.iter().sum::<f64>() / n;
    let variance = residual
        .iter()
        .map(|value| {
            let centered = value - mean;
            centered * centered
        })
        .sum::<f64>()
        / n;
    if !variance.is_finite() || variance <= 1.0e-15 {
        return Err(ResidualStateInformationError::DegenerateResidual { family });
    }
    let std = variance.sqrt();
    Ok(residual.iter().map(|value| (value - mean) / std).collect())
}

fn inject_labelled_bit(residual: &[f64], labels: &[usize], amplitude: f64) -> Vec<f64> {
    residual
        .iter()
        .zip(labels)
        .map(|(&value, &label)| {
            let signed = if label == 0 { -1.0 } else { 1.0 };
            value + amplitude * signed
        })
        .collect()
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

fn fisher_yates_shuffle(values: &mut [usize], rng: &mut SplitMix64) {
    for index in (1..values.len()).rev() {
        let swap_with = unbiased_index(rng, index + 1);
        values.swap(index, swap_with);
    }
}

fn unbiased_index(rng: &mut SplitMix64, upper_exclusive: usize) -> usize {
    debug_assert!(upper_exclusive > 0);
    let bound = upper_exclusive as u64;
    let threshold = bound.wrapping_neg() % bound;
    loop {
        let value = rng.next_u64();
        if value >= threshold {
            return (value % bound) as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::information::InformationError;

    #[test]
    fn panel_with_identity_audit_passes_all_families() {
        let panel = run_residual_state_information_panel(true).unwrap();
        assert!(panel.identity_audit_enabled);
        assert_eq!(panel.families.len(), 4);
        assert!(panel.all_passed);
        for result in &panel.families {
            assert_eq!(
                result.classification,
                ResidualStateKnownAnswerClass::KnownAnswerPassed
            );
            assert!(
                result.positive_mutual_information_bits >= RESIDUAL_STATE_POSITIVE_MI_MIN_BITS,
                "{} positive MI={}",
                result.family.label(),
                result.positive_mutual_information_bits
            );
            assert!(
                result.negative_mutual_information_bits <= RESIDUAL_STATE_NEGATIVE_MI_MAX_BITS,
                "{} negative MI={}",
                result.family.label(),
                result.negative_mutual_information_bits
            );
            let audit = result.identity_audit.as_ref().unwrap();
            assert!(
                audit.information_delta_bits.abs() <= RESIDUAL_STATE_IDENTITY_ABS_DELTA_MAX_BITS
            );
            assert_eq!(result.provenance.module, RESIDUAL_STATE_MODULE);
            assert_eq!(result.provenance.bins, RESIDUAL_STATE_BINS);
            assert_eq!(result.provenance.scirust_revision, U2_SCIRUST_REVISION);
            assert!((result.hidden_state_entropy_bits - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn panel_is_seed_reproducible() {
        let a = run_residual_state_information_panel(true).unwrap();
        let b = run_residual_state_information_panel(true).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn distinct_roots_change_stochastic_family_scores() {
        let a = run_residual_state_information_panel_with_root(1, false).unwrap();
        let b = run_residual_state_information_panel_with_root(2, false).unwrap();
        // Oscillator (index 0) is stochastic; scores or seeds must move.
        assert_ne!(
            a.families[0].provenance.data_seed,
            b.families[0].provenance.data_seed
        );
        assert_ne!(
            (
                a.families[0].positive_mutual_information_bits.to_bits(),
                a.families[0].negative_mutual_information_bits.to_bits()
            ),
            (
                b.families[0].positive_mutual_information_bits.to_bits(),
                b.families[0].negative_mutual_information_bits.to_bits()
            )
        );
    }

    #[test]
    fn positive_control_exceeds_negative_for_each_family() {
        for family in RESIDUAL_STATE_PANEL_FAMILIES {
            let result = calibrate_residual_state_family(family, false).unwrap();
            assert!(
                result.positive_mutual_information_bits
                    > result.negative_mutual_information_bits + 0.1,
                "{}: pos={}, neg={}",
                family.label(),
                result.positive_mutual_information_bits,
                result.negative_mutual_information_bits
            );
        }
    }

    #[test]
    fn injection_and_shuffle_helpers_are_length_preserving() {
        let residual = vec![0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8];
        let labels = declared_balanced_labels(residual.len());
        assert_eq!(labels.len(), residual.len());
        let injected = inject_labelled_bit(&residual, &labels, 2.0);
        assert_eq!(injected.len(), residual.len());
        let mut shuffled = labels.clone();
        fisher_yates_shuffle(&mut shuffled, &mut SplitMix64::new(99));
        assert_eq!(shuffled.len(), labels.len());
        let mut sorted_a = labels.clone();
        let mut sorted_b = shuffled.clone();
        sorted_a.sort_unstable();
        sorted_b.sort_unstable();
        assert_eq!(sorted_a, sorted_b);
    }

    #[test]
    fn malformed_estimator_inputs_fail_closed() {
        assert_eq!(
            histogram_mutual_information_bits(&[], &[], RESIDUAL_STATE_BINS),
            Err(InformationError::EmptyInput)
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0, 1.0], &[0], RESIDUAL_STATE_BINS),
            Err(InformationError::LengthMismatch {
                observed: 2,
                hidden_state: 1,
            })
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0, 1.0], &[0, 1], 1),
            Err(InformationError::TooFewBins { bins: 1 })
        );
    }

    #[test]
    fn length_mismatch_on_finish_is_protocol_error() {
        let err = finish_residual(
            ResidualStateFamily::DrivenDampedOscillator,
            vec![1.0, 2.0],
            "test",
            "test",
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ResidualStateInformationError::LengthMismatch {
                family: ResidualStateFamily::DrivenDampedOscillator,
                got: 2,
                expected: RESIDUAL_STATE_SAMPLES,
            }
        ));
    }

    #[test]
    fn degenerate_residual_fails_closed() {
        let err = finish_residual(
            ResidualStateFamily::BistableLangevin,
            vec![0.0; RESIDUAL_STATE_SAMPLES],
            "test",
            "test",
        )
        .unwrap_err();
        assert!(matches!(
            err,
            ResidualStateInformationError::DegenerateResidual {
                family: ResidualStateFamily::BistableLangevin,
            }
        ));
    }

    #[test]
    fn control_kind_labels_are_stable() {
        assert_eq!(
            ResidualStateControlKind::PositiveInjectedBit.label(),
            "positive_injected_bit"
        );
        assert_eq!(
            ResidualStateControlKind::NegativeShuffledLabel.label(),
            "negative_shuffled_label"
        );
    }
}
