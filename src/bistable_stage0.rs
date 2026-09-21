//! Executable Stage 0 v2 panel for bistable Langevin / Kramers calibration.
//!
//! Materializes the frozen factors from
//! `docs/research/bistable-kramers-stage0-v2.md`, runs the preregistered
//! controls fail-closed, and classifies decisions in protocol language.
//! Smoke mode abbreviates periods and seed count and must never be recorded as
//! scientific evidence.
//!
//! This module does not retune grids, claim Kramers universality, or write a
//! results.md. Negative and inconclusive outcomes are retained.

use crate::evidence::{detect_edge_separated_peak, EdgeSeparatedPeak, EvidenceError};
use crate::langevin::{
    coherent_switching_response, simulate_double_well, DoubleWellLangevin, LangevinError,
    LangevinRun, NoiseResponse,
};
use crate::preregistered::{
    BistableStage0V2, BISTABLE_STAGE0_V2_FALSIFICATION_FREQUENCY_HZ,
    BISTABLE_STAGE0_V2_GRID_FACTORS, BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA,
    BISTABLE_STAGE0_V2_SEEDS,
};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Preregistered total forcing periods (scientific).
pub const BISTABLE_STAGE0_TOTAL_PERIODS: usize = 80;
/// Preregistered burn-in forcing periods (scientific).
pub const BISTABLE_STAGE0_BURN_IN_PERIODS: usize = 20;
/// Edge-separation uncertainty weight (`two` seed-level standard errors).
pub const BISTABLE_STAGE0_UNCERTAINTY_WEIGHT: f64 = 2.0;
/// Factor-of-two Kramers compatibility tolerance (may not be widened post hoc).
pub const BISTABLE_STAGE0_KRAMERS_FACTOR_TOLERANCE: f64 = 2.0;

/// Abbreviated total periods for non-scientific smoke only.
pub const BISTABLE_STAGE0_SMOKE_TOTAL_PERIODS: usize = 4;
/// Abbreviated burn-in periods for non-scientific smoke only.
pub const BISTABLE_STAGE0_SMOKE_BURN_IN_PERIODS: usize = 1;
/// Seed count used only by the non-scientific smoke path.
pub const BISTABLE_STAGE0_SMOKE_SEED_COUNT: usize = 2;

/// Whether this panel invocation may claim scientific Stage 0 evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BistableStage0Mode {
    /// Exact preregistered periods, seeds, grid, and decision contract.
    Scientific,
    /// Explicitly non-scientific CI / developer smoke (abbreviated budget).
    NonScientificSmoke,
}

impl BistableStage0Mode {
    #[must_use]
    pub const fn is_scientific(self) -> bool {
        matches!(self, Self::Scientific)
    }

    #[must_use]
    pub const fn total_periods(self) -> usize {
        match self {
            Self::Scientific => BISTABLE_STAGE0_TOTAL_PERIODS,
            Self::NonScientificSmoke => BISTABLE_STAGE0_SMOKE_TOTAL_PERIODS,
        }
    }

    #[must_use]
    pub const fn burn_in_periods(self) -> usize {
        match self {
            Self::Scientific => BISTABLE_STAGE0_BURN_IN_PERIODS,
            Self::NonScientificSmoke => BISTABLE_STAGE0_SMOKE_BURN_IN_PERIODS,
        }
    }

    #[must_use]
    pub const fn seed_count(self) -> usize {
        match self {
            Self::Scientific => BISTABLE_STAGE0_V2_SEEDS.len(),
            Self::NonScientificSmoke => BISTABLE_STAGE0_SMOKE_SEED_COUNT,
        }
    }
}

/// Configuration frozen before inspecting Stage 0 v2 outcomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BistableStage0Config {
    pub mode: BistableStage0Mode,
}

impl Default for BistableStage0Config {
    fn default() -> Self {
        Self::scientific()
    }
}

impl BistableStage0Config {
    #[must_use]
    pub const fn scientific() -> Self {
        Self {
            mode: BistableStage0Mode::Scientific,
        }
    }

    #[must_use]
    pub const fn non_scientific_smoke() -> Self {
        Self {
            mode: BistableStage0Mode::NonScientificSmoke,
        }
    }
}

/// Analytic controls recorded before any trajectory is generated.
#[derive(Debug, Clone, PartialEq)]
pub struct BistableStage0AnalyticControls {
    pub barrier_height: f64,
    pub kramers_prefactor_rate: f64,
    pub predicted_noise_intensity: f64,
    pub positive_noise_grid: [f64; 7],
    pub grid_factors: [f64; 7],
    pub protocol_blob_sha: &'static str,
    pub forcing_is_subthreshold: bool,
}

/// One noise-grid point with seed-level amplitudes retained for the artifact.
#[derive(Debug, Clone, PartialEq)]
pub struct BistableStage0GridPoint {
    pub noise_intensity: f64,
    pub mean_coherent_amplitude: f64,
    pub sample_stddev: f64,
    pub standard_error: f64,
    pub replicates: usize,
    pub per_seed_amplitudes: Vec<f64>,
}

impl BistableStage0GridPoint {
    fn to_noise_response(&self) -> NoiseResponse {
        NoiseResponse {
            noise_intensity: self.noise_intensity,
            mean_coherent_amplitude: self.mean_coherent_amplitude,
            sample_stddev: self.sample_stddev,
            replicates: self.replicates,
        }
    }
}

/// Outcome of the deterministic `D = 0` control.
#[derive(Debug, Clone, PartialEq)]
pub enum BistableStage0DeterministicControl {
    /// All declared seeds produced identical trajectories.
    SeedIndependent { amplitude: f64 },
    /// Trajectories differed across seeds at zero noise (protocol failure).
    SeedDependent,
    /// Simulation could not complete.
    Failed { detail: String },
}

/// Falsification-regime preflight / execution status.
#[derive(Debug, Clone, PartialEq)]
pub enum BistableStage0FalsificationStatus {
    /// Analytic half-period match has no positive finite solution; no grid run.
    KramersControlInapplicable,
    /// Analytic control existed and the abbreviated/scientific panel was run.
    Executed {
        predicted_noise_intensity: f64,
        decision: BistableStage0Decision,
    },
    Failed {
        detail: String,
    },
}

/// Protocol-language Stage 0 v2 decision for one forced regime.
#[derive(Debug, Clone, PartialEq)]
pub enum BistableStage0Decision {
    /// Finite-grid maximizer only; does not imply resonance or H0 rejection.
    BestSampledPoint {
        noise_intensity: f64,
        mean_coherent_amplitude: f64,
        grid_index: usize,
    },
    /// No uncertainty-separated interior peak; H0 is not rejected.
    H0NotRejected {
        best_sampled: Option<BestSampledSummary>,
    },
    /// Interior peak found, but outside the factor-of-2 Kramers window.
    InteriorPeakWithoutMechanismCompatibility {
        peak: EdgeSeparatedPeak,
        predicted_noise_intensity: f64,
        best_sampled: BestSampledSummary,
    },
    /// Interior peak within factor-of-2 of the independent Kramers control.
    MechanismCompatibleInteriorPeak {
        peak: EdgeSeparatedPeak,
        predicted_noise_intensity: f64,
        best_sampled: BestSampledSummary,
    },
    /// Controls or protocol construction failed closed.
    ControlFailure {
        detail: String,
    },
    ProtocolMismatch {
        detail: String,
    },
}

/// Compact maximizer summary retained alongside peak decisions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BestSampledSummary {
    pub noise_intensity: f64,
    pub mean_coherent_amplitude: f64,
    pub grid_index: usize,
}

/// Panel-level summary after primary, controls, and falsification preflight.
#[derive(Debug, Clone, PartialEq)]
pub enum BistableStage0Summary {
    /// Scientific path: primary regime classified under the frozen contract.
    Classified(BistableStage0Decision),
    /// Smoke path: scientific decision gate is not applied to abbreviated budget.
    NonScientificSmoke,
    /// Panel could not complete required controls.
    Incomplete { detail: String },
}

/// Complete bistable Stage 0 v2 panel result.
#[derive(Debug, Clone, PartialEq)]
pub struct BistableStage0Result {
    pub mode: BistableStage0Mode,
    pub scientific_claim_permitted: bool,
    pub analytic: BistableStage0AnalyticControls,
    pub total_periods: usize,
    pub burn_in_periods: usize,
    pub seed_count: usize,
    pub seeds: Vec<u64>,
    pub deterministic_control: BistableStage0DeterministicControl,
    pub primary_grid: Vec<BistableStage0GridPoint>,
    pub zero_forcing_grid: Vec<BistableStage0GridPoint>,
    pub primary_decision: BistableStage0Decision,
    pub falsification: BistableStage0FalsificationStatus,
    pub summary: BistableStage0Summary,
}

/// Failures that prevent starting the preregistered panel contract.
#[derive(Debug, Clone, PartialEq)]
pub enum BistableStage0Error {
    Langevin(LangevinError),
    Evidence(EvidenceError),
    Protocol { detail: String },
}

impl Display for BistableStage0Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Langevin(error) => write!(f, "bistable Stage 0 Langevin error: {error}"),
            Self::Evidence(error) => write!(f, "bistable Stage 0 evidence error: {error}"),
            Self::Protocol { detail } => write!(f, "bistable Stage 0 protocol error: {detail}"),
        }
    }
}

impl Error for BistableStage0Error {}

impl From<LangevinError> for BistableStage0Error {
    fn from(value: LangevinError) -> Self {
        Self::Langevin(value)
    }
}

impl From<EvidenceError> for BistableStage0Error {
    fn from(value: EvidenceError) -> Self {
        Self::Evidence(value)
    }
}

/// Classify a completed positive-`D` coherent-switching sweep under Stage 0 v2 rules.
///
/// `INTERIOR_RESPONSE_PEAK` requires the finite-grid maximizer to lie on indices
/// 1..=5 (never an endpoint), exceed both immediate neighbors, and remain
/// uncertainty-separated from both scan edges at weight 2. Mechanism
/// compatibility additionally requires the peak coordinate within a factor of
/// two of `predicted_noise_intensity`.
#[must_use]
pub fn classify_bistable_stage0(
    points: &[BistableStage0GridPoint],
    predicted_noise_intensity: f64,
    controls_ok: bool,
    control_failure_detail: Option<&str>,
) -> BistableStage0Decision {
    if !controls_ok {
        return BistableStage0Decision::ControlFailure {
            detail: control_failure_detail
                .unwrap_or("required Stage 0 controls did not complete cleanly")
                .to_string(),
        };
    }
    if points.len() != BISTABLE_STAGE0_V2_GRID_FACTORS.len() {
        return BistableStage0Decision::ProtocolMismatch {
            detail: format!(
                "expected {} positive-D grid points, got {}",
                BISTABLE_STAGE0_V2_GRID_FACTORS.len(),
                points.len()
            ),
        };
    }
    if !(predicted_noise_intensity.is_finite() && predicted_noise_intensity > 0.0) {
        return BistableStage0Decision::ProtocolMismatch {
            detail: "predicted Kramers noise intensity must be finite and positive".to_string(),
        };
    }

    let Some(best) = best_sampled_summary(points) else {
        return BistableStage0Decision::ProtocolMismatch {
            detail: "unable to identify a finite-grid maximizer".to_string(),
        };
    };

    let responses: Vec<NoiseResponse> = points
        .iter()
        .map(BistableStage0GridPoint::to_noise_response)
        .collect();
    let edge = match detect_edge_separated_peak(&responses, BISTABLE_STAGE0_UNCERTAINTY_WEIGHT) {
        Ok(edge) => edge,
        Err(error) => {
            return BistableStage0Decision::ProtocolMismatch {
                detail: format!("edge-separation evidence failed closed: {error}"),
            };
        }
    };

    let interior_ok = best.grid_index >= 1
        && best.grid_index + 1 < points.len()
        && best.grid_index <= points.len().saturating_sub(2)
        && best.grid_index <= 5
        && points[best.grid_index].mean_coherent_amplitude
            > points[best.grid_index - 1].mean_coherent_amplitude
        && points[best.grid_index].mean_coherent_amplitude
            > points[best.grid_index + 1].mean_coherent_amplitude;

    let Some(peak) = edge else {
        return BistableStage0Decision::H0NotRejected {
            best_sampled: Some(best),
        };
    };

    // Edge helper may select a different interior index than the global
    // maximizer; INTERIOR_RESPONSE_PEAK requires the maximizer itself.
    if !interior_ok || peak.coordinate != best.noise_intensity {
        return BistableStage0Decision::H0NotRejected {
            best_sampled: Some(best),
        };
    }

    if within_kramers_factor(peak.coordinate, predicted_noise_intensity) {
        BistableStage0Decision::MechanismCompatibleInteriorPeak {
            peak,
            predicted_noise_intensity,
            best_sampled: best,
        }
    } else {
        BistableStage0Decision::InteriorPeakWithoutMechanismCompatibility {
            peak,
            predicted_noise_intensity,
            best_sampled: best,
        }
    }
}

fn within_kramers_factor(peak: f64, predicted: f64) -> bool {
    if !(peak.is_finite() && predicted.is_finite() && peak > 0.0 && predicted > 0.0) {
        return false;
    }
    let ratio = peak / predicted;
    ratio <= BISTABLE_STAGE0_KRAMERS_FACTOR_TOLERANCE
        && (1.0 / ratio) <= BISTABLE_STAGE0_KRAMERS_FACTOR_TOLERANCE
}

fn best_sampled_summary(points: &[BistableStage0GridPoint]) -> Option<BestSampledSummary> {
    points
        .iter()
        .enumerate()
        .filter(|(_, point)| point.mean_coherent_amplitude.is_finite())
        .max_by(|(_, left), (_, right)| {
            left.mean_coherent_amplitude
                .total_cmp(&right.mean_coherent_amplitude)
        })
        .map(|(grid_index, point)| BestSampledSummary {
            noise_intensity: point.noise_intensity,
            mean_coherent_amplitude: point.mean_coherent_amplitude,
            grid_index,
        })
}

/// Run the bistable Kramers Stage 0 v2 panel (outcome-blind execution path).
pub fn run_bistable_kramers_stage0(
    config: &BistableStage0Config,
) -> Result<BistableStage0Result, BistableStage0Error> {
    let protocol = BistableStage0V2::primary()?;
    if !protocol.model.forcing_is_subthreshold() {
        return Err(BistableStage0Error::Protocol {
            detail: "static subthreshold check failed; sweep must not be interpreted".to_string(),
        });
    }

    let steps_per_period =
        ((1.0 / protocol.model.forcing_frequency_hz) / protocol.run.dt).round() as usize;
    if steps_per_period == 0 {
        return Err(BistableStage0Error::Protocol {
            detail: "steps_per_period resolved to zero".to_string(),
        });
    }

    let total_periods = config.mode.total_periods();
    let burn_in_periods = config.mode.burn_in_periods();
    let run = LangevinRun::new(
        protocol.run.dt,
        total_periods * steps_per_period,
        burn_in_periods * steps_per_period,
    )?;
    let seed_count = config.mode.seed_count();
    if seed_count == 0 || seed_count > protocol.seeds.len() {
        return Err(BistableStage0Error::Protocol {
            detail: format!("invalid seed_count {seed_count}"),
        });
    }
    let seeds = protocol.seeds[..seed_count].to_vec();

    // Freeze analytic controls and derived grid before any trajectory.
    let analytic = BistableStage0AnalyticControls {
        barrier_height: protocol.model.barrier_height(),
        kramers_prefactor_rate: protocol.model.kramers_prefactor_rate(),
        predicted_noise_intensity: protocol.predicted_noise_intensity,
        positive_noise_grid: protocol.positive_noise_grid,
        grid_factors: BISTABLE_STAGE0_V2_GRID_FACTORS,
        protocol_blob_sha: BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA,
        forcing_is_subthreshold: protocol.model.forcing_is_subthreshold(),
    };

    let deterministic_control = run_deterministic_control(protocol.model, run, &seeds);
    let primary_grid =
        match sweep_noise_grid(protocol.model, run, &analytic.positive_noise_grid, &seeds) {
            Ok(grid) => grid,
            Err(error) => {
                return Ok(failed_panel(
                    config,
                    analytic,
                    total_periods,
                    burn_in_periods,
                    seeds,
                    deterministic_control,
                    Vec::new(),
                    Vec::new(),
                    format!("primary grid sweep failed: {error}"),
                ));
            }
        };

    let zero_forcing_model = DoubleWellLangevin::new(
        protocol.model.a,
        protocol.model.b,
        0.0,
        protocol.model.forcing_frequency_hz,
        protocol.model.phase_rad,
    )?;
    let zero_forcing_grid = match sweep_noise_grid(
        zero_forcing_model,
        run,
        &analytic.positive_noise_grid,
        &seeds,
    ) {
        Ok(grid) => grid,
        Err(error) => {
            return Ok(failed_panel(
                config,
                analytic,
                total_periods,
                burn_in_periods,
                seeds,
                deterministic_control,
                primary_grid,
                Vec::new(),
                format!("A=0 control sweep failed: {error}"),
            ));
        }
    };

    let (controls_ok, control_detail) = match &deterministic_control {
        BistableStage0DeterministicControl::SeedIndependent { .. } => (true, None),
        BistableStage0DeterministicControl::SeedDependent => (
            false,
            Some("D=0 trajectories are seed-dependent".to_string()),
        ),
        BistableStage0DeterministicControl::Failed { detail } => (false, Some(detail.clone())),
    };

    let primary_decision = classify_bistable_stage0(
        &primary_grid,
        analytic.predicted_noise_intensity,
        controls_ok,
        control_detail.as_deref(),
    );

    let falsification = match BistableStage0V2::falsification_kramers_control() {
        Ok(None) => BistableStage0FalsificationStatus::KramersControlInapplicable,
        Ok(Some(predicted)) => {
            // Guarded path: current frozen parameters yield None. If a future
            // pin makes a positive control available, execute with the same
            // multiplicative factors rather than inventing a new grid.
            let model = DoubleWellLangevin::new(
                protocol.model.a,
                protocol.model.b,
                protocol.model.forcing_amplitude,
                BISTABLE_STAGE0_V2_FALSIFICATION_FREQUENCY_HZ,
                protocol.model.phase_rad,
            )?;
            let steps_per_period = ((1.0 / model.forcing_frequency_hz) / run.dt).round() as usize;
            let falsification_run = LangevinRun::new(
                run.dt,
                total_periods * steps_per_period,
                burn_in_periods * steps_per_period,
            )?;
            let grid: [f64; 7] = BISTABLE_STAGE0_V2_GRID_FACTORS.map(|factor| factor * predicted);
            match sweep_noise_grid(model, falsification_run, &grid, &seeds) {
                Ok(points) => {
                    let decision = classify_bistable_stage0(
                        &points,
                        predicted,
                        controls_ok,
                        control_detail.as_deref(),
                    );
                    BistableStage0FalsificationStatus::Executed {
                        predicted_noise_intensity: predicted,
                        decision,
                    }
                }
                Err(error) => BistableStage0FalsificationStatus::Failed {
                    detail: error.to_string(),
                },
            }
        }
        Err(error) => BistableStage0FalsificationStatus::Failed {
            detail: error.to_string(),
        },
    };

    let summary = if config.mode.is_scientific() {
        BistableStage0Summary::Classified(primary_decision.clone())
    } else {
        BistableStage0Summary::NonScientificSmoke
    };

    Ok(BistableStage0Result {
        mode: config.mode,
        scientific_claim_permitted: config.mode.is_scientific(),
        analytic,
        total_periods,
        burn_in_periods,
        seed_count,
        seeds,
        deterministic_control,
        primary_grid,
        zero_forcing_grid,
        primary_decision,
        falsification,
        summary,
    })
}

#[allow(clippy::too_many_arguments)]
fn failed_panel(
    config: &BistableStage0Config,
    analytic: BistableStage0AnalyticControls,
    total_periods: usize,
    burn_in_periods: usize,
    seeds: Vec<u64>,
    deterministic_control: BistableStage0DeterministicControl,
    primary_grid: Vec<BistableStage0GridPoint>,
    zero_forcing_grid: Vec<BistableStage0GridPoint>,
    detail: String,
) -> BistableStage0Result {
    let primary_decision = BistableStage0Decision::ControlFailure {
        detail: detail.clone(),
    };
    let summary = if config.mode.is_scientific() {
        BistableStage0Summary::Incomplete { detail }
    } else {
        BistableStage0Summary::NonScientificSmoke
    };
    BistableStage0Result {
        mode: config.mode,
        scientific_claim_permitted: config.mode.is_scientific(),
        analytic,
        total_periods,
        burn_in_periods,
        seed_count: seeds.len(),
        seeds,
        deterministic_control,
        primary_grid,
        zero_forcing_grid,
        primary_decision,
        falsification: BistableStage0FalsificationStatus::KramersControlInapplicable,
        summary,
    }
}

fn run_deterministic_control(
    model: DoubleWellLangevin,
    run: LangevinRun,
    seeds: &[u64],
) -> BistableStage0DeterministicControl {
    if seeds.is_empty() {
        return BistableStage0DeterministicControl::Failed {
            detail: "no seeds supplied for D=0 control".to_string(),
        };
    }
    let mut reference = None;
    let mut amplitude = None;
    for &seed in seeds {
        match simulate_double_well(model, 0.0, -model.well_location(), run, seed) {
            Ok(trajectory) => {
                match coherent_switching_response(model, &trajectory, run.burn_in_steps) {
                    Ok(response) => {
                        if reference.is_none() {
                            reference = Some(trajectory.x.clone());
                            amplitude = Some(response.amplitude);
                        } else if reference.as_ref() != Some(&trajectory.x) {
                            return BistableStage0DeterministicControl::SeedDependent;
                        }
                    }
                    Err(error) => {
                        return BistableStage0DeterministicControl::Failed {
                            detail: error.to_string(),
                        };
                    }
                }
            }
            Err(error) => {
                return BistableStage0DeterministicControl::Failed {
                    detail: error.to_string(),
                };
            }
        }
    }
    BistableStage0DeterministicControl::SeedIndependent {
        amplitude: amplitude.unwrap_or(0.0),
    }
}

fn sweep_noise_grid(
    model: DoubleWellLangevin,
    run: LangevinRun,
    noise_intensities: &[f64],
    seeds: &[u64],
) -> Result<Vec<BistableStage0GridPoint>, LangevinError> {
    if seeds.is_empty() {
        return Err(LangevinError::NoSeeds);
    }
    let mut points = Vec::with_capacity(noise_intensities.len());
    for &noise_intensity in noise_intensities {
        if !noise_intensity.is_finite() {
            return Err(LangevinError::NonFinite("noise_intensity"));
        }
        if noise_intensity < 0.0 {
            return Err(LangevinError::Negative("noise_intensity"));
        }
        let mut amplitudes = Vec::with_capacity(seeds.len());
        for &seed in seeds {
            let trajectory =
                simulate_double_well(model, noise_intensity, -model.well_location(), run, seed)?;
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
        let standard_error = sample_stddev / (amplitudes.len() as f64).sqrt();
        points.push(BistableStage0GridPoint {
            noise_intensity,
            mean_coherent_amplitude: mean,
            sample_stddev,
            standard_error,
            replicates: amplitudes.len(),
            per_seed_amplitudes: amplitudes,
        });
    }
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn point(noise: f64, mean: f64, sd: f64) -> BistableStage0GridPoint {
        BistableStage0GridPoint {
            noise_intensity: noise,
            mean_coherent_amplitude: mean,
            sample_stddev: sd,
            standard_error: sd / 4.0,
            replicates: 16,
            per_seed_amplitudes: vec![mean; 16],
        }
    }

    fn interior_peak_grid() -> Vec<BistableStage0GridPoint> {
        // Indices 0..6; clear maximizer at index 3 with edge separation at w=2.
        vec![
            point(0.07, 0.10, 0.02),
            point(0.12, 0.20, 0.02),
            point(0.19, 0.35, 0.02),
            point(0.30, 0.55, 0.02),
            point(0.47, 0.34, 0.02),
            point(0.74, 0.22, 0.02),
            point(1.18, 0.12, 0.02),
        ]
    }

    #[test]
    fn classifies_mechanism_compatible_interior_peak() {
        let grid = interior_peak_grid();
        let predicted = 0.30;
        let decision = classify_bistable_stage0(&grid, predicted, true, None);
        match decision {
            BistableStage0Decision::MechanismCompatibleInteriorPeak {
                peak,
                predicted_noise_intensity,
                best_sampled,
            } => {
                assert_eq!(peak.coordinate, 0.30);
                assert_eq!(predicted_noise_intensity, predicted);
                assert_eq!(best_sampled.grid_index, 3);
            }
            other => panic!("unexpected decision: {other:?}"),
        }
    }

    #[test]
    fn retains_interior_peak_without_kramers_compatibility() {
        let grid = interior_peak_grid();
        let predicted = 0.05; // far from peak at 0.30
        let decision = classify_bistable_stage0(&grid, predicted, true, None);
        assert!(matches!(
            decision,
            BistableStage0Decision::InteriorPeakWithoutMechanismCompatibility { .. }
        ));
    }

    #[test]
    fn endpoint_maximizer_does_not_reject_h0() {
        let mut grid = interior_peak_grid();
        grid[0].mean_coherent_amplitude = 0.90;
        grid[0].per_seed_amplitudes = vec![0.90; 16];
        let decision = classify_bistable_stage0(&grid, 0.30, true, None);
        assert!(matches!(
            decision,
            BistableStage0Decision::H0NotRejected { .. }
        ));
    }

    #[test]
    fn control_failure_fail_closed() {
        let grid = interior_peak_grid();
        let decision = classify_bistable_stage0(&grid, 0.30, false, Some("D=0 seed-dependent"));
        assert_eq!(
            decision,
            BistableStage0Decision::ControlFailure {
                detail: "D=0 seed-dependent".to_string()
            }
        );
    }

    #[test]
    fn wrong_grid_length_is_protocol_mismatch() {
        let decision = classify_bistable_stage0(&interior_peak_grid()[..3], 0.3, true, None);
        assert!(matches!(
            decision,
            BistableStage0Decision::ProtocolMismatch { .. }
        ));
    }

    #[test]
    fn smoke_panel_is_deterministic_complete_and_non_scientific() {
        let config = BistableStage0Config::non_scientific_smoke();
        let first = run_bistable_kramers_stage0(&config).unwrap();
        let second = run_bistable_kramers_stage0(&config).unwrap();
        assert_eq!(first, second);
        assert!(!first.scientific_claim_permitted);
        assert_eq!(first.mode, BistableStage0Mode::NonScientificSmoke);
        assert_eq!(first.summary, BistableStage0Summary::NonScientificSmoke);
        assert_eq!(first.total_periods, BISTABLE_STAGE0_SMOKE_TOTAL_PERIODS);
        assert_eq!(first.burn_in_periods, BISTABLE_STAGE0_SMOKE_BURN_IN_PERIODS);
        assert_eq!(first.seed_count, BISTABLE_STAGE0_SMOKE_SEED_COUNT);
        assert_eq!(
            first.primary_grid.len(),
            BISTABLE_STAGE0_V2_GRID_FACTORS.len()
        );
        assert_eq!(
            first.zero_forcing_grid.len(),
            BISTABLE_STAGE0_V2_GRID_FACTORS.len()
        );
        assert!(matches!(
            first.deterministic_control,
            BistableStage0DeterministicControl::SeedIndependent { .. }
        ));
        assert_eq!(
            first.falsification,
            BistableStage0FalsificationStatus::KramersControlInapplicable
        );
        assert!((first.analytic.predicted_noise_intensity - 0.2958709422522067).abs() < 1e-12);
        assert!(first.analytic.forcing_is_subthreshold);
    }

    #[test]
    fn scientific_config_exposes_preregistered_budget_only() {
        let config = BistableStage0Config::scientific();
        assert!(config.mode.is_scientific());
        assert_eq!(config.mode.total_periods(), BISTABLE_STAGE0_TOTAL_PERIODS);
        assert_eq!(
            config.mode.burn_in_periods(),
            BISTABLE_STAGE0_BURN_IN_PERIODS
        );
        assert_eq!(config.mode.seed_count(), BISTABLE_STAGE0_V2_SEEDS.len());
    }
}
