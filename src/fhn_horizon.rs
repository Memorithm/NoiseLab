//! Deterministic gate and executable panel for the preregistered FHN
//! observation-horizon control.
//!
//! The classifier does not run trajectories or change any frozen parameter. The
//! panel runner materializes preregistered Stage 0 factors, executes the three
//! horizons fail-closed, and applies [`classify_fhn_horizon_robustness`] only
//! when all scientific calibrations complete. Smoke mode uses abbreviated
//! horizons and must not be recorded as Stage 0 evidence.
//!
//! Protocol: `docs/research/fhn-observation-horizon-stage0.md`.

use crate::fhn::{calibrate_coherence_resonance, CoherenceCalibration, FhnError, FhnRun};
use crate::fhn_stage0::{classify_fhn_stage0, FhnStage0Decision};
use crate::preregistered::{
    FhnCoherenceStage0, FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES, FHN_COHERENCE_STAGE0_SEEDS,
    FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Preregistered total-step budgets for the scientific horizon control.
pub const FHN_HORIZON_STAGE0_STEPS: [usize; 3] = [40_000, 80_000, 160_000];

/// Parent-protocol burn-in retained for every scientific horizon.
pub const FHN_HORIZON_STAGE0_BURN_IN_STEPS: usize = 10_000;

/// Abbreviated total-step budgets for non-scientific smoke only.
///
/// These are **not** the preregistered horizons. Smoke output must not be
/// recorded as a Stage 0 robustness result.
pub const FHN_HORIZON_SMOKE_STEPS: [usize; 3] = [2_000, 4_000, 8_000];

/// Burn-in used only by the non-scientific smoke path.
pub const FHN_HORIZON_SMOKE_BURN_IN_STEPS: usize = 500;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FhnHorizonRobustness {
    ProtocolMismatch,
    NotRobust,
    StableNoInteriorMinimum,
    StableAccepted { noise_amplitudes: [f64; 3] },
}

#[must_use]
pub fn classify_fhn_horizon_robustness(
    observations: &[(usize, CoherenceCalibration)],
) -> FhnHorizonRobustness {
    if observations.len() != FHN_HORIZON_STAGE0_STEPS.len() {
        return FhnHorizonRobustness::ProtocolMismatch;
    }

    let mut accepted = [0.0; 3];
    let mut accepted_count = 0usize;
    let mut no_minimum_count = 0usize;

    for (index, ((steps, calibration), expected_steps)) in observations
        .iter()
        .zip(FHN_HORIZON_STAGE0_STEPS)
        .enumerate()
    {
        if *steps != expected_steps {
            return FhnHorizonRobustness::ProtocolMismatch;
        }

        match classify_fhn_stage0(calibration) {
            FhnStage0Decision::Accepted { noise_amplitude } => {
                accepted[index] = noise_amplitude;
                accepted_count += 1;
            }
            FhnStage0Decision::NoInteriorMinimum => {
                no_minimum_count += 1;
            }
            FhnStage0Decision::OutsideAcceptance { .. } => {
                return FhnHorizonRobustness::NotRobust;
            }
            FhnStage0Decision::ProtocolMismatch => {
                return FhnHorizonRobustness::ProtocolMismatch;
            }
        }
    }

    if accepted_count == FHN_HORIZON_STAGE0_STEPS.len() {
        return FhnHorizonRobustness::StableAccepted {
            noise_amplitudes: accepted,
        };
    }
    if no_minimum_count == FHN_HORIZON_STAGE0_STEPS.len() {
        return FhnHorizonRobustness::StableNoInteriorMinimum;
    }

    FhnHorizonRobustness::NotRobust
}

/// Whether this panel invocation may claim scientific horizon evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FhnHorizonStage0Mode {
    /// Exact preregistered horizons, grid, seeds, and acceptance contract.
    Scientific,
    /// Explicitly non-scientific CI / developer smoke (abbreviated horizons).
    NonScientificSmoke,
}

impl FhnHorizonStage0Mode {
    #[must_use]
    pub const fn is_scientific(self) -> bool {
        matches!(self, Self::Scientific)
    }

    #[must_use]
    pub const fn horizons(self) -> [usize; 3] {
        match self {
            Self::Scientific => FHN_HORIZON_STAGE0_STEPS,
            Self::NonScientificSmoke => FHN_HORIZON_SMOKE_STEPS,
        }
    }

    #[must_use]
    pub const fn burn_in_steps(self) -> usize {
        match self {
            Self::Scientific => FHN_HORIZON_STAGE0_BURN_IN_STEPS,
            Self::NonScientificSmoke => FHN_HORIZON_SMOKE_BURN_IN_STEPS,
        }
    }
}

/// Configuration frozen before inspecting horizon-robustness outcomes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FhnHorizonStage0Config {
    pub mode: FhnHorizonStage0Mode,
}

impl Default for FhnHorizonStage0Config {
    fn default() -> Self {
        Self::scientific()
    }
}

impl FhnHorizonStage0Config {
    /// Preregistered scientific panel (exact horizons and Stage 0 factors).
    #[must_use]
    pub const fn scientific() -> Self {
        Self {
            mode: FhnHorizonStage0Mode::Scientific,
        }
    }

    /// Explicitly labeled non-scientific smoke configuration.
    #[must_use]
    pub const fn non_scientific_smoke() -> Self {
        Self {
            mode: FhnHorizonStage0Mode::NonScientificSmoke,
        }
    }
}

/// Per-horizon outcome after fail-closed calibration + Stage 0 classification.
#[derive(Debug, Clone, PartialEq)]
pub enum FhnHorizonStage0HorizonDecision {
    /// Calibration completed; parent Stage 0 gate applied.
    Stage0(FhnStage0Decision),
    /// A replicate produced fewer than `min_spikes` after burn-in.
    InsufficientSpikes {
        noise_amplitude: f64,
        seed: u64,
        observed: usize,
        required: usize,
    },
    /// Any other calibration / protocol construction failure.
    ProtocolFailure { detail: String },
}

/// One horizon's recorded observation.
#[derive(Debug, Clone, PartialEq)]
pub struct FhnHorizonStage0Observation {
    pub steps: usize,
    pub decision: FhnHorizonStage0HorizonDecision,
    /// Present only when calibration completed without error.
    pub calibration: Option<CoherenceCalibration>,
}

/// Panel-level summary after all three horizons are attempted.
#[derive(Debug, Clone, PartialEq)]
pub enum FhnHorizonStage0Summary {
    /// All three calibrations completed; [`classify_fhn_horizon_robustness`] applied.
    Classified(FhnHorizonRobustness),
    /// All three horizons failed closed on insufficient spikes (same decision class).
    StableInsufficientSpikes,
    /// Horizons did not all complete calibrations, or classes mixed before the gate.
    Incomplete,
    /// Smoke path: robustness gate is not applied to abbreviated horizons.
    NonScientificSmoke,
}

/// Complete FHN observation-horizon Stage 0 panel result.
#[derive(Debug, Clone, PartialEq)]
pub struct FhnHorizonStage0Result {
    pub mode: FhnHorizonStage0Mode,
    pub scientific_claim_permitted: bool,
    pub observations: Vec<FhnHorizonStage0Observation>,
    pub summary: FhnHorizonStage0Summary,
}

/// Failures that prevent even starting the preregistered panel contract.
#[derive(Debug, Clone, PartialEq)]
pub enum FhnHorizonStage0Error {
    Protocol(FhnError),
    HorizonCount { got: usize, expected: usize },
}

impl Display for FhnHorizonStage0Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Protocol(error) => write!(f, "FHN horizon Stage 0 protocol error: {error}"),
            Self::HorizonCount { got, expected } => {
                write!(f, "expected {expected} horizons, got {got}")
            }
        }
    }
}

impl Error for FhnHorizonStage0Error {}

impl From<FhnError> for FhnHorizonStage0Error {
    fn from(value: FhnError) -> Self {
        Self::Protocol(value)
    }
}

/// Run the FHN observation-horizon Stage 0 panel (outcome-blind execution path).
///
/// Uses preregistered FHN Stage 0 factors (noise grid, seeds, uncertainty weight,
/// acceptance interval, model, `dt`, spike threshold, `min_spikes`) unchanged.
/// The only mode-dependent knobs are total steps and burn-in (smoke only).
///
/// Fail-closed: insufficient spikes and other calibration errors are recorded per
/// horizon and never retuned. Scientific robustness classification is applied only
/// when all three preregistered horizons produce complete calibrations.
pub fn run_fhn_horizon_stage0(
    config: &FhnHorizonStage0Config,
) -> Result<FhnHorizonStage0Result, FhnHorizonStage0Error> {
    let protocol = FhnCoherenceStage0::materialize()?;
    let horizons = config.mode.horizons();
    let burn_in_steps = config.mode.burn_in_steps();

    if horizons.len() != FHN_HORIZON_STAGE0_STEPS.len() {
        return Err(FhnHorizonStage0Error::HorizonCount {
            got: horizons.len(),
            expected: FHN_HORIZON_STAGE0_STEPS.len(),
        });
    }

    // Freeze factors before any horizon outcome is inspected.
    let noise_amplitudes = FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES;
    let seeds = FHN_COHERENCE_STAGE0_SEEDS;
    let uncertainty_weight = FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT;
    debug_assert_eq!(protocol.noise_amplitudes, noise_amplitudes);
    debug_assert_eq!(protocol.seeds, seeds);
    debug_assert_eq!(protocol.uncertainty_weight, uncertainty_weight);

    let mut observations = Vec::with_capacity(horizons.len());
    for steps in horizons {
        let observation = run_one_horizon(&protocol, steps, burn_in_steps)?;
        observations.push(observation);
    }

    let summary = summarize_panel(config.mode, &observations);
    Ok(FhnHorizonStage0Result {
        mode: config.mode,
        scientific_claim_permitted: config.mode.is_scientific(),
        observations,
        summary,
    })
}

fn run_one_horizon(
    protocol: &FhnCoherenceStage0,
    steps: usize,
    burn_in_steps: usize,
) -> Result<FhnHorizonStage0Observation, FhnHorizonStage0Error> {
    let run = match FhnRun::new(
        protocol.run.dt,
        steps,
        burn_in_steps,
        protocol.run.spike_threshold,
        protocol.run.min_spikes,
    ) {
        Ok(run) => run,
        Err(error) => {
            return Ok(FhnHorizonStage0Observation {
                steps,
                decision: FhnHorizonStage0HorizonDecision::ProtocolFailure {
                    detail: error.to_string(),
                },
                calibration: None,
            });
        }
    };

    match calibrate_coherence_resonance(
        protocol.model,
        run,
        &protocol.noise_amplitudes,
        &protocol.seeds,
        protocol.uncertainty_weight,
    ) {
        Ok(calibration) => {
            let decision =
                FhnHorizonStage0HorizonDecision::Stage0(classify_fhn_stage0(&calibration));
            Ok(FhnHorizonStage0Observation {
                steps,
                decision,
                calibration: Some(calibration),
            })
        }
        Err(FhnError::InsufficientSpikes {
            noise_amplitude,
            seed,
            observed,
            required,
        }) => Ok(FhnHorizonStage0Observation {
            steps,
            decision: FhnHorizonStage0HorizonDecision::InsufficientSpikes {
                noise_amplitude,
                seed,
                observed,
                required,
            },
            calibration: None,
        }),
        Err(error) => Ok(FhnHorizonStage0Observation {
            steps,
            decision: FhnHorizonStage0HorizonDecision::ProtocolFailure {
                detail: error.to_string(),
            },
            calibration: None,
        }),
    }
}

fn summarize_panel(
    mode: FhnHorizonStage0Mode,
    observations: &[FhnHorizonStage0Observation],
) -> FhnHorizonStage0Summary {
    if !mode.is_scientific() {
        return FhnHorizonStage0Summary::NonScientificSmoke;
    }

    if observations.iter().all(|obs| {
        matches!(
            obs.decision,
            FhnHorizonStage0HorizonDecision::InsufficientSpikes { .. }
        )
    }) {
        return FhnHorizonStage0Summary::StableInsufficientSpikes;
    }

    let mut pairs: Vec<(usize, CoherenceCalibration)> = Vec::with_capacity(observations.len());
    for observation in observations {
        match &observation.calibration {
            Some(calibration) => pairs.push((observation.steps, calibration.clone())),
            None => return FhnHorizonStage0Summary::Incomplete,
        }
    }

    FhnHorizonStage0Summary::Classified(classify_fhn_horizon_robustness(&pairs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fhn::{CoherenceMinimum, CoherenceResponse};
    use crate::preregistered::{
        FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES, FHN_COHERENCE_STAGE0_SEEDS,
        FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
    };

    fn calibration(minimum_noise: Option<f64>) -> CoherenceCalibration {
        let responses: Vec<_> = FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
            .iter()
            .copied()
            .map(|noise_amplitude| CoherenceResponse {
                noise_amplitude,
                mean_cv: if minimum_noise == Some(noise_amplitude) {
                    0.5
                } else {
                    1.0
                },
                sample_stddev_cv: 0.1,
                mean_interval: 1.0,
                replicates: FHN_COHERENCE_STAGE0_SEEDS.len(),
            })
            .collect();

        let interior_minimum = minimum_noise.map(|noise_amplitude| {
            let selected = responses
                .iter()
                .find(|response| response.noise_amplitude == noise_amplitude)
                .unwrap();
            let standard_error = selected.sample_stddev_cv / (selected.replicates as f64).sqrt();
            let conservative_minimum =
                selected.mean_cv + FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * standard_error;
            let edge_standard_error =
                responses[0].sample_stddev_cv / (responses[0].replicates as f64).sqrt();
            let conservative_edge = responses[0].mean_cv
                - FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * edge_standard_error;
            CoherenceMinimum {
                noise_amplitude,
                mean_cv: selected.mean_cv,
                standard_error,
                conservative_minimum,
                conservative_edge,
                conservative_separation: conservative_edge - conservative_minimum,
                uncertainty_weight: FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
            }
        });

        CoherenceCalibration {
            responses,
            interior_minimum,
        }
    }

    #[test]
    fn accepts_stable_interior_minimum_across_all_preregistered_horizons() {
        let observations = FHN_HORIZON_STAGE0_STEPS.map(|steps| (steps, calibration(Some(0.075))));
        assert_eq!(
            classify_fhn_horizon_robustness(&observations),
            FhnHorizonRobustness::StableAccepted {
                noise_amplitudes: [0.075, 0.075, 0.075]
            }
        );
    }

    #[test]
    fn preserves_stable_negative_result() {
        let observations = FHN_HORIZON_STAGE0_STEPS.map(|steps| (steps, calibration(None)));
        assert_eq!(
            classify_fhn_horizon_robustness(&observations),
            FhnHorizonRobustness::StableNoInteriorMinimum
        );
    }

    #[test]
    fn rejects_mixed_decision_classes_as_not_robust() {
        let observations = [
            (40_000, calibration(Some(0.075))),
            (80_000, calibration(None)),
            (160_000, calibration(Some(0.075))),
        ];
        assert_eq!(
            classify_fhn_horizon_robustness(&observations),
            FhnHorizonRobustness::NotRobust
        );
    }

    #[test]
    fn rejects_horizon_drift_before_reading_outcomes() {
        let observations = [
            (39_999, calibration(Some(0.075))),
            (80_000, calibration(Some(0.075))),
            (160_000, calibration(Some(0.075))),
        ];
        assert_eq!(
            classify_fhn_horizon_robustness(&observations),
            FhnHorizonRobustness::ProtocolMismatch
        );
    }

    #[test]
    fn smoke_panel_is_deterministic_complete_and_non_scientific() {
        let config = FhnHorizonStage0Config::non_scientific_smoke();
        let first = run_fhn_horizon_stage0(&config).unwrap();
        let second = run_fhn_horizon_stage0(&config).unwrap();
        assert_eq!(first, second);
        assert!(!first.scientific_claim_permitted);
        assert_eq!(first.mode, FhnHorizonStage0Mode::NonScientificSmoke);
        assert_eq!(first.summary, FhnHorizonStage0Summary::NonScientificSmoke);
        assert_eq!(first.observations.len(), FHN_HORIZON_SMOKE_STEPS.len());
        for (observation, expected_steps) in first.observations.iter().zip(FHN_HORIZON_SMOKE_STEPS)
        {
            assert_eq!(observation.steps, expected_steps);
        }
    }

    #[test]
    fn scientific_config_exposes_preregistered_horizons_only() {
        let config = FhnHorizonStage0Config::scientific();
        assert!(config.mode.is_scientific());
        assert_eq!(config.mode.horizons(), FHN_HORIZON_STAGE0_STEPS);
        assert_eq!(
            config.mode.burn_in_steps(),
            FHN_HORIZON_STAGE0_BURN_IN_STEPS
        );
    }
}
