//! Deterministic gate for the preregistered FHN observation-horizon control.
//!
//! This module does not run trajectories or change any frozen parameter. It
//! validates completed Stage-0 calibration summaries at the three preregistered
//! horizons and classifies only the robustness decision declared in
//! `docs/research/fhn-observation-horizon-stage0.md`.

use crate::fhn::CoherenceCalibration;
use crate::fhn_stage0::{classify_fhn_stage0, FhnStage0Decision};

pub const FHN_HORIZON_STAGE0_STEPS: [usize; 3] = [40_000, 80_000, 160_000];

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
            let conservative_minimum = selected.mean_cv
                + FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * standard_error;
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
}
