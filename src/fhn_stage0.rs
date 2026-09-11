//! Fail-closed decision gate for the frozen FHN coherence-resonance Stage 0.
//!
//! This module does not run the stochastic sweep. It only checks that a supplied
//! calibration was produced on the preregistered coordinate grid with the
//! preregistered replicate count before applying the frozen acceptance interval.

use crate::fhn::CoherenceCalibration;
use crate::preregistered::{
    FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL, FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES,
    FHN_COHERENCE_STAGE0_SEEDS,
};

/// Outcome of applying the frozen Stage 0 decision gate to a calibration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FhnStage0Decision {
    /// The supplied calibration does not match the preregistered grid/replicate contract.
    ProtocolMismatch,
    /// The calibrated responses contain no uncertainty-separated interior minimum.
    NoInteriorMinimum,
    /// An interior minimum exists but lies outside the preregistered acceptance interval.
    OutsideAcceptance { noise_amplitude: f64 },
    /// An interior minimum exists and lies inside the preregistered acceptance interval.
    Accepted { noise_amplitude: f64 },
}

/// Apply the frozen Stage 0 decision rule without running or retuning an experiment.
#[must_use]
pub fn classify_fhn_stage0(calibration: &CoherenceCalibration) -> FhnStage0Decision {
    if calibration.responses.len() != FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES.len() {
        return FhnStage0Decision::ProtocolMismatch;
    }

    for (response, expected_noise) in calibration
        .responses
        .iter()
        .zip(FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES)
    {
        if response.noise_amplitude != expected_noise
            || response.replicates != FHN_COHERENCE_STAGE0_SEEDS.len()
        {
            return FhnStage0Decision::ProtocolMismatch;
        }
    }

    let Some(minimum) = calibration.interior_minimum else {
        return FhnStage0Decision::NoInteriorMinimum;
    };
    if !minimum.noise_amplitude.is_finite()
        || !FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
            [1..FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES.len() - 1]
            .contains(&minimum.noise_amplitude)
    {
        return FhnStage0Decision::ProtocolMismatch;
    }

    let [lower, upper] = FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL;
    if minimum.noise_amplitude < lower || minimum.noise_amplitude > upper {
        FhnStage0Decision::OutsideAcceptance {
            noise_amplitude: minimum.noise_amplitude,
        }
    } else {
        FhnStage0Decision::Accepted {
            noise_amplitude: minimum.noise_amplitude,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fhn::{CoherenceMinimum, CoherenceResponse};

    fn response(noise_amplitude: f64) -> CoherenceResponse {
        CoherenceResponse {
            noise_amplitude,
            mean_cv: 1.0,
            sample_stddev_cv: 0.1,
            mean_interval: 1.0,
            replicates: FHN_COHERENCE_STAGE0_SEEDS.len(),
        }
    }

    fn calibration(minimum_noise: Option<f64>) -> CoherenceCalibration {
        let responses = FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
            .iter()
            .copied()
            .map(response)
            .collect();
        let interior_minimum = minimum_noise.map(|noise_amplitude| CoherenceMinimum {
            noise_amplitude,
            mean_cv: 0.5,
            standard_error: 0.01,
            conservative_minimum: 0.52,
            conservative_edge: 0.8,
            conservative_separation: 0.28,
            uncertainty_weight: 2.0,
        });
        CoherenceCalibration {
            responses,
            interior_minimum,
        }
    }

    #[test]
    fn accepts_only_preregistered_interior_coordinate_inside_acceptance_interval() {
        assert_eq!(
            classify_fhn_stage0(&calibration(Some(0.075))),
            FhnStage0Decision::Accepted {
                noise_amplitude: 0.075
            }
        );
    }

    #[test]
    fn preserves_negative_result_when_no_interior_minimum_exists() {
        assert_eq!(
            classify_fhn_stage0(&calibration(None)),
            FhnStage0Decision::NoInteriorMinimum
        );
    }

    #[test]
    fn rejects_protocol_drift_before_acceptance_logic() {
        let mut drifted = calibration(Some(0.075));
        drifted.responses[2].noise_amplitude = 0.051;
        assert_eq!(
            classify_fhn_stage0(&drifted),
            FhnStage0Decision::ProtocolMismatch
        );
    }

    #[test]
    fn reports_preregistered_grid_point_outside_acceptance_interval() {
        assert_eq!(
            classify_fhn_stage0(&calibration(Some(0.40))),
            FhnStage0Decision::OutsideAcceptance {
                noise_amplitude: 0.40
            }
        );
    }
}
