//! Fail-closed decision gate for the frozen FHN coherence-resonance Stage 0.
//!
//! This module does not run the stochastic sweep. It only checks that a supplied
//! calibration was produced on the preregistered coordinate grid with the
//! preregistered replicate count before applying the frozen acceptance interval.

use crate::fhn::{CoherenceCalibration, CoherenceResponse};
use crate::preregistered::{
    FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL, FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES,
    FHN_COHERENCE_STAGE0_SEEDS, FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
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
        if !valid_response(response, expected_noise) {
            return FhnStage0Decision::ProtocolMismatch;
        }
    }

    let responses = &calibration.responses;
    let (minimum_index, minimum_response) = responses[1..responses.len() - 1]
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| left.mean_cv.total_cmp(&right.mean_cv))
        .map(|(offset, response)| (offset + 1, response))
        .expect("preregistered FHN grid has interior points");

    let minimum_se = standard_error(minimum_response);
    let conservative_minimum =
        minimum_response.mean_cv + FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * minimum_se;
    let conservative_left = responses[0].mean_cv
        - FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * standard_error(&responses[0]);
    let conservative_right = responses[responses.len() - 1].mean_cv
        - FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * standard_error(&responses[responses.len() - 1]);
    let conservative_edge = conservative_left.min(conservative_right);
    let conservative_separation = conservative_edge - conservative_minimum;
    let strict_local_minimum = minimum_response.mean_cv < responses[minimum_index - 1].mean_cv
        && minimum_response.mean_cv < responses[minimum_index + 1].mean_cv;
    let expected_minimum = strict_local_minimum && conservative_separation > 0.0;

    let Some(minimum) = calibration.interior_minimum else {
        return if expected_minimum {
            FhnStage0Decision::ProtocolMismatch
        } else {
            FhnStage0Decision::NoInteriorMinimum
        };
    };
    if !expected_minimum
        || minimum.noise_amplitude != minimum_response.noise_amplitude
        || !approximately_equal(minimum.mean_cv, minimum_response.mean_cv)
        || !approximately_equal(minimum.standard_error, minimum_se)
        || !approximately_equal(minimum.conservative_minimum, conservative_minimum)
        || !approximately_equal(minimum.conservative_edge, conservative_edge)
        || !approximately_equal(minimum.conservative_separation, conservative_separation)
        || minimum.uncertainty_weight != FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT
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

fn valid_response(response: &CoherenceResponse, expected_noise: f64) -> bool {
    response.noise_amplitude == expected_noise
        && response.mean_cv.is_finite()
        && response.mean_cv >= 0.0
        && response.sample_stddev_cv.is_finite()
        && response.sample_stddev_cv >= 0.0
        && response.mean_interval.is_finite()
        && response.mean_interval > 0.0
        && response.replicates == FHN_COHERENCE_STAGE0_SEEDS.len()
}

fn standard_error(response: &CoherenceResponse) -> f64 {
    response.sample_stddev_cv / (response.replicates as f64).sqrt()
}

fn approximately_equal(left: f64, right: f64) -> bool {
    if !left.is_finite() || !right.is_finite() {
        return false;
    }
    let scale = 1.0_f64.max(left.abs()).max(right.abs());
    (left - right).abs() <= 1e-12 * scale
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fhn::{CoherenceMinimum, CoherenceResponse};

    fn response(noise_amplitude: f64, mean_cv: f64) -> CoherenceResponse {
        CoherenceResponse {
            noise_amplitude,
            mean_cv,
            sample_stddev_cv: 0.1,
            mean_interval: 1.0,
            replicates: FHN_COHERENCE_STAGE0_SEEDS.len(),
        }
    }

    fn calibration(minimum_noise: Option<f64>) -> CoherenceCalibration {
        let responses: Vec<_> = FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
            .iter()
            .copied()
            .map(|noise_amplitude| {
                let mean_cv = if minimum_noise == Some(noise_amplitude) {
                    0.5
                } else {
                    1.0
                };
                response(noise_amplitude, mean_cv)
            })
            .collect();
        let interior_minimum = minimum_noise.map(|noise_amplitude| {
            let selected = responses
                .iter()
                .find(|response| response.noise_amplitude == noise_amplitude)
                .unwrap();
            let selected_standard_error = standard_error(selected);
            let conservative_minimum = selected.mean_cv
                + FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * selected_standard_error;
            let conservative_edge = (responses[0].mean_cv
                - FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT * standard_error(&responses[0]))
            .min(
                responses[responses.len() - 1].mean_cv
                    - FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT
                        * standard_error(&responses[responses.len() - 1]),
            );
            CoherenceMinimum {
                noise_amplitude,
                mean_cv: selected.mean_cv,
                standard_error: selected_standard_error,
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
    fn rejects_non_finite_response_summary() {
        let mut malformed = calibration(Some(0.075));
        malformed.responses[0].mean_cv = f64::NAN;
        assert_eq!(
            classify_fhn_stage0(&malformed),
            FhnStage0Decision::ProtocolMismatch
        );
    }

    #[test]
    fn rejects_fabricated_minimum_derived_fields() {
        let mut fabricated = calibration(Some(0.075));
        fabricated
            .interior_minimum
            .as_mut()
            .unwrap()
            .conservative_separation += 0.1;
        assert_eq!(
            classify_fhn_stage0(&fabricated),
            FhnStage0Decision::ProtocolMismatch
        );
    }

    #[test]
    fn rejects_missing_positive_minimum_evidence() {
        let mut incomplete = calibration(Some(0.075));
        incomplete.interior_minimum = None;
        assert_eq!(
            classify_fhn_stage0(&incomplete),
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
