//! Conservative evidence helpers for response sweeps.
//!
//! These helpers deliberately avoid calling every sampled maximum a
//! "resonance". They quantify narrower statements that can be checked from a
//! finite experiment, such as whether the best interior response remains above
//! both scanned boundary regimes after an explicit uncertainty penalty.

use crate::langevin::NoiseResponse;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Evidence that an interior mean response is separated from both scan edges.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EdgeSeparatedPeak {
    /// Noise intensity of the strongest sampled interior mean response.
    pub coordinate: f64,
    /// Mean response at that coordinate.
    pub response: f64,
    /// Standard error of the selected interior mean.
    pub standard_error: f64,
    /// Lower penalized value `mean - w * SE` at the interior maximum.
    pub conservative_peak: f64,
    /// Larger of the two edge upper values `mean + w * SE`.
    pub conservative_edge: f64,
    /// `conservative_peak - conservative_edge`; positive by construction for
    /// returned evidence.
    pub conservative_prominence: f64,
    /// Multiplier `w` used for the standard-error penalty.
    pub uncertainty_weight: f64,
}

/// Invalid response evidence request.
#[derive(Debug, Clone, PartialEq)]
pub enum EvidenceError {
    /// At least three points are required to have an interior candidate.
    SweepTooShort,
    /// The uncertainty multiplier was negative or non-finite.
    InvalidUncertaintyWeight,
    /// A response contains a non-finite coordinate or statistic.
    NonFiniteResponse { index: usize },
    /// A response was produced without any replicate.
    ZeroReplicates { index: usize },
}

impl Display for EvidenceError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SweepTooShort => {
                formatter.write_str("response sweep needs at least three points")
            }
            Self::InvalidUncertaintyWeight => {
                formatter.write_str("uncertainty weight must be finite and non-negative")
            }
            Self::NonFiniteResponse { index } => {
                write!(formatter, "response {index} contains a non-finite value")
            }
            Self::ZeroReplicates { index } => {
                write!(formatter, "response {index} has zero replicates")
            }
        }
    }
}

impl Error for EvidenceError {}

/// Test whether the strongest sampled interior mean response is conservatively
/// separated from both ends of a noise-intensity scan.
///
/// For response `i`, `SE_i = sample_stddev_i / sqrt(n_i)`. The strongest
/// interior mean is retained only when
///
/// `mean_peak - w*SE_peak > max(mean_left + w*SE_left, mean_right + w*SE_right)`.
///
/// This is intentionally described as an uncertainty-penalized heuristic, not
/// a formal confidence interval or a proof of a resonance mechanism. It is
/// useful for broad stochastic-resonance peaks where immediate-neighbor
/// prominence can be arbitrarily small merely because the parameter grid is
/// dense near the maximum.
pub fn detect_edge_separated_peak(
    responses: &[NoiseResponse],
    uncertainty_weight: f64,
) -> Result<Option<EdgeSeparatedPeak>, EvidenceError> {
    if responses.len() < 3 {
        return Err(EvidenceError::SweepTooShort);
    }
    if !uncertainty_weight.is_finite() || uncertainty_weight < 0.0 {
        return Err(EvidenceError::InvalidUncertaintyWeight);
    }

    for (index, response) in responses.iter().enumerate() {
        if !response.noise_intensity.is_finite()
            || !response.mean_coherent_amplitude.is_finite()
            || !response.sample_stddev.is_finite()
        {
            return Err(EvidenceError::NonFiniteResponse { index });
        }
        if response.replicates == 0 {
            return Err(EvidenceError::ZeroReplicates { index });
        }
    }

    let (interior_index, peak) = responses[1..responses.len() - 1]
        .iter()
        .enumerate()
        .max_by(|(_, left), (_, right)| {
            left.mean_coherent_amplitude
                .total_cmp(&right.mean_coherent_amplitude)
        })
        .map(|(offset, response)| (offset + 1, response))
        .expect("slice length is guaranteed non-zero");

    let peak_se = standard_error(peak);
    let conservative_peak = peak.mean_coherent_amplitude - uncertainty_weight * peak_se;
    let left = &responses[0];
    let right = &responses[responses.len() - 1];
    let conservative_left =
        left.mean_coherent_amplitude + uncertainty_weight * standard_error(left);
    let conservative_right =
        right.mean_coherent_amplitude + uncertainty_weight * standard_error(right);
    let conservative_edge = conservative_left.max(conservative_right);
    let conservative_prominence = conservative_peak - conservative_edge;

    if conservative_prominence <= 0.0 {
        return Ok(None);
    }

    // Keep the selected index explicit in the implementation so any future
    // extension that reports scan provenance cannot accidentally reinterpret a
    // boundary point as an interior candidate.
    debug_assert!(interior_index > 0 && interior_index + 1 < responses.len());

    Ok(Some(EdgeSeparatedPeak {
        coordinate: peak.noise_intensity,
        response: peak.mean_coherent_amplitude,
        standard_error: peak_se,
        conservative_peak,
        conservative_edge,
        conservative_prominence,
        uncertainty_weight,
    }))
}

fn standard_error(response: &NoiseResponse) -> f64 {
    response.sample_stddev / (response.replicates as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn response(noise: f64, mean: f64, sd: f64, replicates: usize) -> NoiseResponse {
        NoiseResponse {
            noise_intensity: noise,
            mean_coherent_amplitude: mean,
            sample_stddev: sd,
            replicates,
        }
    }

    #[test]
    fn broad_peak_can_be_strong_even_when_local_top_is_flat() {
        let responses = [
            response(0.05, 0.07, 0.04, 8),
            response(0.20, 0.426, 0.082, 8),
            response(0.30, 0.430, 0.087, 8),
            response(0.45, 0.372, 0.107, 8),
            response(1.00, 0.262, 0.083, 8),
        ];
        let evidence = detect_edge_separated_peak(&responses, 2.0)
            .unwrap()
            .expect("interior response should remain separated from both edges");
        assert_eq!(evidence.coordinate, 0.30);
        assert!(evidence.conservative_prominence > 0.0);
    }

    #[test]
    fn monotone_or_uncertain_sweep_is_not_promoted() {
        let responses = [
            response(0.1, 0.10, 0.10, 4),
            response(0.2, 0.15, 0.10, 4),
            response(0.3, 0.20, 0.10, 4),
        ];
        assert_eq!(detect_edge_separated_peak(&responses, 2.0).unwrap(), None);
    }

    #[test]
    fn malformed_evidence_requests_fail_closed() {
        let short = [response(0.1, 0.1, 0.0, 1), response(0.2, 0.2, 0.0, 1)];
        assert_eq!(
            detect_edge_separated_peak(&short, 1.0),
            Err(EvidenceError::SweepTooShort)
        );
        assert_eq!(
            detect_edge_separated_peak(
                &[
                    response(0.1, 0.1, 0.0, 1),
                    response(0.2, 0.2, 0.0, 0),
                    response(0.3, 0.1, 0.0, 1),
                ],
                1.0,
            ),
            Err(EvidenceError::ZeroReplicates { index: 1 })
        );
    }
}
