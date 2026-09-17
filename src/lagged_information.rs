//! Descriptive non-circular lagged mutual-information diagnostics.
//!
//! The full real-valued observation is quantized once under a declared
//! equal-width bin count, then each requested lag is evaluated only on its
//! overlapping samples. The scan is descriptive: it does not establish
//! causality or authorize post-hoc lag selection in confirmatory experiments.

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, PartialEq)]
pub enum LaggedInformationError {
    EmptyInput,
    LengthMismatch {
        observed: usize,
        hidden_state: usize,
    },
    TooFewBins {
        bins: usize,
    },
    EmptyLagSet,
    InvalidLag {
        lag: i64,
        samples: usize,
    },
    DuplicateLag {
        lag: i64,
    },
    NonFiniteObservation {
        index: usize,
        value: f64,
    },
    UnrepresentableObservationRange {
        minimum: f64,
        maximum: f64,
    },
}

impl Display for LaggedInformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "lagged information requires at least one sample"),
            Self::LengthMismatch { observed, hidden_state } => write!(f, "observation/state length mismatch: observed={observed}, hidden_state={hidden_state}"),
            Self::TooFewBins { bins } => write!(f, "lagged histogram mutual information requires at least 2 bins, got {bins}"),
            Self::EmptyLagSet => write!(f, "lagged information requires at least one declared lag"),
            Self::InvalidLag { lag, samples } => write!(f, "lag {lag} leaves no aligned samples for a {samples}-sample series"),
            Self::DuplicateLag { lag } => write!(f, "lag {lag} is duplicated"),
            Self::NonFiniteObservation { index, value } => write!(f, "observation at index {index} is not finite: {value}"),
            Self::UnrepresentableObservationRange { minimum, maximum } => write!(
                f,
                "observation range cannot be represented safely: minimum={minimum}, maximum={maximum}"
            ),
        }
    }
}

impl Error for LaggedInformationError {}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LaggedMutualInformationPoint {
    pub lag: i64,
    pub aligned_samples: usize,
    pub mutual_information_bits: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LaggedMutualInformationScan {
    pub bins: usize,
    pub points: Vec<LaggedMutualInformationPoint>,
}

/// Measure histogram mutual information across a preregistered lag set.
///
/// The equal-width observation bins are frozen from the complete observation
/// before lag alignment. For `lag > 0`, `observed[t]` is paired with
/// `hidden_state[t + lag]`; for `lag < 0`, `observed[t - lag]` is paired with
/// `hidden_state[t]`. No circular wrapping or imputation occurs.
///
/// This reports association only. Confirmatory experiments must freeze the lag
/// set before outcome inspection or separately account for lag selection.
pub fn lagged_histogram_mutual_information_bits(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
    lags: &[i64],
) -> Result<LaggedMutualInformationScan, LaggedInformationError> {
    validate_inputs(observed, hidden_state, bins)?;
    if lags.is_empty() {
        return Err(LaggedInformationError::EmptyLagSet);
    }

    let samples = observed.len();
    let quantized = quantize_full_observation(observed, bins)?;
    let mut seen = BTreeSet::new();
    let mut points = Vec::with_capacity(lags.len());

    for &lag in lags {
        if !seen.insert(lag) {
            return Err(LaggedInformationError::DuplicateLag { lag });
        }
        let magnitude = lag
            .checked_abs()
            .and_then(|value| usize::try_from(value).ok())
            .ok_or(LaggedInformationError::InvalidLag { lag, samples })?;
        if magnitude >= samples {
            return Err(LaggedInformationError::InvalidLag { lag, samples });
        }

        let (observation_slice, state_slice) = if lag >= 0 {
            (
                &quantized[..samples - magnitude],
                &hidden_state[magnitude..],
            )
        } else {
            (
                &quantized[magnitude..],
                &hidden_state[..samples - magnitude],
            )
        };

        points.push(LaggedMutualInformationPoint {
            lag,
            aligned_samples: observation_slice.len(),
            mutual_information_bits: mutual_information_from_quantized(
                observation_slice,
                state_slice,
            ),
        });
    }

    Ok(LaggedMutualInformationScan { bins, points })
}

fn validate_inputs(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
) -> Result<(), LaggedInformationError> {
    if observed.is_empty() || hidden_state.is_empty() {
        return Err(LaggedInformationError::EmptyInput);
    }
    if observed.len() != hidden_state.len() {
        return Err(LaggedInformationError::LengthMismatch {
            observed: observed.len(),
            hidden_state: hidden_state.len(),
        });
    }
    if bins < 2 {
        return Err(LaggedInformationError::TooFewBins { bins });
    }
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(LaggedInformationError::NonFiniteObservation { index, value });
        }
    }
    Ok(())
}

fn quantize_full_observation(
    observed: &[f64],
    bins: usize,
) -> Result<Vec<usize>, LaggedInformationError> {
    let mut minimum = f64::INFINITY;
    let mut maximum = f64::NEG_INFINITY;
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(LaggedInformationError::NonFiniteObservation { index, value });
        }
        minimum = minimum.min(value);
        maximum = maximum.max(value);
    }
    if minimum == maximum {
        return Ok(vec![0; observed.len()]);
    }
    // Normalize before subtraction so a finite span such as [-1e308, 1e308]
    // cannot overflow to infinity. If distinct finite endpoints collapse after
    // normalization, fail closed instead of silently assigning the wrong bin.
    let scale = minimum.abs().max(maximum.abs());
    if !scale.is_finite() || scale == 0.0 {
        return Err(LaggedInformationError::UnrepresentableObservationRange { minimum, maximum });
    }
    let normalized_minimum = minimum / scale;
    let normalized_maximum = maximum / scale;
    let normalized_span = normalized_maximum - normalized_minimum;
    if !normalized_span.is_finite() || normalized_span <= 0.0 {
        return Err(LaggedInformationError::UnrepresentableObservationRange { minimum, maximum });
    }

    Ok(observed
        .iter()
        .map(|&value| {
            let normalized = value / scale;
            let scaled =
                ((normalized - normalized_minimum) / normalized_span * bins as f64).floor();
            if scaled <= 0.0 {
                0
            } else if scaled >= bins as f64 {
                bins - 1
            } else {
                scaled as usize
            }
        })
        .collect())
}

fn mutual_information_from_quantized(observed: &[usize], hidden_state: &[usize]) -> f64 {
    let mut observation_counts = BTreeMap::<usize, usize>::new();
    let mut state_counts = BTreeMap::<usize, usize>::new();
    let mut joint_counts = BTreeMap::<(usize, usize), usize>::new();
    for (&observation, &state) in observed.iter().zip(hidden_state) {
        *observation_counts.entry(observation).or_insert(0) += 1;
        *state_counts.entry(state).or_insert(0) += 1;
        *joint_counts.entry((observation, state)).or_insert(0) += 1;
    }
    let samples = observed.len() as f64;
    joint_counts
        .into_iter()
        .map(|((observation, state), joint_count)| {
            let p_joint = joint_count as f64 / samples;
            let p_observation = observation_counts[&observation] as f64 / samples;
            let p_state = state_counts[&state] as f64 / samples;
            p_joint * (p_joint / (p_observation * p_state)).log2()
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64) {
        assert!((actual - expected).abs() <= 1e-12, "{actual} != {expected}");
    }

    #[test]
    fn declaration_order_and_overlap_are_retained() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let scan =
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[0, 2, -2]).unwrap();
        assert_eq!(
            scan.points
                .iter()
                .map(|point| point.lag)
                .collect::<Vec<_>>(),
            vec![0, 2, -2]
        );
        assert_eq!(
            scan.points
                .iter()
                .map(|point| point.aligned_samples)
                .collect::<Vec<_>>(),
            vec![8, 6, 6]
        );
        assert!(scan
            .points
            .iter()
            .all(|point| (point.mutual_information_bits - 1.0).abs() <= 1e-12));
    }

    #[test]
    fn noncircular_alignment_detects_declared_state_lead() {
        let hidden = [0, 0, 1, 1, 0, 0, 1, 1];
        let observed = [-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0];
        let scan =
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[0, 1]).unwrap();
        assert_close(scan.points[0].mutual_information_bits, 0.0);
        assert!(scan.points[1].mutual_information_bits > 0.5);
        assert_eq!(scan.points[1].aligned_samples, 7);
    }

    #[test]
    fn extreme_finite_span_does_not_overflow_quantization() {
        let hidden = [0, 1];
        let observed = [-1.0e308, 1.0e308];
        let scan = lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[0]).unwrap();
        assert_close(scan.points[0].mutual_information_bits, 1.0);
    }

    #[test]
    fn constant_observation_has_zero_information() {
        let hidden = [0, 1, 0, 1, 0];
        let observed = [3.0; 5];
        let scan =
            lagged_histogram_mutual_information_bits(&observed, &hidden, 4, &[-2, 0, 2]).unwrap();
        assert!(scan
            .points
            .iter()
            .all(|point| point.mutual_information_bits == 0.0));
    }

    #[test]
    fn invalid_or_duplicate_lags_fail_closed() {
        let observed = [0.0, 1.0, 2.0, 3.0];
        let hidden = [0, 0, 1, 1];
        assert_eq!(
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[]),
            Err(LaggedInformationError::EmptyLagSet)
        );
        assert_eq!(
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[1, 1]),
            Err(LaggedInformationError::DuplicateLag { lag: 1 })
        );
        assert_eq!(
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[4]),
            Err(LaggedInformationError::InvalidLag { lag: 4, samples: 4 })
        );
        assert_eq!(
            lagged_histogram_mutual_information_bits(&observed, &hidden, 2, &[-4]),
            Err(LaggedInformationError::InvalidLag {
                lag: -4,
                samples: 4
            })
        );
    }

    #[test]
    fn malformed_inputs_fail_closed() {
        assert_eq!(
            lagged_histogram_mutual_information_bits(&[], &[], 2, &[0]),
            Err(LaggedInformationError::EmptyInput)
        );
        assert_eq!(
            lagged_histogram_mutual_information_bits(&[0.0], &[0, 1], 2, &[0]),
            Err(LaggedInformationError::LengthMismatch {
                observed: 1,
                hidden_state: 2
            })
        );
        assert_eq!(
            lagged_histogram_mutual_information_bits(&[0.0], &[0], 1, &[0]),
            Err(LaggedInformationError::TooFewBins { bins: 1 })
        );
    }
}
