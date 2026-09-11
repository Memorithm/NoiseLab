//! Diagnostics for testing whether an apparently noisy component carries
//! information about a declared hidden or system state.
//!
//! These routines intentionally implement a small, dependency-free histogram
//! estimator. They are suitable for controlled and preregistered experiments,
//! not for claiming exact continuous mutual information.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Errors returned by information-preservation diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub enum InformationError {
    EmptyInput,
    LengthMismatch {
        observed: usize,
        hidden_state: usize,
    },
    TooFewBins {
        bins: usize,
    },
    NonFiniteObservation {
        index: usize,
        value: f64,
    },
}

impl Display for InformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyInput => write!(f, "information diagnostics require at least one sample"),
            Self::LengthMismatch {
                observed,
                hidden_state,
            } => write!(
                f,
                "observation/state length mismatch: observed={observed}, hidden_state={hidden_state}"
            ),
            Self::TooFewBins { bins } => write!(
                f,
                "histogram mutual information requires at least 2 bins, got {bins}"
            ),
            Self::NonFiniteObservation { index, value } => write!(
                f,
                "observation at index {index} is not finite: {value}"
            ),
        }
    }
}

impl Error for InformationError {}

/// Information carried by an apparently noisy component before and after a
/// declared transformation such as denoising, smoothing or thresholding.
///
/// Values are empirical plug-in estimates in bits. `retained_fraction_of_raw`
/// is `None` when the raw estimate is zero because the ratio is undefined.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoiseInformationAudit {
    pub hidden_state_entropy_bits: f64,
    pub raw_mutual_information_bits: f64,
    pub transformed_mutual_information_bits: f64,
    pub raw_fraction_of_hidden_entropy: f64,
    pub transformed_fraction_of_hidden_entropy: f64,
    /// Signed change: transformed minus raw mutual information.
    pub information_delta_bits: f64,
    pub retained_fraction_of_raw: Option<f64>,
}

/// Empirical entropy of discrete labels, in bits.
pub fn discrete_entropy_bits(labels: &[usize]) -> Result<f64, InformationError> {
    if labels.is_empty() {
        return Err(InformationError::EmptyInput);
    }

    let mut counts = BTreeMap::<usize, usize>::new();
    for &label in labels {
        *counts.entry(label).or_insert(0) += 1;
    }

    let n = labels.len() as f64;
    Ok(counts
        .values()
        .map(|&count| {
            let p = count as f64 / n;
            -p * p.log2()
        })
        .sum())
}

/// Estimate mutual information between an observation and a discrete
/// hidden/system state using equal-width observation bins.
///
/// This is a histogram plug-in estimator. The result depends on `bins`, sample
/// count and observation range. Experiments should preregister binning and
/// report sensitivity analyses rather than interpreting this as exact
/// continuous mutual information.
pub fn histogram_mutual_information_bits(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
) -> Result<f64, InformationError> {
    validate_inputs(observed, hidden_state, bins)?;
    let quantized = quantize_equal_width(observed, bins)?;

    let mut observation_counts = vec![0usize; bins];
    let mut state_counts = BTreeMap::<usize, usize>::new();
    let mut joint_counts = BTreeMap::<(usize, usize), usize>::new();

    for (&bin, &state) in quantized.iter().zip(hidden_state) {
        observation_counts[bin] += 1;
        *state_counts.entry(state).or_insert(0) += 1;
        *joint_counts.entry((bin, state)).or_insert(0) += 1;
    }

    let n = observed.len() as f64;
    let mut mutual_information = 0.0;

    for (&(bin, state), &joint_count) in &joint_counts {
        let p_joint = joint_count as f64 / n;
        let p_observed = observation_counts[bin] as f64 / n;
        let p_state = state_counts[&state] as f64 / n;
        mutual_information += p_joint * (p_joint / (p_observed * p_state)).log2();
    }

    Ok(mutual_information.max(0.0))
}

/// Compare how much declared hidden-state information is present in an
/// apparently noisy component before and after a transformation.
///
/// Raw and transformed series are quantized independently with the same number
/// of equal-width bins. This makes simple affine rescaling benign while still
/// exposing transformations that collapse state-dependent structure.
pub fn audit_noise_information(
    raw_component: &[f64],
    transformed_component: &[f64],
    hidden_state: &[usize],
    bins: usize,
) -> Result<NoiseInformationAudit, InformationError> {
    validate_inputs(raw_component, hidden_state, bins)?;
    validate_inputs(transformed_component, hidden_state, bins)?;

    let hidden_state_entropy_bits = discrete_entropy_bits(hidden_state)?;
    let raw_mutual_information_bits =
        histogram_mutual_information_bits(raw_component, hidden_state, bins)?;
    let transformed_mutual_information_bits =
        histogram_mutual_information_bits(transformed_component, hidden_state, bins)?;

    let raw_fraction_of_hidden_entropy = if hidden_state_entropy_bits > 0.0 {
        raw_mutual_information_bits / hidden_state_entropy_bits
    } else {
        0.0
    };
    let transformed_fraction_of_hidden_entropy = if hidden_state_entropy_bits > 0.0 {
        transformed_mutual_information_bits / hidden_state_entropy_bits
    } else {
        0.0
    };
    let information_delta_bits = transformed_mutual_information_bits - raw_mutual_information_bits;
    let retained_fraction_of_raw = if raw_mutual_information_bits > 0.0 {
        Some(transformed_mutual_information_bits / raw_mutual_information_bits)
    } else {
        None
    };

    Ok(NoiseInformationAudit {
        hidden_state_entropy_bits,
        raw_mutual_information_bits,
        transformed_mutual_information_bits,
        raw_fraction_of_hidden_entropy,
        transformed_fraction_of_hidden_entropy,
        information_delta_bits,
        retained_fraction_of_raw,
    })
}

fn validate_inputs(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
) -> Result<(), InformationError> {
    if observed.is_empty() || hidden_state.is_empty() {
        return Err(InformationError::EmptyInput);
    }
    if observed.len() != hidden_state.len() {
        return Err(InformationError::LengthMismatch {
            observed: observed.len(),
            hidden_state: hidden_state.len(),
        });
    }
    if bins < 2 {
        return Err(InformationError::TooFewBins { bins });
    }
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(InformationError::NonFiniteObservation { index, value });
        }
    }
    Ok(())
}

fn quantize_equal_width(observed: &[f64], bins: usize) -> Result<Vec<usize>, InformationError> {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(InformationError::NonFiniteObservation { index, value });
        }
        min = min.min(value);
        max = max.max(value);
    }

    if min == max {
        return Ok(vec![0; observed.len()]);
    }

    let width = (max - min) / bins as f64;
    Ok(observed
        .iter()
        .map(|&value| {
            let raw_bin = ((value - min) / width).floor() as usize;
            raw_bin.min(bins - 1)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual}, expected={expected}, tolerance={tolerance}"
        );
    }

    #[test]
    fn balanced_binary_state_has_one_bit_entropy() {
        let labels = [0, 1, 0, 1, 0, 1, 0, 1];
        assert_close(discrete_entropy_bits(&labels).unwrap(), 1.0, 1e-12);
    }

    #[test]
    fn known_state_coded_component_carries_one_bit() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let mi = histogram_mutual_information_bits(&observed, &hidden, 2).unwrap();
        assert_close(mi, 1.0, 1e-12);
    }

    #[test]
    fn balanced_independent_control_has_zero_mutual_information() {
        let hidden = [0, 0, 1, 1, 0, 0, 1, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let mi = histogram_mutual_information_bits(&observed, &hidden, 2).unwrap();
        assert_close(mi, 0.0, 1e-12);
    }

    #[test]
    fn destructive_denoising_can_erase_hidden_state_information() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let raw = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let denoised = [0.0; 8];

        let audit = audit_noise_information(&raw, &denoised, &hidden, 2).unwrap();
        assert_close(audit.hidden_state_entropy_bits, 1.0, 1e-12);
        assert_close(audit.raw_mutual_information_bits, 1.0, 1e-12);
        assert_close(audit.transformed_mutual_information_bits, 0.0, 1e-12);
        assert_close(audit.information_delta_bits, -1.0, 1e-12);
        assert_eq!(audit.retained_fraction_of_raw, Some(0.0));
    }

    #[test]
    fn identity_transform_retains_information() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let raw = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];

        let audit = audit_noise_information(&raw, &raw, &hidden, 2).unwrap();
        assert_close(audit.raw_mutual_information_bits, 1.0, 1e-12);
        assert_close(audit.transformed_mutual_information_bits, 1.0, 1e-12);
        assert_close(audit.information_delta_bits, 0.0, 1e-12);
        assert_eq!(audit.retained_fraction_of_raw, Some(1.0));
    }

    #[test]
    fn malformed_inputs_fail_closed() {
        assert_eq!(
            histogram_mutual_information_bits(&[], &[], 2),
            Err(InformationError::EmptyInput)
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0], &[0, 1], 2),
            Err(InformationError::LengthMismatch {
                observed: 1,
                hidden_state: 2,
            })
        );
        assert_eq!(
            histogram_mutual_information_bits(&[0.0], &[0], 1),
            Err(InformationError::TooFewBins { bins: 1 })
        );
        assert!(matches!(
            histogram_mutual_information_bits(&[f64::NAN], &[0], 2),
            Err(InformationError::NonFiniteObservation { index: 0, .. })
        ));
    }
}
