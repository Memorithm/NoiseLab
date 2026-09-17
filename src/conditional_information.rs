//! Histogram conditional-mutual-information calibration primitive.
//!
//! This is a finite-sample plug-in estimator for a real-valued observation X,
//! a discrete declared state Y, and a discrete conditioning variable Z.  It is
//! intended for controlled diagnostics of whether an apparent X-Y association
//! remains after conditioning on a declared measured state.  It is not a causal
//! estimator and does not identify an unobserved mechanism.

use crate::information::{quantize_equal_width, InformationError};
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Errors returned by the conditional-information calibration primitive.
#[derive(Debug, Clone, PartialEq)]
pub enum ConditionalInformationError {
    /// Canonical observation/hidden-state validation or quantization failed.
    Information(InformationError),
    /// The conditioning variable must have exactly one label per observation.
    ConditionLengthMismatch { observed: usize, condition: usize },
}

impl Display for ConditionalInformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Information(error) => write!(f, "conditional-information input error: {error}"),
            Self::ConditionLengthMismatch {
                observed,
                condition,
            } => write!(
                f,
                "observation/condition length mismatch: observed={observed}, condition={condition}"
            ),
        }
    }
}

impl Error for ConditionalInformationError {}

impl From<InformationError> for ConditionalInformationError {
    fn from(value: InformationError) -> Self {
        Self::Information(value)
    }
}

/// Estimate `I(X;Y | Z)` in bits from one frozen equal-width quantization of X.
///
/// The empirical plug-in estimate is
///
/// `sum p(x,y,z) log2( p(x,y,z) p(z) / (p(x,z) p(y,z)) )`.
///
/// `bins` controls only the real-valued observation X. The discrete hidden and
/// conditioning labels are used exactly as supplied. Label numeric values carry
/// no metric meaning. Experiments must freeze the binning and conditioning
/// variable before outcome inspection and report finite-sample sensitivity.
///
/// A positive estimate means only that residual empirical association remains
/// under this declared conditioning variable and histogram estimator. It does
/// not establish causal direction, a physical mechanism, sufficiency of Z, or
/// that the residual is useful information rather than estimator bias.
pub fn conditional_histogram_mutual_information_bits(
    observed: &[f64],
    hidden_state: &[usize],
    condition: &[usize],
    bins: usize,
) -> Result<f64, ConditionalInformationError> {
    if observed.len() != condition.len() {
        return Err(ConditionalInformationError::ConditionLengthMismatch {
            observed: observed.len(),
            condition: condition.len(),
        });
    }

    // Reuse the canonical observation quantizer so the unconditional and
    // conditional diagnostics have identical X binning semantics.
    if observed.is_empty() || hidden_state.is_empty() {
        return Err(InformationError::EmptyInput.into());
    }
    if observed.len() != hidden_state.len() {
        return Err(InformationError::LengthMismatch {
            observed: observed.len(),
            hidden_state: hidden_state.len(),
        }
        .into());
    }
    if bins < 2 {
        return Err(InformationError::TooFewBins { bins }.into());
    }
    let quantized = quantize_equal_width(observed, bins)?;

    let mut z_counts = BTreeMap::<usize, usize>::new();
    let mut xz_counts = BTreeMap::<(usize, usize), usize>::new();
    let mut yz_counts = BTreeMap::<(usize, usize), usize>::new();
    let mut xyz_counts = BTreeMap::<(usize, usize, usize), usize>::new();

    for ((&x, &y), &z) in quantized.iter().zip(hidden_state).zip(condition) {
        *z_counts.entry(z).or_insert(0) += 1;
        *xz_counts.entry((x, z)).or_insert(0) += 1;
        *yz_counts.entry((y, z)).or_insert(0) += 1;
        *xyz_counts.entry((x, y, z)).or_insert(0) += 1;
    }

    // Canonicalize summation order by sufficient counts rather than arbitrary
    // numeric labels. Pure relabeling of Y/Z therefore cannot perturb a tie via
    // BTreeMap label order and floating-point addition order.
    let mut terms = Vec::with_capacity(xyz_counts.len());
    for (&(x, y, z), &xyz) in &xyz_counts {
        terms.push((xyz, xz_counts[&(x, z)], yz_counts[&(y, z)], z_counts[&z]));
    }
    terms.sort_unstable();

    let n = observed.len() as f64;
    let mut conditional_information = 0.0;
    for (xyz, xz, yz, z) in terms {
        let p_xyz = xyz as f64 / n;
        // Count-form algebra is equivalent to the probability expression and
        // avoids forming four tiny probabilities independently.
        let ratio = (xyz as f64 * z as f64) / (xz as f64 * yz as f64);
        conditional_information += p_xyz * ratio.log2();
    }

    Ok(conditional_information.max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::histogram_mutual_information_bits;

    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 1e-12,
            "actual={actual}, expected={expected}"
        );
    }

    #[test]
    fn conditioning_on_the_state_itself_removes_state_coded_association() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        assert_close(
            histogram_mutual_information_bits(&observed, &hidden, 2).unwrap(),
            1.0,
        );
        assert_close(
            conditional_histogram_mutual_information_bits(&observed, &hidden, &hidden, 2).unwrap(),
            0.0,
        );
    }

    #[test]
    fn irrelevant_condition_preserves_one_bit_state_coded_association() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let condition = [0, 0, 1, 1, 0, 0, 1, 1];
        assert_close(
            conditional_histogram_mutual_information_bits(&observed, &hidden, &condition, 2)
                .unwrap(),
            1.0,
        );
    }

    #[test]
    fn xor_exposes_conditional_information_hidden_from_unconditional_mi() {
        // X and Y are marginally independent, but Y = X XOR Z. Knowing Z
        // exposes one full bit of association between X and Y.
        let x = [0usize, 0, 1, 1, 0, 0, 1, 1];
        let z = [0usize, 1, 0, 1, 0, 1, 0, 1];
        let y: Vec<usize> = x.iter().zip(&z).map(|(&x, &z)| x ^ z).collect();
        let observed: Vec<f64> = x
            .iter()
            .map(|&value| if value == 0 { -1.0 } else { 1.0 })
            .collect();

        assert_close(
            histogram_mutual_information_bits(&observed, &y, 2).unwrap(),
            0.0,
        );
        assert_close(
            conditional_histogram_mutual_information_bits(&observed, &y, &z, 2).unwrap(),
            1.0,
        );
    }

    #[test]
    fn pure_label_rewrites_are_bit_identical() {
        let observed = [-2.0, -1.0, 1.0, 2.0, -2.0, -1.0, 1.0, 2.0];
        let hidden = [0, 0, 1, 1, 1, 1, 0, 0];
        let condition = [0, 1, 0, 1, 0, 1, 0, 1];
        let hidden_relabel = hidden.map(|v| if v == 0 { 17 } else { 3 });
        let condition_relabel = condition.map(|v| if v == 0 { 42 } else { 9 });
        let original =
            conditional_histogram_mutual_information_bits(&observed, &hidden, &condition, 2)
                .unwrap();
        let relabeled = conditional_histogram_mutual_information_bits(
            &observed,
            &hidden_relabel,
            &condition_relabel,
            2,
        )
        .unwrap();
        assert_eq!(original.to_bits(), relabeled.to_bits());
    }

    #[test]
    fn constant_observation_has_zero_conditional_information() {
        let observed = [3.0; 8];
        let hidden = [0, 1, 0, 1, 1, 0, 1, 0];
        let condition = [0, 0, 1, 1, 0, 0, 1, 1];
        assert_close(
            conditional_histogram_mutual_information_bits(&observed, &hidden, &condition, 4)
                .unwrap(),
            0.0,
        );
    }

    #[test]
    fn condition_length_mismatch_fails_closed() {
        let observed = [0.0, 1.0, 2.0];
        let hidden = [0, 1, 0];
        let condition = [0, 1];
        assert_eq!(
            conditional_histogram_mutual_information_bits(&observed, &hidden, &condition, 2),
            Err(ConditionalInformationError::ConditionLengthMismatch {
                observed: 3,
                condition: 2,
            })
        );
    }

    #[test]
    fn canonical_input_errors_are_preserved() {
        let observed = [0.0, f64::NAN];
        let hidden = [0, 1];
        let condition = [0, 1];
        assert!(matches!(
            conditional_histogram_mutual_information_bits(&observed, &hidden, &condition, 2),
            Err(ConditionalInformationError::Information(
                InformationError::NonFiniteObservation { index: 1, .. }
            ))
        ));
    }
}
