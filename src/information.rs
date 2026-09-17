//! Diagnostics for testing whether an apparently noisy component carries
//! information about a declared hidden or system state.
//!
//! These routines intentionally implement a small histogram estimator plus a
//! deterministic permutation null. They are suitable for controlled and
//! preregistered experiments, not for claiming exact continuous mutual
//! information or causal mechanism.

use scirust_sim::SplitMix64;
use std::collections::BTreeMap;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Hard safety ceiling for a single permutation-null request.
pub const MAX_INFORMATION_PERMUTATIONS: usize = 100_000;

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
    TooFewPermutations {
        permutations: usize,
    },
    TooManyPermutations {
        requested: usize,
        maximum: usize,
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
            Self::TooFewPermutations { permutations } => write!(
                f,
                "permutation null requires at least 1 permutation, got {permutations}"
            ),
            Self::TooManyPermutations { requested, maximum } => write!(
                f,
                "requested {requested} permutations exceeds safety maximum {maximum}"
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

/// Deterministic finite-sample null for histogram mutual information.
///
/// The null preserves the observed values and the exact hidden-state label
/// multiset, then independently permutes the labels for every surrogate. The
/// one-sided p-value uses the standard finite permutation correction
/// `(1 + exceedances) / (1 + permutations)` and therefore never reports zero.
/// This is an association-bias diagnostic, not evidence of causality.
#[derive(Debug, Clone, PartialEq)]
pub struct MutualInformationPermutationNull {
    pub observed_mutual_information_bits: f64,
    pub surrogate_mean_mutual_information_bits: f64,
    pub exceedances_at_or_above_observed: usize,
    pub permutation_p_value: f64,
    pub seed: u64,
    pub surrogate_mutual_information_bits: Vec<f64>,
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
    Ok(mutual_information_from_quantized(
        &quantized,
        hidden_state,
        bins,
    ))
}

/// Compare the observed histogram mutual information with a label-permutation
/// null that preserves the observation vector and exact hidden-state marginal.
///
/// The same quantized observation vector is reused for every surrogate, so the
/// permutation test changes only the declared state association. Each
/// surrogate starts from the original label order and is shuffled using a
/// seeded Fisher-Yates permutation with rejection-sampled indices, avoiding
/// modulo bias. The returned surrogate scores are retained explicitly so
/// downstream reports can preserve the actual finite null distribution.
///
/// Experiments must preregister `bins`, `permutations` and `seed`. A small
/// p-value rejects only this permutation null under the declared estimator and
/// exchangeability assumption; it does not establish mechanism or causality.
pub fn permutation_mutual_information_null(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
    permutations: usize,
    seed: u64,
) -> Result<MutualInformationPermutationNull, InformationError> {
    validate_inputs(observed, hidden_state, bins)?;
    if permutations == 0 {
        return Err(InformationError::TooFewPermutations { permutations });
    }
    if permutations > MAX_INFORMATION_PERMUTATIONS {
        return Err(InformationError::TooManyPermutations {
            requested: permutations,
            maximum: MAX_INFORMATION_PERMUTATIONS,
        });
    }

    let quantized = quantize_equal_width(observed, bins)?;
    let observed_mutual_information_bits =
        mutual_information_from_quantized(&quantized, hidden_state, bins);
    let mut rng = SplitMix64::new(seed);
    let mut shuffled = hidden_state.to_vec();
    let mut surrogate_mutual_information_bits = Vec::new();
    surrogate_mutual_information_bits
        .try_reserve_exact(permutations)
        .map_err(|_| InformationError::TooManyPermutations {
            requested: permutations,
            maximum: MAX_INFORMATION_PERMUTATIONS,
        })?;

    let mut exceedances_at_or_above_observed = 0usize;
    let mut surrogate_sum = 0.0;
    for _ in 0..permutations {
        shuffled.copy_from_slice(hidden_state);
        fisher_yates_shuffle(&mut shuffled, &mut rng);
        let surrogate = mutual_information_from_quantized(&quantized, &shuffled, bins);
        if surrogate >= observed_mutual_information_bits {
            exceedances_at_or_above_observed += 1;
        }
        surrogate_sum += surrogate;
        surrogate_mutual_information_bits.push(surrogate);
    }

    let surrogate_mean_mutual_information_bits = surrogate_sum / permutations as f64;
    let permutation_p_value =
        (exceedances_at_or_above_observed as f64 + 1.0) / (permutations as f64 + 1.0);

    Ok(MutualInformationPermutationNull {
        observed_mutual_information_bits,
        surrogate_mean_mutual_information_bits,
        exceedances_at_or_above_observed,
        permutation_p_value,
        seed,
        surrogate_mutual_information_bits,
    })
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

fn mutual_information_from_quantized(
    quantized: &[usize],
    hidden_state: &[usize],
    bins: usize,
) -> f64 {
    let mut observation_counts = vec![0usize; bins];
    let mut state_counts = BTreeMap::<usize, usize>::new();
    let mut joint_counts = BTreeMap::<(usize, usize), usize>::new();

    for (&bin, &state) in quantized.iter().zip(hidden_state) {
        observation_counts[bin] += 1;
        *state_counts.entry(state).or_insert(0) += 1;
        *joint_counts.entry((bin, state)).or_insert(0) += 1;
    }

    // Canonicalize contribution order by sufficient counts, not state labels.
    // A pure relabeling therefore produces bit-identical MI and cannot turn a
    // mathematical tie into a sub-observed surrogate through addition order.
    let mut terms = Vec::with_capacity(joint_counts.len());
    for (&(bin, state), &joint_count) in &joint_counts {
        terms.push((joint_count, observation_counts[bin], state_counts[&state]));
    }
    terms.sort_unstable();

    let n = quantized.len() as f64;
    let mut mutual_information = 0.0;
    for (joint_count, observation_count, state_count) in terms {
        let p_joint = joint_count as f64 / n;
        let p_observed = observation_count as f64 / n;
        let p_state = state_count as f64 / n;
        mutual_information += p_joint * (p_joint / (p_observed * p_state)).log2();
    }
    mutual_information.max(0.0)
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
    fn mutual_information_is_bit_invariant_to_state_relabeling() {
        let quantized = [0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2];
        let states = [0, 0, 1, 2, 0, 1, 1, 1, 2, 2, 0, 0, 0, 1, 2, 2];
        let relabeled = states.map(|state| match state {
            0 => 17,
            1 => 3,
            2 => 11,
            _ => unreachable!(),
        });

        let original = mutual_information_from_quantized(&quantized, &states, 3);
        let renamed = mutual_information_from_quantized(&quantized, &relabeled, 3);
        assert_eq!(original.to_bits(), renamed.to_bits());
    }

    #[test]
    fn permutation_null_is_seed_reproducible_and_retains_scores() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0];
        let a = permutation_mutual_information_null(&observed, &hidden, 2, 64, 42).unwrap();
        let b = permutation_mutual_information_null(&observed, &hidden, 2, 64, 42).unwrap();
        assert_eq!(a, b);
        assert_eq!(a.surrogate_mutual_information_bits.len(), 64);
        assert_close(a.observed_mutual_information_bits, 1.0, 1e-12);
        assert!(a.permutation_p_value > 0.0);
        assert!(a.permutation_p_value <= 1.0);
    }

    #[test]
    fn permutation_null_preserves_label_marginal_but_breaks_pairing() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [
            -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
        ];
        let null = permutation_mutual_information_null(&observed, &hidden, 2, 256, 7).unwrap();
        assert_close(null.observed_mutual_information_bits, 1.0, 1e-12);
        assert!(null.surrogate_mean_mutual_information_bits < 1.0);
        assert!(null.exceedances_at_or_above_observed < 256);
    }

    #[test]
    fn permutation_null_rejects_unbounded_or_empty_requests() {
        let observed = [0.0, 1.0];
        let hidden = [0, 1];
        assert_eq!(
            permutation_mutual_information_null(&observed, &hidden, 2, 0, 1),
            Err(InformationError::TooFewPermutations { permutations: 0 })
        );
        assert_eq!(
            permutation_mutual_information_null(
                &observed,
                &hidden,
                2,
                MAX_INFORMATION_PERMUTATIONS + 1,
                1,
            ),
            Err(InformationError::TooManyPermutations {
                requested: MAX_INFORMATION_PERMUTATIONS + 1,
                maximum: MAX_INFORMATION_PERMUTATIONS,
            })
        );
    }

    #[test]
    fn constant_observation_has_zero_observed_and_null_information() {
        let hidden = [0, 1, 0, 1, 0, 1, 0, 1];
        let observed = [3.0; 8];
        let null = permutation_mutual_information_null(&observed, &hidden, 4, 32, 9).unwrap();
        assert_close(null.observed_mutual_information_bits, 0.0, 1e-12);
        assert!(null
            .surrogate_mutual_information_bits
            .iter()
            .all(|&score| score == 0.0));
        assert_eq!(null.exceedances_at_or_above_observed, 32);
        assert_close(null.permutation_p_value, 1.0, 1e-12);
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
