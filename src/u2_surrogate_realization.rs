//! Canonical materialization of one preregistered Stage U2 surrogate pair.
//!
//! This module exists to give capture/replay code one bounded, deterministic way
//! to reproduce the exact null transformation attached to a [`U2SurrogateJob`].
//! It does not compute p-values, inspect outcomes, alter preregistered seeds, or
//! promote a surrogate realization into scientific evidence.

use crate::spectral_null::spectral_phase_null;
use crate::u2_plan::{U2NullFamily, U2SurrogateJob};
use crate::universality_u2::SPECTRAL_RIGHT_SEED_TAG;
use scirust_signal::surrogate::SurrogateError;
use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Exact left/right arrays materialized for one preregistered U2 surrogate job.
#[derive(Debug, Clone, PartialEq)]
pub struct U2SurrogateRealization {
    /// Outcome-blind job identity that fixes null family, repetition and seed.
    pub job: U2SurrogateJob,
    /// Realized left-hand surrogate array.
    pub left: Vec<f64>,
    /// Realized right-hand surrogate array.
    pub right: Vec<f64>,
}

/// Fail-closed errors while materializing a U2 surrogate pair.
#[derive(Debug)]
pub enum U2SurrogateRealizationError {
    /// Left/right source arrays have different lengths.
    LengthMismatch { left: usize, right: usize },
    /// An input array contains a non-finite value.
    NonFiniteInput { side: &'static str, index: usize },
    /// Allocation for a shuffled-marginal realization failed.
    AllocationFailed,
    /// SciRust rejected a phase-randomized spectral realization.
    Spectral(SurrogateError),
}

impl Display for U2SurrogateRealizationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LengthMismatch { left, right } => {
                write!(f, "U2 surrogate inputs have different lengths: left={left}, right={right}")
            }
            Self::NonFiniteInput { side, index } => {
                write!(f, "U2 surrogate {side} input contains a non-finite value at index {index}")
            }
            Self::AllocationFailed => f.write_str("unable to allocate U2 surrogate realization"),
            Self::Spectral(error) => write!(f, "U2 spectral-null realization failed: {error}"),
        }
    }
}

impl Error for U2SurrogateRealizationError {}

impl From<SurrogateError> for U2SurrogateRealizationError {
    fn from(value: SurrogateError) -> Self {
        Self::Spectral(value)
    }
}

/// Materialize one exact U2 surrogate pair from a preregistered job.
///
/// Shuffled-marginal jobs use one `SplitMix64` stream, shuffling the left array
/// first and the right array second, exactly as the existing U2 analysis path.
/// Phase-randomized jobs use the frozen left seed and the independently tagged
/// right seed. The function is deterministic for fixed inputs and job identity.
///
/// This is an evidence/replay primitive, not a scientific decision. Persisting
/// its output does not itself establish that a panel completed or that a result
/// is positive.
///
/// # Errors
///
/// Returns an error for mismatched lengths, non-finite inputs, allocation
/// failure, or a rejected spectral surrogate operation.
pub fn materialize_u2_surrogate_pair(
    job: U2SurrogateJob,
    series_a: &[f64],
    series_b: &[f64],
) -> Result<U2SurrogateRealization, U2SurrogateRealizationError> {
    if series_a.len() != series_b.len() {
        return Err(U2SurrogateRealizationError::LengthMismatch {
            left: series_a.len(),
            right: series_b.len(),
        });
    }
    if let Some(index) = series_a.iter().position(|value| !value.is_finite()) {
        return Err(U2SurrogateRealizationError::NonFiniteInput {
            side: "left",
            index,
        });
    }
    if let Some(index) = series_b.iter().position(|value| !value.is_finite()) {
        return Err(U2SurrogateRealizationError::NonFiniteInput {
            side: "right",
            index,
        });
    }

    let (left, right) = match job.null_family {
        U2NullFamily::ShuffledMarginal => {
            let mut left = copy_series(series_a)?;
            let mut right = copy_series(series_b)?;
            let mut rng = SplitMix64::new(job.seed);
            fisher_yates_shuffle(&mut left, &mut rng);
            fisher_yates_shuffle(&mut right, &mut rng);
            (left, right)
        }
        U2NullFamily::PhaseRandomizedSpectrum => (
            spectral_phase_null(series_a, job.seed)?,
            spectral_phase_null(series_b, job.seed ^ SPECTRAL_RIGHT_SEED_TAG)?,
        ),
    };

    Ok(U2SurrogateRealization { job, left, right })
}

fn copy_series(values: &[f64]) -> Result<Vec<f64>, U2SurrogateRealizationError> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| U2SurrogateRealizationError::AllocationFailed)?;
    copy.extend_from_slice(values);
    Ok(copy)
}

fn fisher_yates_shuffle(values: &mut [f64], rng: &mut SplitMix64) {
    for upper in (1..values.len()).rev() {
        let index = unbiased_index(rng, upper + 1);
        values.swap(upper, index);
    }
}

fn unbiased_index(rng: &mut SplitMix64, upper_exclusive: usize) -> usize {
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
    use crate::u2_manifest::materialize_u2_manifest;
    use crate::u2_plan::{U2ExecutionPlan, U2_FROZEN_PAIRS};
    use crate::universality::observed_multiscale_comparison;
    use crate::universality_u2::{analyze_u2_pair, StageU2PairRequest};

    fn synthetic_pair() -> (Vec<f64>, Vec<f64>) {
        let left = (0..256).map(|i| (0.11 * i as f64).sin()).collect();
        let right = (0..256)
            .map(|i| (0.07 * i as f64).cos() + 0.05 * i as f64)
            .collect();
        (left, right)
    }

    #[test]
    fn realization_is_bit_reproducible_for_each_null() {
        let (left, right) = synthetic_pair();
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: 19,
            seed_root: 0xABCD,
        };
        for null_family in [
            U2NullFamily::ShuffledMarginal,
            U2NullFamily::PhaseRandomizedSpectrum,
        ] {
            let job = plan.surrogate_job(0, null_family, 7).unwrap();
            let first = materialize_u2_surrogate_pair(job, &left, &right).unwrap();
            let second = materialize_u2_surrogate_pair(job, &left, &right).unwrap();
            assert_eq!(first, second);
            assert_eq!(first.job, job);
        }
    }

    #[test]
    fn canonical_realizations_match_existing_analysis_scores_exactly() {
        let (left, right) = synthetic_pair();
        let scales = [1, 2, 4];
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: 19,
            seed_root: 0xABCD,
        };
        let jobs = materialize_u2_manifest(plan);
        let analysis = analyze_u2_pair(StageU2PairRequest {
            pair_index: 0,
            pair: U2_FROZEN_PAIRS[0],
            series_a: &left,
            series_b: &right,
            scales: &scales,
            alpha: 0.05,
            jobs: &jobs,
            retain_protocol_failures: false,
        })
        .unwrap();

        for null_family in [
            U2NullFamily::ShuffledMarginal,
            U2NullFamily::PhaseRandomizedSpectrum,
        ] {
            let score = analysis
                .surrogate_scores
                .iter()
                .find(|score| score.job.null_family == null_family && score.job.repetition == 0)
                .unwrap();
            let realization =
                materialize_u2_surrogate_pair(score.job, &left, &right).unwrap();
            let reproduced_score = observed_multiscale_comparison(
                &realization.left,
                &realization.right,
                &scales,
            )
            .unwrap()
            .convergence_score;
            assert_eq!(reproduced_score.to_bits(), score.convergence_score.to_bits());
        }
    }

    #[test]
    fn invalid_inputs_fail_closed_before_realization() {
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: 19,
            seed_root: 7,
        };
        let job = plan
            .surrogate_job(0, U2NullFamily::ShuffledMarginal, 0)
            .unwrap();
        assert!(matches!(
            materialize_u2_surrogate_pair(job, &[1.0, 2.0], &[1.0]),
            Err(U2SurrogateRealizationError::LengthMismatch { .. })
        ));
        assert!(matches!(
            materialize_u2_surrogate_pair(job, &[1.0, f64::NAN], &[1.0, 2.0]),
            Err(U2SurrogateRealizationError::NonFiniteInput {
                side: "left",
                index: 1,
            })
        ));
    }
}
