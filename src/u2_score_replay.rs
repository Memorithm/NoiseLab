//! Deterministic replay of Stage U2 p-values and decisions from retained scores.
//!
//! This module deliberately does not regenerate surrogate arrays. It accepts an
//! outcome-blind [`U2ExecutionPlan`], one frozen pair index, the observed
//! convergence score and the retained per-job [`U2SurrogateScore`] sequence. It
//! then revalidates the exact preregistered job ordering before recomputing the
//! one-sided empirical p-values and applying the frozen U2 decision rule.
//!
//! Successful replay is reproducibility evidence only. It does not establish a
//! scientific result, validate the upstream score-generation mechanism, or
//! authorize use of a protected holdout.

use crate::u2_decision::{classify_stage_u2, StageU2Decision, StageU2DecisionError};
use crate::u2_plan::{U2ExecutionPlan, U2NullFamily};
use crate::universality_u2::U2SurrogateScore;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Exact statistics reconstructed from an already retained U2 score sequence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReplayedU2Statistics {
    /// One-sided empirical p-value for shuffled empirical-marginal controls.
    pub p_shuffle: f64,
    /// One-sided empirical p-value for phase-randomized spectral controls.
    pub p_phase: f64,
    /// Frozen Stage U2 classification produced from the replayed p-values.
    pub decision: StageU2Decision,
    /// Number of shuffled controls at least as extreme as the observed score.
    pub shuffled_at_least_as_extreme: usize,
    /// Number of phase-randomized controls at least as extreme as the observed score.
    pub phase_at_least_as_extreme: usize,
}

/// Fail-closed errors for score-only Stage U2 replay.
#[derive(Debug, Clone, PartialEq)]
pub enum U2ScoreReplayError {
    /// The requested pair is not present in the frozen execution plan.
    PairIndexOutOfRange { pair_index: usize },
    /// A plan with zero surrogate repetitions cannot produce an empirical p-value.
    EmptySurrogateSet,
    /// Computing the expected dual-null score count overflowed `usize`.
    ScoreCountOverflow,
    /// The retained score count does not match the exact dual-null plan.
    ScoreCountMismatch { expected: usize, actual: usize },
    /// A retained score is not bound to the exact preregistered job at its position.
    JobMismatch { score_index: usize },
    /// A retained surrogate convergence score is non-finite.
    NonFiniteScore { score_index: usize, value: f64 },
    /// The frozen decision rule rejected a numerical input.
    Decision(StageU2DecisionError),
}

impl Display for U2ScoreReplayError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PairIndexOutOfRange { pair_index } => {
                write!(f, "U2 replay pair index {pair_index} is out of range")
            }
            Self::EmptySurrogateSet => {
                f.write_str("U2 replay requires at least one surrogate per null")
            }
            Self::ScoreCountOverflow => f.write_str("U2 replay score count overflow"),
            Self::ScoreCountMismatch { expected, actual } => write!(
                f,
                "U2 replay expected {expected} retained scores, got {actual}"
            ),
            Self::JobMismatch { score_index } => write!(
                f,
                "U2 replay score {score_index} is not bound to the expected preregistered job"
            ),
            Self::NonFiniteScore { score_index, value } => write!(
                f,
                "U2 replay score {score_index} must be finite, got {value}"
            ),
            Self::Decision(error) => write!(f, "U2 replay decision error: {error}"),
        }
    }
}

impl Error for U2ScoreReplayError {}

impl From<StageU2DecisionError> for U2ScoreReplayError {
    fn from(value: StageU2DecisionError) -> Self {
        Self::Decision(value)
    }
}

/// Recompute Stage U2 empirical p-values and the frozen decision from retained scores.
///
/// `scores` must contain exactly the jobs for `pair_index` in canonical order:
/// all shuffled-marginal repetitions `0..R`, followed by all phase-randomized
/// repetitions `0..R`, where `R = plan.surrogates_per_null`. Every retained job
/// identity, including its deterministic seed, is checked against `plan` before
/// its numerical score contributes to a p-value.
///
/// The empirical rule is the preregistered one-sided `+1` correction:
/// `(1 + count(score >= observed)) / (R + 1)`.
///
/// # Errors
///
/// Returns an error for an invalid pair index, zero/overflowing repetition count,
/// incomplete or reordered retained scores, non-finite retained scores, or any
/// numerical input rejected by the frozen decision rule.
pub fn replay_u2_statistics_from_scores(
    plan: U2ExecutionPlan,
    pair_index: usize,
    observed_convergence_score: f64,
    alpha: f64,
    scores: &[U2SurrogateScore],
) -> Result<ReplayedU2Statistics, U2ScoreReplayError> {
    if pair_index >= plan.pairs.len() {
        return Err(U2ScoreReplayError::PairIndexOutOfRange { pair_index });
    }
    if plan.surrogates_per_null == 0 {
        return Err(U2ScoreReplayError::EmptySurrogateSet);
    }

    let expected = plan
        .surrogates_per_null
        .checked_mul(2)
        .ok_or(U2ScoreReplayError::ScoreCountOverflow)?;
    if scores.len() != expected {
        return Err(U2ScoreReplayError::ScoreCountMismatch {
            expected,
            actual: scores.len(),
        });
    }

    let shuffled_at_least_as_extreme = count_extreme_scores(
        plan,
        pair_index,
        U2NullFamily::ShuffledMarginal,
        0,
        observed_convergence_score,
        scores,
    )?;
    let phase_at_least_as_extreme = count_extreme_scores(
        plan,
        pair_index,
        U2NullFamily::PhaseRandomizedSpectrum,
        plan.surrogates_per_null,
        observed_convergence_score,
        scores,
    )?;

    let p_shuffle = empirical_p(shuffled_at_least_as_extreme, plan.surrogates_per_null);
    let p_phase = empirical_p(phase_at_least_as_extreme, plan.surrogates_per_null);
    let decision = classify_stage_u2(observed_convergence_score, p_shuffle, p_phase, alpha)?;

    Ok(ReplayedU2Statistics {
        p_shuffle,
        p_phase,
        decision,
        shuffled_at_least_as_extreme,
        phase_at_least_as_extreme,
    })
}

fn count_extreme_scores(
    plan: U2ExecutionPlan,
    pair_index: usize,
    null_family: U2NullFamily,
    offset: usize,
    observed_convergence_score: f64,
    scores: &[U2SurrogateScore],
) -> Result<usize, U2ScoreReplayError> {
    let mut at_least_as_extreme = 0usize;
    for repetition in 0..plan.surrogates_per_null {
        let score_index = offset + repetition;
        let score = scores[score_index];
        let expected_job = plan
            .surrogate_job(pair_index, null_family, repetition)
            .ok_or(U2ScoreReplayError::JobMismatch { score_index })?;
        if score.job != expected_job {
            return Err(U2ScoreReplayError::JobMismatch { score_index });
        }
        if !score.convergence_score.is_finite() {
            return Err(U2ScoreReplayError::NonFiniteScore {
                score_index,
                value: score.convergence_score,
            });
        }
        if score.convergence_score >= observed_convergence_score {
            at_least_as_extreme += 1;
        }
    }
    Ok(at_least_as_extreme)
}

fn empirical_p(at_least_as_extreme: usize, repetitions: usize) -> f64 {
    (at_least_as_extreme + 1) as f64 / (repetitions + 1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::u2_plan::U2_FROZEN_PAIRS;

    fn plan(repetitions: usize) -> U2ExecutionPlan {
        U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: repetitions,
            seed_root: 0x51_52_53_54,
        }
    }

    fn scores_for(values: &[f64], plan: U2ExecutionPlan) -> Vec<U2SurrogateScore> {
        assert_eq!(values.len(), plan.surrogates_per_null * 2);
        let mut scores = Vec::with_capacity(values.len());
        for null_family in [
            U2NullFamily::ShuffledMarginal,
            U2NullFamily::PhaseRandomizedSpectrum,
        ] {
            for repetition in 0..plan.surrogates_per_null {
                let index = scores.len();
                scores.push(U2SurrogateScore {
                    job: plan.surrogate_job(0, null_family, repetition).unwrap(),
                    convergence_score: values[index],
                });
            }
        }
        scores
    }

    #[test]
    fn exact_plus_one_p_values_are_replayed_without_regenerating_surrogates() {
        let plan = plan(3);
        let scores = scores_for(&[0.6, 0.4, 0.8, 0.2, 0.3, 0.4], plan);
        let replay = replay_u2_statistics_from_scores(plan, 0, 0.5, 0.5, &scores).unwrap();

        assert_eq!(replay.shuffled_at_least_as_extreme, 2);
        assert_eq!(replay.phase_at_least_as_extreme, 0);
        assert_eq!(replay.p_shuffle, 0.75);
        assert_eq!(replay.p_phase, 0.25);
        assert_eq!(replay.decision, StageU2Decision::CompatibleWithMarginalNull);
    }

    #[test]
    fn equality_at_alpha_preserves_the_frozen_decision_boundary() {
        let plan = plan(3);
        let scores = scores_for(&[0.6, 0.4, 0.8, 0.2, 0.3, 0.4], plan);
        let replay = replay_u2_statistics_from_scores(plan, 0, 0.5, 0.75, &scores).unwrap();

        assert_eq!(replay.p_shuffle, 0.75);
        assert_eq!(replay.p_phase, 0.25);
        assert_eq!(replay.decision, StageU2Decision::CrossMechanismCandidate);
    }

    #[test]
    fn incomplete_reordered_and_nonfinite_score_archives_fail_closed() {
        let plan = plan(2);
        let scores = scores_for(&[0.5, 0.4, 0.3, 0.2], plan);

        let truncated = &scores[..scores.len() - 1];
        assert!(matches!(
            replay_u2_statistics_from_scores(plan, 0, 0.4, 0.05, truncated),
            Err(U2ScoreReplayError::ScoreCountMismatch { .. })
        ));

        let mut reordered = scores.clone();
        reordered.swap(0, 1);
        assert_eq!(
            replay_u2_statistics_from_scores(plan, 0, 0.4, 0.05, &reordered),
            Err(U2ScoreReplayError::JobMismatch { score_index: 0 })
        );

        let mut nonfinite = scores;
        nonfinite[2].convergence_score = f64::NAN;
        assert!(matches!(
            replay_u2_statistics_from_scores(plan, 0, 0.4, 0.05, &nonfinite),
            Err(U2ScoreReplayError::NonFiniteScore { score_index: 2, .. })
        ));
    }

    #[test]
    fn invalid_plan_shape_and_pair_index_fail_before_replay() {
        let empty = plan(0);
        assert_eq!(
            replay_u2_statistics_from_scores(empty, 0, 0.1, 0.05, &[]),
            Err(U2ScoreReplayError::EmptySurrogateSet)
        );

        let overflowing = plan(usize::MAX);
        assert_eq!(
            replay_u2_statistics_from_scores(overflowing, 0, 0.1, 0.05, &[]),
            Err(U2ScoreReplayError::ScoreCountOverflow)
        );

        let ordinary = plan(1);
        assert_eq!(
            replay_u2_statistics_from_scores(ordinary, U2_FROZEN_PAIRS.len(), 0.1, 0.05, &[]),
            Err(U2ScoreReplayError::PairIndexOutOfRange {
                pair_index: U2_FROZEN_PAIRS.len()
            })
        );
    }
}
