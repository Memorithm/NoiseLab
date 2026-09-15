//! Dual-null Stage U2 analysis for one frozen cross-mechanism residual pair.
//!
//! For each pair the observed multiscale convergence score is computed once.
//! Surrogate p-values are then obtained independently for:
//!
//! - shuffled empirical-marginal controls (Fisher–Yates / `SplitMix64`);
//! - phase-randomized spectral controls via [`spectral_phase_null`].
//!
//! Seeds come from the outcome-blind [`U2SurrogateJob`] manifest. Empirical
//! p-values use the standard one-sided `+1` correction
//! `(1 + count(surr >= observed)) / (R + 1)`. Protocol / numerical failures are
//! retained rather than dropped.

use crate::spectral_null::spectral_phase_null;
use crate::u2_decision::{classify_stage_u2, StageU2Decision, StageU2DecisionError};
use crate::u2_plan::{U2NullFamily, U2Pair, U2SurrogateJob};
use crate::universality::{
    observed_multiscale_comparison, ObservedMultiscaleComparison, UniversalityError,
};
use scirust_signal::surrogate::SurrogateError;
use scirust_sim::SplitMix64;
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Tag mixed into the job seed for the right-hand spectral surrogate.
///
/// Rule: `seed_left = job.seed`, `seed_right = job.seed ^ SPECTRAL_RIGHT_SEED_TAG`.
/// Both series of a pair are transformed independently; the tag keeps the pair
/// deterministic without sharing an identical phase draw.
pub const SPECTRAL_RIGHT_SEED_TAG: u64 = 0x9e37_79b9_7f4a_7c15;

/// Exact score produced by one preregistered surrogate job.
///
/// This is execution evidence, not a decision or a claim. The job preserves the
/// pair/null/repetition/seed identity fixed before outcomes were observed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct U2SurrogateScore {
    pub job: U2SurrogateJob,
    pub convergence_score: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StageU2PairAnalysis {
    pub pair_index: usize,
    pub pair: U2Pair,
    pub fine_scale_distance: f64,
    pub terminal_scale_distance: f64,
    pub convergence_score: f64,
    pub p_shuffle: Option<f64>,
    pub p_phase: Option<f64>,
    pub decision: Option<StageU2Decision>,
    /// Per-job surrogate scores in deterministic order: shuffled then phase.
    pub surrogate_scores: Vec<U2SurrogateScore>,
    /// Non-empty when a protocol / numerical failure prevented a full decision.
    pub protocol_error: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StageU2AnalysisError {
    Universality(UniversalityError),
    Decision(StageU2DecisionError),
    Spectral(SurrogateError),
    MissingJobs {
        pair_index: usize,
        null_family: U2NullFamily,
    },
    EmptyScales,
    AllocationFailed,
}

impl Display for StageU2AnalysisError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Universality(error) => write!(f, "U2 universality error: {error}"),
            Self::Decision(error) => write!(f, "U2 decision error: {error}"),
            Self::Spectral(error) => write!(f, "U2 spectral-null error: {error}"),
            Self::MissingJobs {
                pair_index,
                null_family,
            } => write!(
                f,
                "no surrogate jobs for pair {pair_index} under null {null_family:?}"
            ),
            Self::EmptyScales => f.write_str("U2 analysis requires a non-empty scale grid"),
            Self::AllocationFailed => {
                f.write_str("unable to allocate U2 per-surrogate score evidence")
            }
        }
    }
}

impl Error for StageU2AnalysisError {}

impl From<UniversalityError> for StageU2AnalysisError {
    fn from(value: UniversalityError) -> Self {
        Self::Universality(value)
    }
}

impl From<StageU2DecisionError> for StageU2AnalysisError {
    fn from(value: StageU2DecisionError) -> Self {
        Self::Decision(value)
    }
}

impl From<SurrogateError> for StageU2AnalysisError {
    fn from(value: SurrogateError) -> Self {
        Self::Spectral(value)
    }
}

/// Inputs for one dual-null pair analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StageU2PairRequest<'a> {
    pub pair_index: usize,
    pub pair: U2Pair,
    pub series_a: &'a [f64],
    pub series_b: &'a [f64],
    pub scales: &'a [usize],
    pub alpha: f64,
    pub jobs: &'a [U2SurrogateJob],
    pub retain_protocol_failures: bool,
}

/// Analyze one residual pair under both preregistered null families.
///
/// `jobs` must contain every surrogate job for `pair_index` (both nulls). Jobs
/// for other pairs are ignored. Failures are returned as `protocol_error` on the
/// analysis struct when `retain_protocol_failures` is true; otherwise they
/// propagate as `Err`.
pub fn analyze_u2_pair(
    request: StageU2PairRequest<'_>,
) -> Result<StageU2PairAnalysis, StageU2AnalysisError> {
    let StageU2PairRequest {
        pair_index,
        pair,
        series_a,
        series_b,
        scales,
        alpha,
        jobs,
        retain_protocol_failures,
    } = request;

    if scales.is_empty() {
        return Err(StageU2AnalysisError::EmptyScales);
    }

    let observed = match observed_multiscale_comparison(series_a, series_b, scales) {
        Ok(value) => value,
        Err(error) => {
            return finish_protocol_failure(
                pair_index,
                pair,
                retain_protocol_failures,
                StageU2AnalysisError::Universality(error),
            );
        }
    };

    let shuffle_jobs = jobs_for(pair_index, U2NullFamily::ShuffledMarginal, jobs);
    let phase_jobs = jobs_for(pair_index, U2NullFamily::PhaseRandomizedSpectrum, jobs);
    if shuffle_jobs.is_empty() {
        return finish_protocol_failure(
            pair_index,
            pair,
            retain_protocol_failures,
            StageU2AnalysisError::MissingJobs {
                pair_index,
                null_family: U2NullFamily::ShuffledMarginal,
            },
        );
    }
    if phase_jobs.is_empty() {
        return finish_protocol_failure(
            pair_index,
            pair,
            retain_protocol_failures,
            StageU2AnalysisError::MissingJobs {
                pair_index,
                null_family: U2NullFamily::PhaseRandomizedSpectrum,
            },
        );
    }

    let (p_shuffle, mut surrogate_scores) =
        match shuffle_p_value(series_a, series_b, scales, &observed, &shuffle_jobs) {
            Ok(value) => value,
            Err(error) => {
                return finish_protocol_failure(pair_index, pair, retain_protocol_failures, error);
            }
        };
    let (p_phase, phase_scores) =
        match phase_p_value(series_a, series_b, scales, &observed, &phase_jobs) {
            Ok(value) => value,
            Err(error) => {
                return finish_protocol_failure(pair_index, pair, retain_protocol_failures, error);
            }
        };
    surrogate_scores
        .try_reserve(phase_scores.len())
        .map_err(|_| StageU2AnalysisError::AllocationFailed)?;
    surrogate_scores.extend(phase_scores);

    let decision = match classify_stage_u2(observed.convergence_score, p_shuffle, p_phase, alpha) {
        Ok(value) => value,
        Err(error) => {
            return finish_protocol_failure(
                pair_index,
                pair,
                retain_protocol_failures,
                StageU2AnalysisError::Decision(error),
            );
        }
    };

    Ok(StageU2PairAnalysis {
        pair_index,
        pair,
        fine_scale_distance: observed.fine_scale_distance,
        terminal_scale_distance: observed.terminal_scale_distance,
        convergence_score: observed.convergence_score,
        p_shuffle: Some(p_shuffle),
        p_phase: Some(p_phase),
        decision: Some(decision),
        surrogate_scores,
        protocol_error: None,
    })
}

fn finish_protocol_failure(
    pair_index: usize,
    pair: U2Pair,
    retain: bool,
    error: StageU2AnalysisError,
) -> Result<StageU2PairAnalysis, StageU2AnalysisError> {
    if retain {
        Ok(StageU2PairAnalysis {
            pair_index,
            pair,
            fine_scale_distance: f64::NAN,
            terminal_scale_distance: f64::NAN,
            convergence_score: f64::NAN,
            p_shuffle: None,
            p_phase: None,
            decision: None,
            surrogate_scores: Vec::new(),
            protocol_error: Some(error.to_string()),
        })
    } else {
        Err(error)
    }
}

fn jobs_for(
    pair_index: usize,
    null_family: U2NullFamily,
    jobs: &[U2SurrogateJob],
) -> Vec<U2SurrogateJob> {
    jobs.iter()
        .copied()
        .filter(|job| job.pair_index == pair_index && job.null_family == null_family)
        .collect()
}

fn score_buffer(capacity: usize) -> Result<Vec<U2SurrogateScore>, StageU2AnalysisError> {
    let mut scores = Vec::new();
    scores
        .try_reserve_exact(capacity)
        .map_err(|_| StageU2AnalysisError::AllocationFailed)?;
    Ok(scores)
}

fn shuffle_p_value(
    series_a: &[f64],
    series_b: &[f64],
    scales: &[usize],
    observed: &ObservedMultiscaleComparison,
    jobs: &[U2SurrogateJob],
) -> Result<(f64, Vec<U2SurrogateScore>), StageU2AnalysisError> {
    let mut at_least_as_extreme = 0usize;
    let mut scores = score_buffer(jobs.len())?;
    for job in jobs {
        let mut surrogate_a = series_a.to_vec();
        let mut surrogate_b = series_b.to_vec();
        let mut rng = SplitMix64::new(job.seed);
        fisher_yates_shuffle(&mut surrogate_a, &mut rng);
        fisher_yates_shuffle(&mut surrogate_b, &mut rng);
        let score =
            observed_multiscale_comparison(&surrogate_a, &surrogate_b, scales)?.convergence_score;
        if score >= observed.convergence_score {
            at_least_as_extreme += 1;
        }
        scores.push(U2SurrogateScore {
            job: *job,
            convergence_score: score,
        });
    }
    Ok((empirical_p(at_least_as_extreme, jobs.len()), scores))
}

fn phase_p_value(
    series_a: &[f64],
    series_b: &[f64],
    scales: &[usize],
    observed: &ObservedMultiscaleComparison,
    jobs: &[U2SurrogateJob],
) -> Result<(f64, Vec<U2SurrogateScore>), StageU2AnalysisError> {
    let mut at_least_as_extreme = 0usize;
    let mut scores = score_buffer(jobs.len())?;
    for job in jobs {
        let seed_left = job.seed;
        let seed_right = job.seed ^ SPECTRAL_RIGHT_SEED_TAG;
        let surrogate_a = spectral_phase_null(series_a, seed_left)?;
        let surrogate_b = spectral_phase_null(series_b, seed_right)?;
        let score =
            observed_multiscale_comparison(&surrogate_a, &surrogate_b, scales)?.convergence_score;
        if score >= observed.convergence_score {
            at_least_as_extreme += 1;
        }
        scores.push(U2SurrogateScore {
            job: *job,
            convergence_score: score,
        });
    }
    Ok((empirical_p(at_least_as_extreme, jobs.len()), scores))
}

fn empirical_p(at_least_as_extreme: usize, repetitions: usize) -> f64 {
    (at_least_as_extreme + 1) as f64 / (repetitions + 1) as f64
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
    use crate::u2_plan::{U2ExecutionPlan, U2SourceFamily, U2_FROZEN_PAIRS};

    fn tiny_scales() -> Vec<usize> {
        vec![1, 2, 4]
    }

    fn synthetic_pair() -> (Vec<f64>, Vec<f64>) {
        let a: Vec<f64> = (0..256).map(|i| (0.11 * i as f64).sin()).collect();
        let b: Vec<f64> = (0..256)
            .map(|i| (0.07 * i as f64).cos() + 0.05 * (i as f64))
            .collect();
        (a, b)
    }

    #[test]
    fn dual_null_is_seed_reproducible_and_p_in_unit_interval() {
        let (a, b) = synthetic_pair();
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: 19,
            seed_root: 0xABCDu64,
        };
        let jobs = materialize_u2_manifest(plan);
        let request = StageU2PairRequest {
            pair_index: 0,
            pair: U2_FROZEN_PAIRS[0],
            series_a: &a,
            series_b: &b,
            scales: &tiny_scales(),
            alpha: 0.05,
            jobs: &jobs,
            retain_protocol_failures: true,
        };
        let first = analyze_u2_pair(request).unwrap();
        let second = analyze_u2_pair(request).unwrap();
        assert_eq!(first, second);
        assert!(first.protocol_error.is_none());
        let p_shuffle = first.p_shuffle.unwrap();
        let p_phase = first.p_phase.unwrap();
        assert!((0.0..=1.0).contains(&p_shuffle));
        assert!((0.0..=1.0).contains(&p_phase));
        assert!(first.decision.is_some());
        assert_eq!(first.pair.left, U2SourceFamily::DrivenDampedOscillator);
        assert_eq!(first.surrogate_scores.len(), 38);
        assert!(
            first.surrogate_scores[..19]
                .iter()
                .all(|score| score.job.null_family == U2NullFamily::ShuffledMarginal)
        );
        assert!(
            first.surrogate_scores[19..]
                .iter()
                .all(|score| score.job.null_family == U2NullFamily::PhaseRandomizedSpectrum)
        );
        assert!(
            first
                .surrogate_scores
                .iter()
                .all(|score| score.convergence_score.is_finite())
        );

        let shuffle_extreme = first.surrogate_scores[..19]
            .iter()
            .filter(|score| score.convergence_score >= first.convergence_score)
            .count();
        let phase_extreme = first.surrogate_scores[19..]
            .iter()
            .filter(|score| score.convergence_score >= first.convergence_score)
            .count();
        assert_eq!(p_shuffle, empirical_p(shuffle_extreme, 19));
        assert_eq!(p_phase, empirical_p(phase_extreme, 19));
    }

    #[test]
    fn spectral_surrogate_same_job_seed_is_deterministic() {
        let signal: Vec<f64> = (0..64).map(|i| (0.2 * i as f64).sin()).collect();
        let a = spectral_phase_null(&signal, 9).unwrap();
        let b = spectral_phase_null(&signal, 9).unwrap();
        assert_eq!(
            a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn fail_closed_inputs_are_retained_as_protocol_errors() {
        let bad = vec![f64::NAN; 256];
        let ok: Vec<f64> = (0..256).map(|i| i as f64).collect();
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: 19,
            seed_root: 1,
        };
        let jobs = materialize_u2_manifest(plan);
        let result = analyze_u2_pair(StageU2PairRequest {
            pair_index: 0,
            pair: U2_FROZEN_PAIRS[0],
            series_a: &bad,
            series_b: &ok,
            scales: &tiny_scales(),
            alpha: 0.05,
            jobs: &jobs,
            retain_protocol_failures: true,
        })
        .unwrap();
        assert!(result.protocol_error.is_some());
        assert!(result.decision.is_none());
        assert!(result.surrogate_scores.is_empty());
    }
}
