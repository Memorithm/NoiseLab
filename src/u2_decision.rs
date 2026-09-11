//! Deterministic decision gate for the preregistered fluctuation-universality U2 panel.
//!
//! This module deliberately does not generate spectral surrogates. U2 execution remains
//! blocked until the SciRust phase-randomized surrogate primitive is merged, requalified,
//! and pinned by immutable revision. The code here only freezes the already-preregistered
//! classification rule so later execution cannot silently reinterpret observed outcomes.

use std::error::Error;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageU2Decision {
    /// The declared first-to-last descriptor distance did not decrease.
    NoObservedConvergence,
    /// Positive convergence remained compatible with shuffled empirical-marginal controls.
    CompatibleWithMarginalNull,
    /// Shuffled controls were separated, but the phase-randomized spectral null was not.
    SpectrumExplainedCandidate,
    /// Both preregistered null families were separated. This remains exploratory only.
    CrossMechanismCandidate,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StageU2DecisionError {
    NonFiniteConvergenceScore(f64),
    InvalidPValue { name: &'static str, value: f64 },
    InvalidAlpha(f64),
}

impl Display for StageU2DecisionError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFiniteConvergenceScore(value) => {
                write!(f, "convergence score must be finite, got {value}")
            }
            Self::InvalidPValue { name, value } => {
                write!(f, "{name} must be finite and in [0, 1], got {value}")
            }
            Self::InvalidAlpha(value) => {
                write!(f, "alpha must be finite and in (0, 1], got {value}")
            }
        }
    }
}

impl Error for StageU2DecisionError {}

/// Apply the frozen U2 Stage-0 classification rule.
///
/// `p_shuffle` is the one-sided empirical p-value against shuffled empirical-marginal
/// controls and `p_phase` is the one-sided empirical p-value against phase-randomized
/// spectral controls. The caller remains responsible for producing those values under the
/// preregistered surrogate protocol. Invalid numerical inputs fail closed rather than being
/// coerced into a scientific label.
pub fn classify_stage_u2(
    convergence_score: f64,
    p_shuffle: f64,
    p_phase: f64,
    alpha: f64,
) -> Result<StageU2Decision, StageU2DecisionError> {
    if !convergence_score.is_finite() {
        return Err(StageU2DecisionError::NonFiniteConvergenceScore(
            convergence_score,
        ));
    }
    validate_p_value("p_shuffle", p_shuffle)?;
    validate_p_value("p_phase", p_phase)?;
    if !alpha.is_finite() || alpha <= 0.0 || alpha > 1.0 {
        return Err(StageU2DecisionError::InvalidAlpha(alpha));
    }

    if convergence_score <= 0.0 {
        Ok(StageU2Decision::NoObservedConvergence)
    } else if p_shuffle > alpha {
        Ok(StageU2Decision::CompatibleWithMarginalNull)
    } else if p_phase > alpha {
        Ok(StageU2Decision::SpectrumExplainedCandidate)
    } else {
        Ok(StageU2Decision::CrossMechanismCandidate)
    }
}

fn validate_p_value(name: &'static str, value: f64) -> Result<(), StageU2DecisionError> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(StageU2DecisionError::InvalidPValue { name, value });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const ALPHA: f64 = 0.05;

    #[test]
    fn non_positive_score_has_priority_over_null_separation() {
        assert_eq!(
            classify_stage_u2(0.0, 0.001, 0.001, ALPHA).unwrap(),
            StageU2Decision::NoObservedConvergence
        );
        assert_eq!(
            classify_stage_u2(-0.1, 0.001, 0.001, ALPHA).unwrap(),
            StageU2Decision::NoObservedConvergence
        );
    }

    #[test]
    fn positive_score_compatible_with_shuffle_stops_at_marginal_null() {
        assert_eq!(
            classify_stage_u2(0.2, 0.051, 0.001, ALPHA).unwrap(),
            StageU2Decision::CompatibleWithMarginalNull
        );
    }

    #[test]
    fn phase_null_prevents_cross_mechanism_label() {
        assert_eq!(
            classify_stage_u2(0.2, ALPHA, 0.051, ALPHA).unwrap(),
            StageU2Decision::SpectrumExplainedCandidate
        );
    }

    #[test]
    fn equality_at_alpha_counts_as_separation_per_preregistration() {
        assert_eq!(
            classify_stage_u2(0.2, ALPHA, ALPHA, ALPHA).unwrap(),
            StageU2Decision::CrossMechanismCandidate
        );
    }

    #[test]
    fn invalid_numerical_inputs_fail_closed() {
        assert!(matches!(
            classify_stage_u2(f64::NAN, 0.1, 0.1, ALPHA),
            Err(StageU2DecisionError::NonFiniteConvergenceScore(_))
        ));
        assert!(matches!(
            classify_stage_u2(0.1, -0.01, 0.1, ALPHA),
            Err(StageU2DecisionError::InvalidPValue {
                name: "p_shuffle",
                ..
            })
        ));
        assert!(matches!(
            classify_stage_u2(0.1, 0.1, 1.01, ALPHA),
            Err(StageU2DecisionError::InvalidPValue {
                name: "p_phase",
                ..
            })
        ));
        assert!(matches!(
            classify_stage_u2(0.1, 0.1, 0.1, 0.0),
            Err(StageU2DecisionError::InvalidAlpha(_))
        ));
    }
}
