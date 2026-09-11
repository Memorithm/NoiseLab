//! Outcome-blind execution plan for the preregistered fluctuation-universality U2 panel.
//!
//! This module fixes source identities, unordered pair ordering and surrogate
//! seed derivation before any U2 outcome exists. It does not generate source
//! series, surrogates, scores or p-values.

use crate::u2_readiness::{
    U2Readiness, U2ReadinessError, U2_MIN_SURROGATES_PER_NULL, U2_UNORDERED_PAIRS,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum U2SourceFamily {
    DrivenDampedOscillator,
    SemiconductorLaser,
    BistableLangevin,
    FitzHughNagumo,
}

pub const U2_FROZEN_SOURCES: [U2SourceFamily; 4] = [
    U2SourceFamily::DrivenDampedOscillator,
    U2SourceFamily::SemiconductorLaser,
    U2SourceFamily::BistableLangevin,
    U2SourceFamily::FitzHughNagumo,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct U2Pair {
    pub left: U2SourceFamily,
    pub right: U2SourceFamily,
}

pub const U2_FROZEN_PAIRS: [U2Pair; U2_UNORDERED_PAIRS] = [
    U2Pair {
        left: U2SourceFamily::DrivenDampedOscillator,
        right: U2SourceFamily::SemiconductorLaser,
    },
    U2Pair {
        left: U2SourceFamily::DrivenDampedOscillator,
        right: U2SourceFamily::BistableLangevin,
    },
    U2Pair {
        left: U2SourceFamily::DrivenDampedOscillator,
        right: U2SourceFamily::FitzHughNagumo,
    },
    U2Pair {
        left: U2SourceFamily::SemiconductorLaser,
        right: U2SourceFamily::BistableLangevin,
    },
    U2Pair {
        left: U2SourceFamily::SemiconductorLaser,
        right: U2SourceFamily::FitzHughNagumo,
    },
    U2Pair {
        left: U2SourceFamily::BistableLangevin,
        right: U2SourceFamily::FitzHughNagumo,
    },
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum U2NullFamily {
    ShuffledMarginal,
    PhaseRandomizedSpectrum,
}

pub const U2_FROZEN_NULLS: [U2NullFamily; 2] = [
    U2NullFamily::ShuffledMarginal,
    U2NullFamily::PhaseRandomizedSpectrum,
];

/// One outcome-blind surrogate job fixed entirely by preregistered indices.
///
/// The descriptor contains no observed series, score, p-value or decision. It
/// can therefore be persisted before execution as an auditable work manifest.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct U2SurrogateJob {
    pub pair_index: usize,
    pub pair: U2Pair,
    pub null_family: U2NullFamily,
    pub repetition: usize,
    pub seed: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct U2ExecutionPlan {
    pub pairs: &'static [U2Pair; U2_UNORDERED_PAIRS],
    pub surrogates_per_null: usize,
    pub seed_root: u64,
}

impl U2ExecutionPlan {
    pub fn preregistered(seed_root: u64) -> Result<Self, U2ReadinessError> {
        U2Readiness::preregistered().validate()?;
        Ok(Self {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: U2_MIN_SURROGATES_PER_NULL,
            seed_root,
        })
    }

    /// Derive one deterministic seed without consulting any source or outcome.
    #[must_use]
    pub const fn surrogate_seed(
        self,
        pair_index: usize,
        null_family: U2NullFamily,
        repetition: usize,
    ) -> Option<u64> {
        if pair_index >= U2_UNORDERED_PAIRS || repetition >= self.surrogates_per_null {
            return None;
        }
        let null_tag = match null_family {
            U2NullFamily::ShuffledMarginal => 0x5348_5546_464c_4501,
            U2NullFamily::PhaseRandomizedSpectrum => 0x5048_4153_4500_0002,
        };
        let pair_tag = (pair_index as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let repetition_tag = (repetition as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        Some(self.seed_root ^ null_tag ^ pair_tag ^ repetition_tag)
    }

    /// Materialize one preregistered surrogate job without observing outcomes.
    #[must_use]
    pub const fn surrogate_job(
        self,
        pair_index: usize,
        null_family: U2NullFamily,
        repetition: usize,
    ) -> Option<U2SurrogateJob> {
        let seed = match self.surrogate_seed(pair_index, null_family, repetition) {
            Some(seed) => seed,
            None => return None,
        };
        Some(U2SurrogateJob {
            pair_index,
            pair: self.pairs[pair_index],
            null_family,
            repetition,
            seed,
        })
    }

    /// Exact number of preregistered surrogate jobs in the frozen U2 panel.
    #[must_use]
    pub const fn surrogate_job_count(self) -> usize {
        U2_UNORDERED_PAIRS * U2_FROZEN_NULLS.len() * self.surrogates_per_null
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn frozen_panel_contains_every_unordered_pair_once() {
        assert_eq!(U2_FROZEN_PAIRS.len(), U2_UNORDERED_PAIRS);
        let unique = U2_FROZEN_PAIRS
            .iter()
            .map(|pair| (pair.left, pair.right))
            .collect::<BTreeSet<_>>();
        assert_eq!(unique.len(), U2_UNORDERED_PAIRS);
        assert!(U2_FROZEN_PAIRS.iter().all(|pair| pair.left < pair.right));
    }

    #[test]
    fn seed_plan_is_deterministic_and_null_specific() {
        let plan = U2ExecutionPlan::preregistered(0x1234_5678).unwrap();
        let a = plan
            .surrogate_seed(2, U2NullFamily::ShuffledMarginal, 17)
            .unwrap();
        let b = plan
            .surrogate_seed(2, U2NullFamily::ShuffledMarginal, 17)
            .unwrap();
        let spectral = plan
            .surrogate_seed(2, U2NullFamily::PhaseRandomizedSpectrum, 17)
            .unwrap();
        assert_eq!(a, b);
        assert_ne!(a, spectral);
    }

    #[test]
    fn job_manifest_is_complete_without_outcome_data() {
        let plan = U2ExecutionPlan::preregistered(7).unwrap();
        assert_eq!(
            plan.surrogate_job_count(),
            U2_UNORDERED_PAIRS * 2 * U2_MIN_SURROGATES_PER_NULL
        );
        let job = plan
            .surrogate_job(5, U2NullFamily::PhaseRandomizedSpectrum, 198)
            .unwrap();
        assert_eq!(job.pair_index, 5);
        assert_eq!(job.pair, U2_FROZEN_PAIRS[5]);
        assert_eq!(job.repetition, 198);
        assert_eq!(
            Some(job.seed),
            plan.surrogate_seed(5, U2NullFamily::PhaseRandomizedSpectrum, 198)
        );
    }

    #[test]
    fn out_of_contract_indices_are_rejected() {
        let plan = U2ExecutionPlan::preregistered(1).unwrap();
        assert_eq!(
            plan.surrogate_seed(U2_UNORDERED_PAIRS, U2NullFamily::ShuffledMarginal, 0),
            None
        );
        assert_eq!(
            plan.surrogate_seed(
                0,
                U2NullFamily::ShuffledMarginal,
                U2_MIN_SURROGATES_PER_NULL
            ),
            None
        );
        assert_eq!(
            plan.surrogate_job(U2_UNORDERED_PAIRS, U2NullFamily::ShuffledMarginal, 0),
            None
        );
    }
}
