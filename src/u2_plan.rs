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
    }
}
