//! Outcome-blind readiness checks for the preregistered fluctuation-universality U2 stage.
//!
//! This module never runs U2 and never inspects U2 outcomes. It only verifies
//! that frozen structural prerequisites are satisfied before an executor may
//! start the preregistered six-pair panel.

use core::fmt;

pub const U2_SCIRUST_REVISION: &str = "0e2eaccac631b689f97c242c47bad11d433847d9";
pub const U2_SOURCE_FAMILIES: usize = 4;
pub const U2_UNORDERED_PAIRS: usize = 6;
pub const U2_MIN_SURROGATES_PER_NULL: usize = 199;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum U2ReadinessError {
    WrongScirustRevision,
    SpectralNullNotQualified,
    WrongSourceFamilyCount,
    WrongPairCount,
    TooFewSurrogates,
}

impl fmt::Display for U2ReadinessError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            Self::WrongScirustRevision => "U2 requires the exact preregistered SciRust revision",
            Self::SpectralNullNotQualified => "NoiseLab spectral-null consumption is not qualified",
            Self::WrongSourceFamilyCount => "U2 requires exactly four frozen source families",
            Self::WrongPairCount => "U2 requires all six unordered source-family pairs",
            Self::TooFewSurrogates => "U2 requires at least 199 repetitions per null family",
        };
        formatter.write_str(message)
    }
}

impl std::error::Error for U2ReadinessError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct U2Readiness {
    pub scirust_revision: &'static str,
    pub spectral_null_qualified: bool,
    pub source_family_count: usize,
    pub unordered_pair_count: usize,
    pub surrogates_per_null: usize,
}

impl U2Readiness {
    #[must_use]
    pub const fn preregistered() -> Self {
        Self {
            scirust_revision: U2_SCIRUST_REVISION,
            // PR #26 qualified the NoiseLab adapter against the immutable
            // revision above. This is a code/provenance fact, not U2 evidence.
            spectral_null_qualified: true,
            source_family_count: U2_SOURCE_FAMILIES,
            unordered_pair_count: U2_UNORDERED_PAIRS,
            surrogates_per_null: U2_MIN_SURROGATES_PER_NULL,
        }
    }

    pub fn validate(self) -> Result<(), U2ReadinessError> {
        if self.scirust_revision != U2_SCIRUST_REVISION {
            return Err(U2ReadinessError::WrongScirustRevision);
        }
        if !self.spectral_null_qualified {
            return Err(U2ReadinessError::SpectralNullNotQualified);
        }
        if self.source_family_count != U2_SOURCE_FAMILIES {
            return Err(U2ReadinessError::WrongSourceFamilyCount);
        }
        if self.unordered_pair_count != U2_UNORDERED_PAIRS {
            return Err(U2ReadinessError::WrongPairCount);
        }
        if self.surrogates_per_null < U2_MIN_SURROGATES_PER_NULL {
            return Err(U2ReadinessError::TooFewSurrogates);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preregistered_contract_is_ready_without_running_u2() {
        assert_eq!(U2Readiness::preregistered().validate(), Ok(()));
    }

    #[test]
    fn revision_drift_fails_closed() {
        let mut readiness = U2Readiness::preregistered();
        readiness.scirust_revision = "deadbeef";
        assert_eq!(readiness.validate(), Err(U2ReadinessError::WrongScirustRevision));
    }

    #[test]
    fn missing_qualification_fails_closed() {
        let mut readiness = U2Readiness::preregistered();
        readiness.spectral_null_qualified = false;
        assert_eq!(readiness.validate(), Err(U2ReadinessError::SpectralNullNotQualified));
    }

    #[test]
    fn underpowered_null_family_fails_closed() {
        let mut readiness = U2Readiness::preregistered();
        readiness.surrogates_per_null = U2_MIN_SURROGATES_PER_NULL - 1;
        assert_eq!(readiness.validate(), Err(U2ReadinessError::TooFewSurrogates));
    }
}
