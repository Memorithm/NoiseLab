//! Outcome-blind catalog of immutable protocol provenance.
//!
//! This module records only identities needed to reproduce or audit frozen
//! protocols. It does not execute an experiment, inspect an outcome, or imply
//! that a calibration or U2 hypothesis has been supported.

use crate::preregistered::{
    BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA, FHN_COHERENCE_STAGE0_PROTOCOL_BLOB_SHA,
};
use crate::u2_readiness::U2_SCIRUST_REVISION;

/// Kind of immutable Git identity carried by a protocol dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GitIdentityKind {
    /// SHA-1 Git blob identity for a frozen protocol document.
    Blob,
    /// Full Git commit identity for an external code dependency.
    Commit,
}

/// One outcome-blind immutable dependency of a preregistered protocol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProtocolProvenance {
    /// Stable human-readable protocol/dependency key.
    pub key: &'static str,
    /// Repository that owns the referenced object.
    pub repository: &'static str,
    /// Git object kind.
    pub kind: GitIdentityKind,
    /// Exact lowercase hexadecimal Git SHA-1 identity.
    pub sha: &'static str,
}

impl ProtocolProvenance {
    /// Validate structural provenance without reading referenced outcomes.
    pub fn validate(self) -> Result<(), ProtocolProvenanceError> {
        if self.key.is_empty() || self.repository.is_empty() {
            return Err(ProtocolProvenanceError::EmptyIdentityField);
        }
        if self.sha.len() != 40
            || !self
                .sha
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(ProtocolProvenanceError::InvalidGitSha1);
        }
        Ok(())
    }
}

/// Structural provenance validation errors.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProtocolProvenanceError {
    /// A stable key or repository identity is absent.
    EmptyIdentityField,
    /// A declared immutable Git identity is not a full lowercase SHA-1.
    InvalidGitSha1,
}

/// Frozen outcome-blind provenance catalog for the active calibration ladder and U2 preflight.
pub const PROTOCOL_PROVENANCE: [ProtocolProvenance; 3] = [
    ProtocolProvenance {
        key: "bistable-stage0-v2",
        repository: "Memorithm/NoiseLab",
        kind: GitIdentityKind::Blob,
        sha: BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA,
    },
    ProtocolProvenance {
        key: "fhn-coherence-stage0",
        repository: "Memorithm/NoiseLab",
        kind: GitIdentityKind::Blob,
        sha: FHN_COHERENCE_STAGE0_PROTOCOL_BLOB_SHA,
    },
    ProtocolProvenance {
        key: "u2-scirust-spectral-null",
        repository: "Memorithm/scirust",
        kind: GitIdentityKind::Commit,
        sha: U2_SCIRUST_REVISION,
    },
];

/// Validate every immutable identity in the active outcome-blind catalog.
pub fn validate_protocol_provenance() -> Result<(), ProtocolProvenanceError> {
    for record in PROTOCOL_PROVENANCE {
        record.validate()?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_catalog_is_structurally_valid() {
        assert_eq!(validate_protocol_provenance(), Ok(()));
        assert_eq!(PROTOCOL_PROVENANCE.len(), 3);
    }

    #[test]
    fn catalog_pins_expected_external_scirust_revision() {
        let u2 = PROTOCOL_PROVENANCE
            .iter()
            .find(|record| record.key == "u2-scirust-spectral-null")
            .unwrap();
        assert_eq!(u2.repository, "Memorithm/scirust");
        assert_eq!(u2.kind, GitIdentityKind::Commit);
        assert_eq!(u2.sha, U2_SCIRUST_REVISION);
    }

    #[test]
    fn malformed_git_identity_fails_closed() {
        let malformed = ProtocolProvenance {
            key: "bad",
            repository: "Memorithm/NoiseLab",
            kind: GitIdentityKind::Blob,
            sha: "DEADBEEF",
        };
        assert_eq!(
            malformed.validate(),
            Err(ProtocolProvenanceError::InvalidGitSha1)
        );
    }
}
