//! Provenance-bearing receipts for the NoiseLab calibration ladder.
//!
//! A receipt states where calibration evidence was produced. It does not verify
//! the scientific result and cannot by itself turn an experiment into proof.

use core::fmt;

use crate::calibration_ladder::{CalibrationEvidence, CalibrationStage};

const MAX_EVIDENCE_REF_BYTES: usize = 2048;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CalibrationReceipt {
    pub stage: CalibrationStage,
    pub repository: String,
    pub commit_sha: String,
    pub evidence_ref: String,
}

impl CalibrationReceipt {
    pub fn new(
        stage: CalibrationStage,
        repository: impl Into<String>,
        commit_sha: impl Into<String>,
        evidence_ref: impl Into<String>,
    ) -> Result<Self, CalibrationReceiptError> {
        let repository = repository.into();
        let commit_sha = commit_sha.into();
        let evidence_ref = evidence_ref.into();

        validate_repository(&repository)?;
        validate_commit(&commit_sha)?;
        if evidence_ref.trim().is_empty() || evidence_ref.len() > MAX_EVIDENCE_REF_BYTES {
            return Err(CalibrationReceiptError::InvalidEvidenceRef);
        }

        Ok(Self {
            stage,
            repository,
            commit_sha,
            evidence_ref,
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationReceiptError {
    InvalidRepository,
    InvalidCommitSha,
    InvalidEvidenceRef,
    DuplicateStage(CalibrationStage),
}

impl fmt::Display for CalibrationReceiptError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRepository => formatter.write_str("invalid owner/repository identity"),
            Self::InvalidCommitSha => formatter.write_str("invalid full Git commit SHA"),
            Self::InvalidEvidenceRef => {
                formatter.write_str("invalid calibration evidence reference")
            }
            Self::DuplicateStage(stage) => write!(formatter, "duplicate receipt for {stage:?}"),
        }
    }
}

impl std::error::Error for CalibrationReceiptError {}

/// Derive the existing outcome-blind gate flags from provenance-bearing receipts.
///
/// This validates receipt structure and uniqueness only. It deliberately does
/// not inspect or endorse the referenced scientific evidence.
pub fn evidence_from_receipts(
    receipts: &[CalibrationReceipt],
) -> Result<CalibrationEvidence, CalibrationReceiptError> {
    let mut evidence = CalibrationEvidence::none();
    for receipt in receipts {
        let slot = match receipt.stage {
            CalibrationStage::DrivenOscillator => &mut evidence.driven_oscillator,
            CalibrationStage::SemiconductorLaser => &mut evidence.semiconductor_laser,
            CalibrationStage::BistableKramers => &mut evidence.bistable_kramers,
            CalibrationStage::FitzHughNagumo => &mut evidence.fitzhugh_nagumo,
        };
        if *slot {
            return Err(CalibrationReceiptError::DuplicateStage(receipt.stage));
        }
        *slot = true;
    }
    Ok(evidence)
}

fn validate_repository(repository: &str) -> Result<(), CalibrationReceiptError> {
    let Some((owner, name)) = repository.split_once('/') else {
        return Err(CalibrationReceiptError::InvalidRepository);
    };
    if owner.is_empty()
        || name.is_empty()
        || name.contains('/')
        || !owner.bytes().all(valid_repo_byte)
        || !name.bytes().all(valid_repo_byte)
    {
        return Err(CalibrationReceiptError::InvalidRepository);
    }
    Ok(())
}

const fn valid_repo_byte(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.')
}

fn validate_commit(commit: &str) -> Result<(), CalibrationReceiptError> {
    if !matches!(commit.len(), 40 | 64) || !commit.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(CalibrationReceiptError::InvalidCommitSha);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receipt(stage: CalibrationStage) -> CalibrationReceipt {
        CalibrationReceipt::new(
            stage,
            "Memorithm/NoiseLab",
            "0123456789abcdef0123456789abcdef01234567",
            "results/calibration.json",
        )
        .expect("valid receipt")
    }

    #[test]
    fn receipts_drive_existing_outcome_blind_gate_flags() {
        let evidence = evidence_from_receipts(&[
            receipt(CalibrationStage::DrivenOscillator),
            receipt(CalibrationStage::SemiconductorLaser),
        ])
        .expect("evidence");
        assert!(evidence.driven_oscillator);
        assert!(evidence.semiconductor_laser);
        assert!(!evidence.bistable_kramers);
        assert!(!evidence.fitzhugh_nagumo);
    }

    #[test]
    fn duplicate_stage_fails_closed() {
        assert_eq!(
            evidence_from_receipts(&[
                receipt(CalibrationStage::DrivenOscillator),
                receipt(CalibrationStage::DrivenOscillator),
            ]),
            Err(CalibrationReceiptError::DuplicateStage(
                CalibrationStage::DrivenOscillator
            ))
        );
    }

    #[test]
    fn malformed_provenance_is_rejected() {
        assert_eq!(
            CalibrationReceipt::new(
                CalibrationStage::DrivenOscillator,
                "NoiseLab",
                "0123456789abcdef0123456789abcdef01234567",
                "evidence",
            ),
            Err(CalibrationReceiptError::InvalidRepository)
        );
        assert_eq!(
            CalibrationReceipt::new(
                CalibrationStage::DrivenOscillator,
                "Memorithm/NoiseLab",
                "deadbeef",
                "evidence",
            ),
            Err(CalibrationReceiptError::InvalidCommitSha)
        );
    }
}
