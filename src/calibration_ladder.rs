//! Outcome-blind calibration-order gate for NoiseLab research programmes.
//!
//! The gate records only whether prerequisite calibration evidence has been
//! supplied by the caller. It does not inspect experimental outcomes, infer a
//! resonance mechanism, or promote code availability into scientific evidence.

use core::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationStage {
    DrivenOscillator,
    SemiconductorLaser,
    BistableKramers,
    FitzHughNagumo,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CalibrationEvidence {
    pub driven_oscillator: bool,
    pub semiconductor_laser: bool,
    pub bistable_kramers: bool,
    pub fitzhugh_nagumo: bool,
}

impl CalibrationEvidence {
    #[must_use]
    pub const fn none() -> Self {
        Self {
            driven_oscillator: false,
            semiconductor_laser: false,
            bistable_kramers: false,
            fitzhugh_nagumo: false,
        }
    }

    #[must_use]
    pub const fn completed_through(self, stage: CalibrationStage) -> bool {
        match stage {
            CalibrationStage::DrivenOscillator => self.driven_oscillator,
            CalibrationStage::SemiconductorLaser => {
                self.driven_oscillator && self.semiconductor_laser
            }
            CalibrationStage::BistableKramers => {
                self.driven_oscillator && self.semiconductor_laser && self.bistable_kramers
            }
            CalibrationStage::FitzHughNagumo => {
                self.driven_oscillator
                    && self.semiconductor_laser
                    && self.bistable_kramers
                    && self.fitzhugh_nagumo
            }
        }
    }

    pub fn validate_complex_system_entry(self) -> Result<(), CalibrationGateError> {
        if !self.driven_oscillator {
            return Err(CalibrationGateError::Missing(
                CalibrationStage::DrivenOscillator,
            ));
        }
        if !self.semiconductor_laser {
            return Err(CalibrationGateError::Missing(
                CalibrationStage::SemiconductorLaser,
            ));
        }
        if !self.bistable_kramers {
            return Err(CalibrationGateError::Missing(
                CalibrationStage::BistableKramers,
            ));
        }
        if !self.fitzhugh_nagumo {
            return Err(CalibrationGateError::Missing(
                CalibrationStage::FitzHughNagumo,
            ));
        }
        Ok(())
    }
}

impl Default for CalibrationEvidence {
    fn default() -> Self {
        Self::none()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationGateError {
    Missing(CalibrationStage),
}

impl fmt::Display for CalibrationGateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Missing(stage) => write!(formatter, "missing calibration evidence for {stage:?}"),
        }
    }
}

impl std::error::Error for CalibrationGateError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complex_system_entry_is_closed_by_default() {
        assert_eq!(
            CalibrationEvidence::none().validate_complex_system_entry(),
            Err(CalibrationGateError::Missing(
                CalibrationStage::DrivenOscillator
            ))
        );
    }

    #[test]
    fn gate_reports_first_missing_stage_in_frozen_order() {
        let evidence = CalibrationEvidence {
            driven_oscillator: true,
            semiconductor_laser: true,
            bistable_kramers: false,
            fitzhugh_nagumo: true,
        };
        assert_eq!(
            evidence.validate_complex_system_entry(),
            Err(CalibrationGateError::Missing(
                CalibrationStage::BistableKramers
            ))
        );
    }

    #[test]
    fn all_four_explicit_prerequisites_allow_complex_entry() {
        let evidence = CalibrationEvidence {
            driven_oscillator: true,
            semiconductor_laser: true,
            bistable_kramers: true,
            fitzhugh_nagumo: true,
        };
        assert_eq!(evidence.validate_complex_system_entry(), Ok(()));
        assert!(evidence.completed_through(CalibrationStage::FitzHughNagumo));
    }

    #[test]
    fn later_flags_do_not_skip_earlier_stages() {
        let evidence = CalibrationEvidence {
            driven_oscillator: false,
            semiconductor_laser: true,
            bistable_kramers: true,
            fitzhugh_nagumo: true,
        };
        assert!(!evidence.completed_through(CalibrationStage::SemiconductorLaser));
        assert!(!evidence.completed_through(CalibrationStage::FitzHughNagumo));
    }
}
