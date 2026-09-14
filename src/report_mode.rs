//! Explicit workload selection for NoiseLab report executables.
//!
//! Unset flags select non-scientific smoke. Invalid or contradictory flags fail
//! before any experiment runs. This selector does not authorize a scientific
//! claim, change a protocol, or grant access to a protected holdout.
// Copyright 2026 Tarek Zekriti
// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0

use std::ffi::OsStr;
use std::fmt::{Display, Formatter};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportMode {
    NonScientificSmoke,
    Scientific,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportModeError {
    ConflictingFlags,
    InvalidFlag(&'static str),
    InvalidEnvironmentKeys,
}

impl Display for ReportModeError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ConflictingFlags => write!(f, "FULL and SMOKE cannot both be enabled"),
            Self::InvalidFlag(flag) => {
                write!(f, "invalid {flag} flag: expected 1/0 or true/false")
            }
            Self::InvalidEnvironmentKeys => {
                write!(f, "FULL and SMOKE require distinct valid environment keys")
            }
        }
    }
}

impl std::error::Error for ReportModeError {}

impl ReportMode {
    /// Read flags without treating non-Unicode values as missing variables.
    pub fn from_environment(full_key: &str, smoke_key: &str) -> Result<Self, ReportModeError> {
        if full_key == smoke_key
            || [full_key, smoke_key]
                .iter()
                .any(|key| key.is_empty() || key.contains(['=', '\0']))
        {
            return Err(ReportModeError::InvalidEnvironmentKeys);
        }
        let full = std::env::var_os(full_key);
        let smoke = std::env::var_os(smoke_key);
        Self::from_flags(full.as_deref(), smoke.as_deref())
    }

    /// Pure selector, testable without changing the process environment.
    pub fn from_flags(
        full: Option<&OsStr>,
        smoke: Option<&OsStr>,
    ) -> Result<Self, ReportModeError> {
        let full = parse_flag(full, "FULL")?;
        let smoke = parse_flag(smoke, "SMOKE")?;
        if full && smoke {
            return Err(ReportModeError::ConflictingFlags);
        }
        Ok(if full {
            Self::Scientific
        } else {
            Self::NonScientificSmoke
        })
    }

    #[must_use]
    pub const fn is_scientific(self) -> bool {
        matches!(self, Self::Scientific)
    }
}

fn parse_flag(value: Option<&OsStr>, name: &'static str) -> Result<bool, ReportModeError> {
    let Some(value) = value else {
        return Ok(false);
    };
    let value = value.to_str().ok_or(ReportModeError::InvalidFlag(name))?;
    if value == "1" || value.eq_ignore_ascii_case("true") {
        Ok(true)
    } else if value == "0" || value.eq_ignore_ascii_case("false") {
        Ok(false)
    } else {
        Err(ReportModeError::InvalidFlag(name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_supported_flag_combinations_have_explicit_semantics() {
        let values = [
            (None, false),
            (Some("0"), false),
            (Some("false"), false),
            (Some("FaLsE"), false),
            (Some("1"), true),
            (Some("true"), true),
            (Some("TRUE"), true),
        ];
        for (full_value, full) in values {
            for (smoke_value, smoke) in values {
                let actual =
                    ReportMode::from_flags(full_value.map(OsStr::new), smoke_value.map(OsStr::new));
                let expected = if full && smoke {
                    Err(ReportModeError::ConflictingFlags)
                } else if full {
                    Ok(ReportMode::Scientific)
                } else {
                    Ok(ReportMode::NonScientificSmoke)
                };
                assert_eq!(actual, expected);
            }
        }
        assert!(ReportMode::Scientific.is_scientific());
        assert!(!ReportMode::NonScientificSmoke.is_scientific());
    }

    #[test]
    fn invalid_flags_never_silently_select_a_workload() {
        for value in ["", "yes", "2", " true", "false ", "\n1"] {
            assert_eq!(
                ReportMode::from_flags(Some(OsStr::new(value)), None),
                Err(ReportModeError::InvalidFlag("FULL"))
            );
            assert_eq!(
                ReportMode::from_flags(Some(OsStr::new("1")), Some(OsStr::new(value))),
                Err(ReportModeError::InvalidFlag("SMOKE"))
            );
        }
    }

    #[test]
    fn invalid_environment_keys_fail_before_reading_the_environment() {
        let cases = [
            ("", "SMOKE"),
            ("A=B", "SMOKE"),
            ("FULL", "A\0B"),
            ("A", "A"),
        ];
        for (full, smoke) in cases {
            assert_eq!(
                ReportMode::from_environment(full, smoke),
                Err(ReportModeError::InvalidEnvironmentKeys)
            );
        }
    }

    #[cfg(unix)]
    #[test]
    fn non_unicode_values_are_invalid_not_absent() {
        use std::os::unix::ffi::OsStrExt;
        let invalid = OsStr::from_bytes(&[0xff]);
        assert_eq!(
            ReportMode::from_flags(Some(invalid), None),
            Err(ReportModeError::InvalidFlag("FULL"))
        );
        assert_eq!(
            ReportMode::from_flags(None, Some(invalid)),
            Err(ReportModeError::InvalidFlag("SMOKE"))
        );
    }
}
