//! Replay of captured Stage U2 inputs without regenerating source trajectories.
//!
//! This module consumes the exact residual arrays and surrogate-job manifest that
//! an evidence bundle has already loaded and integrity-checked. It replays the
//! existing dual-null analysis; it does not authenticate the bundle, regenerate
//! sources, retune the protocol, or promote an exploratory U2 outcome into proof.

use crate::u2_manifest::materialize_u2_manifest;
use crate::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_PAIRS, U2_FROZEN_SOURCES};
use crate::u2_readiness::{U2Readiness, U2ReadinessError, U2_SCIRUST_REVISION};
use crate::u2_sources::{
    residual_for_family, U2ResidualSeries, U2_SOURCE_BURN_IN, U2_SOURCE_SAMPLES,
};
use crate::universality_u2::{
    analyze_u2_pair, StageU2AnalysisError, StageU2PairAnalysis, StageU2PairRequest,
};
use crate::universality_u2_panel::{StageU2PanelConfig, StageU2PanelMode};
use std::error::Error;
use std::fmt::{Display, Formatter};

const DATA_SEED_STREAM_MULTIPLIER: u64 = 0x9e37_79b9_7f4a_7c15;

/// Replay output for one already-captured U2 panel input set.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU2ReplayResult {
    /// Frozen workload mode whose configuration was validated before replay.
    pub mode: StageU2PanelMode,
    /// Six pair analyses in the canonical preregistered pair order.
    pub pairs: Vec<StageU2PairAnalysis>,
}

/// Fail-closed validation or analysis errors for captured U2 replay.
#[derive(Debug, Clone, PartialEq)]
pub enum StageU2ReplayError {
    /// The caller supplied a configuration other than the canonical configuration
    /// for the declared mode.
    ConfigurationDrift,
    /// A residual occupies the wrong canonical array slot.
    ResidualFamilyMismatch {
        index: usize,
        expected: crate::u2_plan::U2SourceFamily,
        actual: crate::u2_plan::U2SourceFamily,
    },
    /// The residual and its provenance disagree about source family.
    ProvenanceFamilyMismatch { index: usize },
    /// The capture names a SciRust revision other than the frozen U2 revision.
    ScirustRevisionMismatch { index: usize },
    /// A captured residual provenance record does not contain the per-family seed
    /// implied by the canonical panel data-seed root.
    DataSeedMismatch {
        index: usize,
        expected: u64,
        actual: u64,
    },
    /// The retained residual length or recorded sample count is not canonical.
    SampleCountMismatch {
        index: usize,
        residual: usize,
        recorded: usize,
        expected: usize,
    },
    /// The recorded burn-in count differs from the frozen source contract.
    BurnInMismatch {
        index: usize,
        actual: usize,
        expected: usize,
    },
    /// A captured residual contains a non-finite binary64 value.
    NonFiniteResidual { index: usize, sample: usize },
    /// The captured surrogate-job manifest is not byte-semantic-equivalent to the
    /// canonical outcome-blind manifest for this mode.
    JobManifestMismatch { expected: usize, actual: usize },
    /// The preregistered execution plan is not currently admissible.
    Readiness(U2ReadinessError),
    /// Existing pair analysis rejected an input or numerical operation.
    Analysis(StageU2AnalysisError),
}

impl Display for StageU2ReplayError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ConfigurationDrift => {
                f.write_str("captured U2 replay configuration differs from the canonical protocol")
            }
            Self::ResidualFamilyMismatch {
                index,
                expected,
                actual,
            } => write!(
                f,
                "captured U2 residual {index} has family {actual:?}, expected {expected:?}"
            ),
            Self::ProvenanceFamilyMismatch { index } => write!(
                f,
                "captured U2 residual {index} disagrees with its provenance family"
            ),
            Self::ScirustRevisionMismatch { index } => write!(
                f,
                "captured U2 residual {index} does not name the frozen SciRust revision"
            ),
            Self::DataSeedMismatch {
                index,
                expected,
                actual,
            } => write!(
                f,
                "captured U2 residual {index} records data seed {actual:#018x}, expected {expected:#018x}"
            ),
            Self::SampleCountMismatch {
                index,
                residual,
                recorded,
                expected,
            } => write!(
                f,
                "captured U2 residual {index} has {residual} values / {recorded} recorded, expected {expected}"
            ),
            Self::BurnInMismatch {
                index,
                actual,
                expected,
            } => write!(
                f,
                "captured U2 residual {index} records burn-in {actual}, expected {expected}"
            ),
            Self::NonFiniteResidual { index, sample } => write!(
                f,
                "captured U2 residual {index} contains a non-finite value at sample {sample}"
            ),
            Self::JobManifestMismatch { expected, actual } => write!(
                f,
                "captured U2 job manifest has {actual} entries or altered contents/order; canonical manifest has {expected}"
            ),
            Self::Readiness(error) => write!(f, "U2 replay readiness failed: {error}"),
            Self::Analysis(error) => write!(f, "U2 replay analysis failed: {error}"),
        }
    }
}

impl Error for StageU2ReplayError {}

impl From<U2ReadinessError> for StageU2ReplayError {
    fn from(value: U2ReadinessError) -> Self {
        Self::Readiness(value)
    }
}

impl From<StageU2AnalysisError> for StageU2ReplayError {
    fn from(value: StageU2AnalysisError) -> Self {
        Self::Analysis(value)
    }
}

/// Replay the frozen U2 pair analyses from already-captured inputs.
///
/// `residuals` and `jobs` are expected to come from an independently
/// integrity-checked capture bundle. This function validates the canonical mode,
/// residual structure/provenance boundary and exact outcome-blind job manifest,
/// then invokes the same [`analyze_u2_pair`] implementation used by the live
/// panel. It deliberately does not regenerate the four source trajectories.
///
/// A successful replay establishes deterministic analysis agreement for the
/// supplied inputs. It is not authentication, independent physical evidence, a
/// novelty result, or proof of a common noise origin.
///
/// # Errors
///
/// Returns an error for protocol/configuration drift, malformed residual arrays,
/// altered job manifests, readiness failure, or pair-analysis failure.
pub fn replay_stage_u2_captured_inputs(
    config: &StageU2PanelConfig,
    residuals: &[U2ResidualSeries; 4],
    jobs: &[U2SurrogateJob],
) -> Result<StageU2ReplayResult, StageU2ReplayError> {
    let expected_config = match config.mode {
        StageU2PanelMode::Scientific => StageU2PanelConfig::scientific(),
        StageU2PanelMode::NonScientificSmoke => StageU2PanelConfig::non_scientific_smoke(),
    };
    if config != &expected_config {
        return Err(StageU2ReplayError::ConfigurationDrift);
    }

    U2Readiness::preregistered().validate()?;
    validate_residuals(residuals, config.data_seed)?;

    let plan = match config.mode {
        StageU2PanelMode::Scientific => U2ExecutionPlan::preregistered(config.surrogate_seed_root)?,
        StageU2PanelMode::NonScientificSmoke => U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: config.mode.surrogates_per_null(),
            seed_root: config.surrogate_seed_root,
        },
    };
    let expected_jobs = materialize_u2_manifest(plan);
    if jobs != expected_jobs {
        return Err(StageU2ReplayError::JobManifestMismatch {
            expected: expected_jobs.len(),
            actual: jobs.len(),
        });
    }

    let mut pairs = Vec::with_capacity(U2_FROZEN_PAIRS.len());
    for (pair_index, pair) in U2_FROZEN_PAIRS.iter().copied().enumerate() {
        let left = residual_for_family(residuals, pair.left);
        let right = residual_for_family(residuals, pair.right);
        pairs.push(analyze_u2_pair(StageU2PairRequest {
            pair_index,
            pair,
            series_a: &left.residual,
            series_b: &right.residual,
            scales: &config.scales,
            alpha: config.alpha,
            jobs,
            retain_protocol_failures: true,
        })?);
    }

    Ok(StageU2ReplayResult {
        mode: config.mode,
        pairs,
    })
}

fn validate_residuals(
    residuals: &[U2ResidualSeries; 4],
    data_seed_root: u64,
) -> Result<(), StageU2ReplayError> {
    for (index, expected) in U2_FROZEN_SOURCES.iter().copied().enumerate() {
        let series = &residuals[index];
        if series.family != expected {
            return Err(StageU2ReplayError::ResidualFamilyMismatch {
                index,
                expected,
                actual: series.family,
            });
        }
        if series.provenance.family != series.family {
            return Err(StageU2ReplayError::ProvenanceFamilyMismatch { index });
        }
        if series.provenance.scirust_revision != U2_SCIRUST_REVISION {
            return Err(StageU2ReplayError::ScirustRevisionMismatch { index });
        }
        let expected_seed = derive_source_seed(data_seed_root, index);
        if series.provenance.data_seed != expected_seed {
            return Err(StageU2ReplayError::DataSeedMismatch {
                index,
                expected: expected_seed,
                actual: series.provenance.data_seed,
            });
        }
        if series.residual.len() != U2_SOURCE_SAMPLES
            || series.provenance.samples != U2_SOURCE_SAMPLES
        {
            return Err(StageU2ReplayError::SampleCountMismatch {
                index,
                residual: series.residual.len(),
                recorded: series.provenance.samples,
                expected: U2_SOURCE_SAMPLES,
            });
        }
        if series.provenance.burn_in_samples != U2_SOURCE_BURN_IN {
            return Err(StageU2ReplayError::BurnInMismatch {
                index,
                actual: series.provenance.burn_in_samples,
                expected: U2_SOURCE_BURN_IN,
            });
        }
        if let Some(sample) = series.residual.iter().position(|value| !value.is_finite()) {
            return Err(StageU2ReplayError::NonFiniteResidual { index, sample });
        }
    }
    Ok(())
}

fn derive_source_seed(root: u64, family_index: usize) -> u64 {
    let stream = u64::try_from(family_index + 1).expect("four U2 families fit in u64");
    root ^ stream.wrapping_mul(DATA_SEED_STREAM_MULTIPLIER)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::u2_sources::generate_u2_residuals;
    use crate::universality_u2_panel::run_stage_u2_panel;

    fn smoke_inputs() -> (
        StageU2PanelConfig,
        [U2ResidualSeries; 4],
        Vec<U2SurrogateJob>,
    ) {
        let config = StageU2PanelConfig::non_scientific_smoke();
        let residuals = generate_u2_residuals(config.data_seed).unwrap();
        let plan = U2ExecutionPlan {
            pairs: &U2_FROZEN_PAIRS,
            surrogates_per_null: config.mode.surrogates_per_null(),
            seed_root: config.surrogate_seed_root,
        };
        let jobs = materialize_u2_manifest(plan);
        (config, residuals, jobs)
    }

    #[test]
    fn captured_replay_matches_live_smoke_analysis_exactly() {
        let (config, residuals, jobs) = smoke_inputs();
        let live = run_stage_u2_panel(&config).unwrap();
        let replay = replay_stage_u2_captured_inputs(&config, &residuals, &jobs).unwrap();
        assert_eq!(replay.mode, live.mode);
        assert_eq!(replay.pairs, live.pairs);
    }

    #[test]
    fn altered_job_order_is_rejected_before_analysis() {
        let (config, residuals, mut jobs) = smoke_inputs();
        jobs.swap(0, 1);
        assert_eq!(
            replay_stage_u2_captured_inputs(&config, &residuals, &jobs),
            Err(StageU2ReplayError::JobManifestMismatch {
                expected: jobs.len(),
                actual: jobs.len(),
            })
        );
    }

    #[test]
    fn non_finite_captured_value_is_rejected() {
        let (config, mut residuals, jobs) = smoke_inputs();
        residuals[2].residual[17] = f64::NAN;
        assert_eq!(
            replay_stage_u2_captured_inputs(&config, &residuals, &jobs),
            Err(StageU2ReplayError::NonFiniteResidual {
                index: 2,
                sample: 17,
            })
        );
    }

    #[test]
    fn noncanonical_data_seed_is_rejected() {
        let (config, _, jobs) = smoke_inputs();
        let residuals = generate_u2_residuals(config.data_seed ^ 1).unwrap();
        let expected = derive_source_seed(config.data_seed, 0);
        assert_eq!(
            replay_stage_u2_captured_inputs(&config, &residuals, &jobs),
            Err(StageU2ReplayError::DataSeedMismatch {
                index: 0,
                expected,
                actual: residuals[0].provenance.data_seed,
            })
        );
    }

    #[test]
    fn retuned_config_is_rejected_without_replaying() {
        let (mut config, residuals, jobs) = smoke_inputs();
        config.alpha = 0.1;
        assert_eq!(
            replay_stage_u2_captured_inputs(&config, &residuals, &jobs),
            Err(StageU2ReplayError::ConfigurationDrift)
        );
    }
}
