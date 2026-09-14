//! Stage U2 panel runner (outcome-blind execution path).
//!
//! Execution order is fixed:
//!
//! 1. [`U2Readiness::preregistered().validate()`] **before** any pair outcome;
//! 2. validate the actual panel configuration against the frozen contract;
//! 3. generate the four frozen residual series and materialize the manifest;
//! 4. analyze all six unordered pairs under both null families;
//! 5. emit the complete matrix, retaining protocol failures.
//!
//! Smoke mode (`surrogates_per_null = 19`) is explicitly **non-scientific**. It
//! exercises the readiness → manifest → dual-null → classify path only. The
//! scientific default remains `199` surrogates per null family.

use crate::u2_manifest::materialize_u2_manifest;
use crate::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_PAIRS};
use crate::u2_readiness::{U2Readiness, U2ReadinessError, U2_MIN_SURROGATES_PER_NULL};
use crate::u2_sources::{
    generate_u2_residuals, residual_for_family, U2ResidualProvenance, U2ResidualSeries,
    U2SourceError, U2_DATA_SEED_ROOT,
};
use crate::universality_u2::{
    analyze_u2_pair, StageU2AnalysisError, StageU2PairAnalysis, StageU2PairRequest,
};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Minimum surrogate repetitions accepted for a labeled non-scientific smoke run.
pub const U2_SMOKE_SURROGATES_PER_NULL: usize = 19;
/// Science seed root for surrogate-job derivation (distinct from U1).
pub const U2_SURROGATE_SEED_ROOT: u64 = 0x5532_554e_4956_324c;

/// Whether this panel invocation uses the preregistered scientific workload.
///
/// A scientific workload is not by itself evidence of universality. The complete
/// matrix, including failures, still needs provenance and independent review.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageU2PanelMode {
    /// Preregistered scientific load: ≥ 199 surrogates per null per pair.
    Scientific,
    /// Explicitly non-scientific CI / developer smoke (19 surrogates).
    NonScientificSmoke,
}

impl StageU2PanelMode {
    #[must_use]
    pub const fn surrogates_per_null(self) -> usize {
        match self {
            Self::Scientific => U2_MIN_SURROGATES_PER_NULL,
            Self::NonScientificSmoke => U2_SMOKE_SURROGATES_PER_NULL,
        }
    }

    #[must_use]
    pub const fn is_scientific(self) -> bool {
        matches!(self, Self::Scientific)
    }
}

/// Configuration frozen before inspecting Stage U2 outcomes.
///
/// Fields remain public for compatibility, but the canonical panel rejects
/// changes to scales, alpha or the implementation's recorded seed roots before
/// generating any source. Custom experiments need a separately declared protocol
/// and must not be labeled as executions of this frozen panel.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU2PanelConfig {
    pub mode: StageU2PanelMode,
    pub scales: Vec<usize>,
    pub alpha: f64,
    pub data_seed: u64,
    pub surrogate_seed_root: u64,
}

impl Default for StageU2PanelConfig {
    fn default() -> Self {
        Self::scientific()
    }
}

impl StageU2PanelConfig {
    /// Preregistered scientific panel (199 surrogates / null / pair).
    #[must_use]
    pub fn scientific() -> Self {
        Self {
            mode: StageU2PanelMode::Scientific,
            scales: vec![1, 2, 4, 8, 16],
            alpha: 0.05,
            data_seed: U2_DATA_SEED_ROOT,
            surrogate_seed_root: U2_SURROGATE_SEED_ROOT,
        }
    }

    /// Explicitly labeled non-scientific smoke configuration.
    #[must_use]
    pub fn non_scientific_smoke() -> Self {
        Self {
            mode: StageU2PanelMode::NonScientificSmoke,
            scales: vec![1, 2, 4, 8, 16],
            alpha: 0.05,
            data_seed: U2_DATA_SEED_ROOT,
            surrogate_seed_root: U2_SURROGATE_SEED_ROOT,
        }
    }
}

/// Complete Stage U2 panel matrix.
#[derive(Debug, Clone, PartialEq)]
pub struct StageU2PanelResult {
    pub mode: StageU2PanelMode,
    /// The frozen scientific workload was selected, not a positive-result gate.
    pub scientific_claim_permitted: bool,
    pub residual_provenance: Vec<U2ResidualProvenance>,
    pub pairs: Vec<StageU2PairAnalysis>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StageU2PanelError {
    Readiness(U2ReadinessError),
    Source(U2SourceError),
    Analysis(StageU2AnalysisError),
    /// The input sink failed before any pair analysis was started.
    Capture(String),
    InvalidAlpha(f64),
    InvalidSmokeSurrogateCount(usize),
    /// The named field differs from the existing frozen panel configuration.
    PreregistrationDrift(&'static str),
}

impl Display for StageU2PanelError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Readiness(error) => write!(f, "U2 readiness failed: {error}"),
            Self::Source(error) => write!(f, "U2 source generation failed: {error}"),
            Self::Analysis(error) => write!(f, "U2 pair analysis failed: {error}"),
            Self::Capture(error) => write!(f, "U2 input capture failed: {error}"),
            Self::InvalidAlpha(alpha) => {
                write!(f, "alpha must be finite and in (0, 1], got {alpha}")
            }
            Self::InvalidSmokeSurrogateCount(count) => write!(
                f,
                "non-scientific smoke requires exactly {U2_SMOKE_SURROGATES_PER_NULL} surrogates, got {count}"
            ),
            Self::PreregistrationDrift(field) => {
                write!(f, "U2 configuration differs from the frozen protocol: {field}")
            }
        }
    }
}

impl Error for StageU2PanelError {}

impl From<U2ReadinessError> for StageU2PanelError {
    fn from(value: U2ReadinessError) -> Self {
        Self::Readiness(value)
    }
}

impl From<U2SourceError> for StageU2PanelError {
    fn from(value: U2SourceError) -> Self {
        Self::Source(value)
    }
}

impl From<StageU2AnalysisError> for StageU2PanelError {
    fn from(value: StageU2AnalysisError) -> Self {
        Self::Analysis(value)
    }
}

/// Execute the Stage U2 panel.
///
/// Readiness and the actual configuration are validated before generating any
/// residuals. Smoke runs never flip `scientific_claim_permitted` to true.
pub fn run_stage_u2_panel(
    config: &StageU2PanelConfig,
) -> Result<StageU2PanelResult, StageU2PanelError> {
    run_stage_u2_panel_with_capture(config, |_, _| Ok(()))
}

/// Capture the exact generated residuals and job manifest before pair analysis.
///
/// The sink receives immutable references to the same arrays subsequently used
/// by both null families. Sources are not regenerated. Readiness/configuration
/// validation precedes the sink, which is called once after source generation.
/// An error from the sink stops execution before any pair result is produced.
/// Capturing inputs does not authorize or validate a scientific claim.
///
/// # Examples
///
/// ```
/// use noiselab::universality_u2_panel::{
///     run_stage_u2_panel_with_capture, StageU2PanelConfig, StageU2PanelError,
/// };
/// let mut config = StageU2PanelConfig::non_scientific_smoke();
/// config.data_seed ^= 1;
/// let result = run_stage_u2_panel_with_capture(&config, |_, _| {
///     panic!("invalid configuration must not reach the sink")
/// });
/// assert!(matches!(result, Err(StageU2PanelError::PreregistrationDrift(_))));
/// ```
pub fn run_stage_u2_panel_with_capture<F>(
    config: &StageU2PanelConfig,
    capture: F,
) -> Result<StageU2PanelResult, StageU2PanelError>
where
    F: FnOnce(&[U2ResidualSeries; 4], &[U2SurrogateJob]) -> Result<(), String>,
{
    // Gate 1 — outcome-blind readiness before any pair outcome exists.
    U2Readiness::preregistered().validate()?;
    validate_panel_config(config)?;

    let residuals = generate_u2_residuals(config.data_seed)?;
    let plan = execution_plan(config)?;
    let jobs = materialize_u2_manifest(plan);
    capture(&residuals, &jobs).map_err(StageU2PanelError::Capture)?;

    let mut pairs = Vec::with_capacity(U2_FROZEN_PAIRS.len());
    for (pair_index, pair) in U2_FROZEN_PAIRS.iter().copied().enumerate() {
        let left = residual_for_family(&residuals, pair.left);
        let right = residual_for_family(&residuals, pair.right);
        let analysis = analyze_u2_pair(StageU2PairRequest {
            pair_index,
            pair,
            series_a: &left.residual,
            series_b: &right.residual,
            scales: &config.scales,
            alpha: config.alpha,
            jobs: &jobs,
            retain_protocol_failures: true,
        })?;
        pairs.push(analysis);
    }

    Ok(StageU2PanelResult {
        mode: config.mode,
        scientific_claim_permitted: config.mode.is_scientific(),
        residual_provenance: residuals
            .iter()
            .map(|series: &U2ResidualSeries| series.provenance.clone())
            .collect(),
        pairs,
    })
}

fn execution_plan(config: &StageU2PanelConfig) -> Result<U2ExecutionPlan, StageU2PanelError> {
    match config.mode {
        StageU2PanelMode::Scientific => {
            U2ExecutionPlan::preregistered(config.surrogate_seed_root).map_err(Into::into)
        }
        StageU2PanelMode::NonScientificSmoke => {
            if config.mode.surrogates_per_null() != U2_SMOKE_SURROGATES_PER_NULL {
                return Err(StageU2PanelError::InvalidSmokeSurrogateCount(
                    config.mode.surrogates_per_null(),
                ));
            }
            Ok(U2ExecutionPlan {
                pairs: &U2_FROZEN_PAIRS,
                surrogates_per_null: U2_SMOKE_SURROGATES_PER_NULL,
                seed_root: config.surrogate_seed_root,
            })
        }
    }
}

fn validate_panel_config(config: &StageU2PanelConfig) -> Result<(), StageU2PanelError> {
    if !config.alpha.is_finite() || config.alpha <= 0.0 || config.alpha > 1.0 {
        return Err(StageU2PanelError::InvalidAlpha(config.alpha));
    }
    // Validate supplied values, not just a newly constructed readiness object.
    // These are the original protocol values and original implementation seeds;
    // no outcome is consulted and neither scientific nor smoke defaults change.
    if config.scales != [1, 2, 4, 8, 16] {
        return Err(StageU2PanelError::PreregistrationDrift("scales"));
    }
    if config.alpha.to_bits() != 0.05_f64.to_bits() {
        return Err(StageU2PanelError::PreregistrationDrift("alpha"));
    }
    if config.data_seed != U2_DATA_SEED_ROOT {
        return Err(StageU2PanelError::PreregistrationDrift("data_seed"));
    }
    if config.surrogate_seed_root != U2_SURROGATE_SEED_ROOT {
        return Err(StageU2PanelError::PreregistrationDrift(
            "surrogate_seed_root",
        ));
    }
    if config.mode == StageU2PanelMode::NonScientificSmoke
        && config.mode.surrogates_per_null() != U2_SMOKE_SURROGATES_PER_NULL
    {
        return Err(StageU2PanelError::InvalidSmokeSurrogateCount(
            config.mode.surrogates_per_null(),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::u2_decision::StageU2Decision;

    #[test]
    fn smoke_panel_is_deterministic_and_complete() {
        let config = StageU2PanelConfig::non_scientific_smoke();
        let first = run_stage_u2_panel(&config).unwrap();
        let second = run_stage_u2_panel(&config).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.pairs.len(), 6);
        assert!(!first.scientific_claim_permitted);
        assert_eq!(first.residual_provenance.len(), 4);
        assert!(first.pairs.iter().all(|pair| {
            pair.protocol_error.is_some()
                || (pair.decision.is_some()
                    && pair.p_shuffle.is_some()
                    && pair.p_phase.is_some()
                    && pair.convergence_score.is_finite())
        }));
    }

    #[test]
    fn invalid_alpha_prevents_pair_outcomes() {
        let mut config = StageU2PanelConfig::non_scientific_smoke();
        config.alpha = 0.0;
        assert!(matches!(
            run_stage_u2_panel(&config),
            Err(StageU2PanelError::InvalidAlpha(_))
        ));
    }

    #[test]
    fn canonical_modes_validate_without_generating_outcomes() {
        for config in [
            StageU2PanelConfig::scientific(),
            StageU2PanelConfig::non_scientific_smoke(),
        ] {
            assert_eq!(validate_panel_config(&config), Ok(()));
            assert_eq!(config.scales, [1, 2, 4, 8, 16]);
            assert_eq!(config.alpha.to_bits(), 0.05_f64.to_bits());
        }
    }

    #[test]
    fn finite_but_retuned_alpha_is_rejected_before_execution() {
        for mode in [
            StageU2PanelMode::Scientific,
            StageU2PanelMode::NonScientificSmoke,
        ] {
            for alpha in [0.01, 0.1, 1.0, f64::from_bits(0.05_f64.to_bits() + 1)] {
                let mut config = StageU2PanelConfig::scientific();
                config.mode = mode;
                config.alpha = alpha;
                assert_eq!(
                    run_stage_u2_panel(&config),
                    Err(StageU2PanelError::PreregistrationDrift("alpha"))
                );
            }
        }
    }

    #[test]
    fn changed_missing_reordered_or_duplicated_scales_are_rejected() {
        for mode in [
            StageU2PanelMode::Scientific,
            StageU2PanelMode::NonScientificSmoke,
        ] {
            for scales in [
                vec![],
                vec![1, 2, 4, 8],
                vec![1, 2, 4, 8, 32],
                vec![1, 4, 2, 8, 16],
                vec![1, 2, 4, 8, 16, 16],
            ] {
                let mut config = StageU2PanelConfig::scientific();
                config.mode = mode;
                config.scales = scales;
                assert_eq!(
                    run_stage_u2_panel(&config),
                    Err(StageU2PanelError::PreregistrationDrift("scales"))
                );
            }
        }
    }

    #[test]
    fn seed_search_cannot_masquerade_as_the_frozen_panel() {
        for mode in [
            StageU2PanelMode::Scientific,
            StageU2PanelMode::NonScientificSmoke,
        ] {
            let mut config = StageU2PanelConfig::scientific();
            config.mode = mode;
            config.data_seed ^= 1;
            assert_eq!(
                run_stage_u2_panel(&config),
                Err(StageU2PanelError::PreregistrationDrift("data_seed"))
            );
            config.data_seed = U2_DATA_SEED_ROOT;
            config.surrogate_seed_root ^= 1;
            assert_eq!(
                run_stage_u2_panel(&config),
                Err(StageU2PanelError::PreregistrationDrift(
                    "surrogate_seed_root"
                ))
            );
        }
    }

    #[test]
    fn nonfinite_or_out_of_range_alpha_stays_invalid() {
        for alpha in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.1, 0.0, 1.1] {
            let mut config = StageU2PanelConfig::scientific();
            config.alpha = alpha;
            assert!(matches!(
                run_stage_u2_panel(&config),
                Err(StageU2PanelError::InvalidAlpha(_))
            ));
        }
    }

    #[test]
    fn classify_path_is_exercised_when_analysis_succeeds() {
        let result = run_stage_u2_panel(&StageU2PanelConfig::non_scientific_smoke()).unwrap();
        let decided = result
            .pairs
            .iter()
            .filter_map(|pair| pair.decision)
            .collect::<Vec<_>>();
        // At least the successful pairs must land in the frozen label set.
        for decision in decided {
            assert!(matches!(
                decision,
                StageU2Decision::NoObservedConvergence
                    | StageU2Decision::CompatibleWithMarginalNull
                    | StageU2Decision::SpectrumExplainedCandidate
                    | StageU2Decision::CrossMechanismCandidate
            ));
        }
    }
}
