//! NoiseLab research substrate.
//!
//! NoiseLab studies reproducible stochastic perturbations and their measured
//! effects on dynamical, numerical, physical and machine-learning systems.
//! Scientific claims must be tied to executable experiments; a numerical
//! optimum is not automatically a universal law or a proof.

#![forbid(unsafe_code)]

pub mod attention;
pub mod calibration;
pub mod calibration_evidence;
pub mod calibration_ladder;
pub mod evidence;
pub mod fhn;
pub mod fhn_horizon;
pub mod fhn_stage0;
pub mod information;
pub mod lagged_information;
pub mod langevin;
pub mod laser_calibration;
pub mod preregistered;
pub mod protocol_provenance;
pub mod report_mode;
pub mod resonance;
pub mod scirust_bridge;
pub mod spectral_null;
pub mod u2_decision;
pub mod u2_manifest;
pub mod u2_plan;
pub mod u2_readiness;
pub mod u2_replay;
pub mod u2_score_replay;
pub mod u2_sources;
pub mod u2_surrogate_realization;
pub mod universality;
pub mod universality_panel;
pub mod universality_u2;
pub mod universality_u2_panel;

pub use attention::{
    evaluate_flat_rope_gaussian_perturbation, flat_rope_control, AttentionExperimentError,
    AttentionNoiseSite, AttentionPerturbationResponse, AttentionPerturbationSpec,
};
pub use calibration::{
    measure_steady_state_displacement_amplitude, sweep_linear_resonance, CalibrationError,
    DrivenLinearOscillator, LinearResonanceCalibration, ResponseMeasurement,
};
pub use calibration_evidence::{
    evidence_from_receipts, CalibrationReceipt, CalibrationReceiptError,
};
pub use calibration_ladder::{CalibrationEvidence, CalibrationGateError, CalibrationStage};
pub use evidence::{detect_edge_separated_peak, EdgeSeparatedPeak, EvidenceError};
pub use fhn::{
    calibrate_coherence_resonance, inter_spike_stats, simulate_fhn_spikes, CoherenceCalibration,
    CoherenceMinimum, CoherenceResponse, FhnError, FhnRun, FhnSpikeTrain, FitzHughNagumo,
    InterSpikeStats,
};
pub use fhn_horizon::{
    classify_fhn_horizon_robustness, run_fhn_horizon_stage0, FhnHorizonRobustness,
    FhnHorizonStage0Config, FhnHorizonStage0Error, FhnHorizonStage0HorizonDecision,
    FhnHorizonStage0Mode, FhnHorizonStage0Observation, FhnHorizonStage0Result,
    FhnHorizonStage0Summary, FHN_HORIZON_SMOKE_BURN_IN_STEPS, FHN_HORIZON_SMOKE_STEPS,
    FHN_HORIZON_STAGE0_BURN_IN_STEPS, FHN_HORIZON_STAGE0_STEPS,
};
pub use fhn_stage0::{classify_fhn_stage0, FhnStage0Decision};
pub use information::{
    audit_noise_information, cyclic_shift_mutual_information_null, discrete_entropy_bits,
    histogram_mutual_information_bits, permutation_mutual_information_null,
    spectral_surrogate_mutual_information_null, InformationError, MutualInformationCyclicShiftNull,
    MutualInformationPermutationNull, MutualInformationSpectralNull, NoiseInformationAudit,
    MAX_INFORMATION_CYCLIC_SHIFTS, MAX_INFORMATION_PERMUTATIONS,
    MAX_INFORMATION_SPECTRAL_SURROGATES,
};
pub use lagged_information::{
    lagged_histogram_mutual_information_bits, LaggedInformationError, LaggedMutualInformationPoint,
    LaggedMutualInformationScan,
};
pub use langevin::{
    calibrate_kramers_matching, coherent_switching_response, simulate_double_well,
    CoherentResponse, DoubleWellLangevin, KramersCalibration, LangevinError, LangevinRun,
    LangevinTrajectory, NoiseResponse,
};
pub use laser_calibration::{
    calibrate_laser_relaxation, LaserCalibrationError, LaserRelaxationCalibration,
    LaserRingMeasurement,
};
pub use preregistered::{
    BistableStage0V2, FhnCoherenceStage0, BISTABLE_STAGE0_V2_GRID_FACTORS,
    BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA, BISTABLE_STAGE0_V2_SEEDS,
    FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL, FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES,
    FHN_COHERENCE_STAGE0_PROTOCOL_BLOB_SHA, FHN_COHERENCE_STAGE0_SEEDS,
    FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
};
pub use protocol_provenance::{
    validate_protocol_provenance, GitIdentityKind, ProtocolProvenance, ProtocolProvenanceError,
    PROTOCOL_PROVENANCE,
};
pub use resonance::{
    detect_interior_response_peak, kramers_escape_rate,
    linear_displacement_resonance_angular_frequency, matched_kramers_noise_intensity,
    search_operating_points, stochastic_resonance_target_rate, CandidateEstimate, InteriorPeak,
    OperatingPoint, ResonanceInputError, SearchConfig, SearchError, SearchResult, SweepSample,
};
pub use scirust_bridge::{
    gaussian_white_noise, ornstein_uhlenbeck_path, spectral_signature, GaussianNoise,
    NoiseInputError, ScirustSpectralSignature,
};
pub use spectral_null::spectral_phase_null;
pub use u2_decision::{classify_stage_u2, StageU2Decision, StageU2DecisionError};
pub use u2_manifest::materialize_u2_manifest;
pub use u2_plan::{
    U2ExecutionPlan, U2NullFamily, U2Pair, U2SourceFamily, U2_FROZEN_PAIRS, U2_FROZEN_SOURCES,
};
pub use u2_readiness::{
    U2Readiness, U2ReadinessError, U2_MIN_SURROGATES_PER_NULL, U2_SCIRUST_REVISION,
    U2_SOURCE_FAMILIES, U2_UNORDERED_PAIRS,
};
pub use u2_replay::{replay_stage_u2_captured_inputs, StageU2ReplayError, StageU2ReplayResult};
pub use u2_score_replay::{
    replay_u2_statistics_from_scores, ReplayedU2Statistics, U2ScoreReplayError,
};
pub use u2_surrogate_realization::{
    materialize_u2_surrogate_pair, U2SurrogateRealization, U2SurrogateRealizationError,
};
pub use universality::{
    multiscale_trace, observed_multiscale_comparison, test_fluctuation_universality,
    FluctuationSignature, MultiscaleTrace, ObservedMultiscaleComparison, ScaleDistance,
    UniversalityError, UniversalityEvidence, UniversalityTestConfig, UniversalityTestResult,
};
pub use universality_panel::{
    run_stage_u1_panel, StageU1Comparison, StageU1Error, StageU1PairPurpose, StageU1PanelConfig,
    StageU1PanelResult, StageU1Source,
};

pub use u2_sources::{
    generate_u2_residuals, residual_for_family, U2ResidualProvenance, U2ResidualSeries,
    U2SourceError, U2_DATA_SEED_ROOT, U2_SOURCE_BURN_IN, U2_SOURCE_SAMPLES,
};
pub use universality_u2::{
    analyze_u2_pair, StageU2AnalysisError, StageU2PairAnalysis, StageU2PairRequest,
    U2SurrogateScore, SPECTRAL_RIGHT_SEED_TAG,
};
pub use universality_u2_panel::{
    run_stage_u2_panel, StageU2PanelConfig, StageU2PanelError, StageU2PanelMode,
    StageU2PanelResult, U2_SMOKE_SURROGATES_PER_NULL, U2_SURROGATE_SEED_ROOT,
};
