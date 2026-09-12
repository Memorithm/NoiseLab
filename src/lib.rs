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
pub mod fhn_stage0;
pub mod information;
pub mod langevin;
pub mod laser_calibration;
pub mod preregistered;
pub mod protocol_provenance;
pub mod resonance;
pub mod scirust_bridge;
pub mod spectral_null;
pub mod u2_decision;
pub mod u2_manifest;
pub mod u2_plan;
pub mod u2_readiness;
pub mod universality;
pub mod universality_panel;

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
pub use fhn_stage0::{classify_fhn_stage0, FhnStage0Decision};
pub use information::{
    audit_noise_information, discrete_entropy_bits, histogram_mutual_information_bits,
    InformationError, NoiseInformationAudit,
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
pub use universality::{
    multiscale_trace, test_fluctuation_universality, FluctuationSignature, MultiscaleTrace,
    ScaleDistance, UniversalityError, UniversalityEvidence, UniversalityTestConfig,
    UniversalityTestResult,
};
pub use universality_panel::{
    run_stage_u1_panel, StageU1Comparison, StageU1Error, StageU1PairPurpose, StageU1PanelConfig,
    StageU1PanelResult, StageU1Source,
};
