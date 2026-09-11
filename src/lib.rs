//! NoiseLab research substrate.
//!
//! NoiseLab studies reproducible stochastic perturbations and their measured
//! effects on dynamical, numerical, physical and machine-learning systems.
//! Scientific claims must be tied to executable experiments; a numerical
//! optimum is not automatically a universal law or a proof.

#![forbid(unsafe_code)]

pub mod attention;
pub mod calibration;
pub mod evidence;
pub mod fhn;
pub mod fhn_stage0;
pub mod information;
pub mod langevin;
pub mod laser_calibration;
pub mod preregistered;
pub mod resonance;
pub mod scirust_bridge;
pub mod universality;

pub use attention::{
    evaluate_flat_rope_gaussian_perturbation, flat_rope_control, AttentionExperimentError,
    AttentionNoiseSite, AttentionPerturbationResponse, AttentionPerturbationSpec,
};
pub use calibration::{
    measure_steady_state_displacement_amplitude, sweep_linear_resonance, CalibrationError,
    DrivenLinearOscillator, LinearResonanceCalibration, ResponseMeasurement,
};
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
pub use universality::{
    multiscale_trace, test_fluctuation_universality, FluctuationSignature, MultiscaleTrace,
    ScaleDistance, UniversalityError, UniversalityEvidence, UniversalityTestConfig,
    UniversalityTestResult,
};
