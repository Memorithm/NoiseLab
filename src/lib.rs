//! NoiseLab research substrate.
//!
//! NoiseLab studies reproducible stochastic perturbations and their measured
//! effects on dynamical, numerical, physical and machine-learning systems.
//! Scientific claims must be tied to executable experiments; a numerical
//! optimum is not automatically a universal law or a proof.

#![forbid(unsafe_code)]

pub mod calibration;
pub mod evidence;
pub mod langevin;
pub mod resonance;
pub mod scirust_bridge;

pub use calibration::{
    measure_steady_state_displacement_amplitude, sweep_linear_resonance, CalibrationError,
    DrivenLinearOscillator, LinearResonanceCalibration, ResponseMeasurement,
};
pub use evidence::{detect_edge_separated_peak, EdgeSeparatedPeak, EvidenceError};
pub use langevin::{
    calibrate_kramers_matching, coherent_switching_response, simulate_double_well,
    CoherentResponse, DoubleWellLangevin, KramersCalibration, LangevinError, LangevinRun,
    LangevinTrajectory, NoiseResponse,
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
