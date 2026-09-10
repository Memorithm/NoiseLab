//! NoiseLab research substrate.
//!
//! NoiseLab studies reproducible stochastic perturbations and their measured
//! effects on dynamical, numerical, physical and machine-learning systems.
//! Scientific claims must be tied to executable experiments; a numerical
//! optimum is not automatically a universal law or a proof.

#![forbid(unsafe_code)]

pub mod calibration;
pub mod resonance;
pub mod scirust_bridge;

pub use calibration::{
    CalibrationError, DrivenLinearOscillator, LinearResonanceCalibration, ResponseMeasurement,
    measure_steady_state_displacement_amplitude, sweep_linear_resonance,
};
pub use resonance::{
    CandidateEstimate, InteriorPeak, OperatingPoint, ResonanceInputError, SearchConfig, SearchError,
    SearchResult, SweepSample, detect_interior_response_peak, kramers_escape_rate,
    linear_displacement_resonance_angular_frequency, matched_kramers_noise_intensity,
    search_operating_points, stochastic_resonance_target_rate,
};
pub use scirust_bridge::{
    GaussianNoise, NoiseInputError, ScirustSpectralSignature, gaussian_white_noise,
    ornstein_uhlenbeck_path, spectral_signature,
};
