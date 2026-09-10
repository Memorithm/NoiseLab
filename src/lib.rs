//! NoiseLab research substrate.
//!
//! NoiseLab studies reproducible stochastic perturbations and their measured
//! effects on dynamical, numerical, physical and machine-learning systems.
//! Scientific claims must be tied to executable experiments; a numerical
//! optimum is not automatically a universal law or a proof.

#![forbid(unsafe_code)]

pub mod scirust_bridge;

pub use scirust_bridge::{
    gaussian_white_noise, ornstein_uhlenbeck_path, spectral_signature, GaussianNoise,
    NoiseInputError, ScirustSpectralSignature,
};
