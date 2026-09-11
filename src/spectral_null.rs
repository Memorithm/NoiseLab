//! Spectral null-model adapter for preregistered universality controls.
//!
//! This module consumes SciRust's phase-randomized surrogate primitive at the
//! immutable revision pinned in `Cargo.toml`. It does not execute Stage U2 or
//! interpret a surrogate as evidence for universality; it only exposes the
//! deterministic null transformation needed by that later experiment.

use scirust_signal::surrogate::{phase_randomized_surrogate, SurrogateError};

/// Construct a deterministic phase-randomized spectral null.
///
/// SciRust preserves the input Fourier magnitudes (up to floating-point
/// roundoff), DC and Nyquist components while randomizing the remaining phase
/// relationships from an explicit seed. The upstream validation contract is
/// deliberately retained rather than duplicated here.
pub fn spectral_phase_null(signal: &[f64], seed: u64) -> Result<Vec<f64>, SurrogateError> {
    phase_randomized_surrogate(signal, seed)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_is_reproducible_through_noiselab_boundary() {
        let signal: Vec<f64> = (0..64)
            .map(|index| {
                let x = index as f64;
                (0.17 * x).sin() + 0.25 * (0.41 * x).cos()
            })
            .collect();

        let a = spectral_phase_null(&signal, 0x5eed).unwrap();
        let b = spectral_phase_null(&signal, 0x5eed).unwrap();
        assert_eq!(
            a.iter().map(|value| value.to_bits()).collect::<Vec<_>>(),
            b.iter().map(|value| value.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn malformed_input_is_rejected_by_upstream_contract() {
        assert!(spectral_phase_null(&[0.0; 6], 1).is_err());
    }

    #[test]
    fn wrapper_does_not_silently_modify_length() {
        let signal: Vec<f64> = (0..128).map(|index| (0.09 * index as f64).sin()).collect();
        let surrogate = spectral_phase_null(&signal, 9).unwrap();
        assert_eq!(surrogate.len(), signal.len());
    }
}
