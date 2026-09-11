//! Executable, outcome-blind protocol definitions.
//!
//! These helpers encode preregistered inputs and preflight checks. They do not
//! execute scientific sweeps and do not promote any empirical result.

use crate::fhn::{FhnError, FhnRun, FitzHughNagumo};
use crate::langevin::{DoubleWellLangevin, LangevinError, LangevinRun};
use crate::resonance::matched_kramers_noise_intensity;

/// Git blob SHA of `docs/research/bistable-kramers-stage0-v2.md` at protocol freeze.
pub const BISTABLE_STAGE0_V2_PROTOCOL_BLOB_SHA: &str = "54cb4bfb9880e33117df3383864cb03e7bc706e2";

/// Paired seeds frozen by the Stage 0 v2 preregistration.
pub const BISTABLE_STAGE0_V2_SEEDS: [u64; 16] = [
    101, 211, 307, 401, 503, 601, 701, 809, 907, 1009, 1103, 1201, 1301, 1409, 1511, 1601,
];

/// Frozen multiplicative grid around the independent Kramers control.
pub const BISTABLE_STAGE0_V2_GRID_FACTORS: [f64; 7] = [0.25, 0.40, 0.63, 1.00, 1.58, 2.50, 4.00];

/// Frozen forcing frequency of the separately preregistered falsification regime.
pub const BISTABLE_STAGE0_V2_FALSIFICATION_FREQUENCY_HZ: f64 = 0.20;

/// Git blob SHA of `docs/research/coherence-resonance.md` at protocol freeze.
pub const FHN_COHERENCE_STAGE0_PROTOCOL_BLOB_SHA: &str = "dafd7b7b4cf9764a96c1cdcc0c2fbd458a462445";

/// Paired deterministic seeds frozen by the FHN coherence calibration protocol.
pub const FHN_COHERENCE_STAGE0_SEEDS: [u64; 8] = [11, 23, 37, 41, 53, 67, 79, 97];

/// Frozen non-adaptive FHN noise-amplitude scan.
pub const FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES: [f64; 8] =
    [0.02, 0.03, 0.05, 0.075, 0.10, 0.20, 0.40, 0.70];

/// Frozen uncertainty multiplier used by the calibration decision heuristic.
pub const FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT: f64 = 2.0;

/// Broad preregistered acceptance interval for the sampled calibration coordinate.
pub const FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL: [f64; 2] = [0.03, 0.20];

/// Fully materialized Stage 0 v2 primary-regime inputs.
#[derive(Debug, Clone, PartialEq)]
pub struct BistableStage0V2 {
    pub model: DoubleWellLangevin,
    pub run: LangevinRun,
    pub seeds: [u64; 16],
    pub predicted_noise_intensity: f64,
    pub positive_noise_grid: [f64; 7],
}

impl BistableStage0V2 {
    /// Materialize and validate the frozen primary protocol before any trajectory is generated.
    pub fn primary() -> Result<Self, LangevinError> {
        let model = DoubleWellLangevin::new(1.2, 1.0, 0.15, 0.04, 0.0)?;
        if !model.forcing_is_subthreshold() {
            return Err(LangevinError::NonPositive("subthreshold forcing margin"));
        }

        let dt = 0.02;
        let total_periods = 80usize;
        let burn_in_periods = 20usize;
        let steps_per_period = ((1.0 / model.forcing_frequency_hz) / dt).round() as usize;
        let run = LangevinRun::new(
            dt,
            total_periods * steps_per_period,
            burn_in_periods * steps_per_period,
        )?;

        let predicted_noise_intensity = matched_kramers_noise_intensity(
            model.barrier_height(),
            model.kramers_prefactor_rate(),
            model.forcing_frequency_hz,
        )
        .map_err(LangevinError::Resonance)?
        .ok_or(LangevinError::NonPositive(
            "Kramers matched noise intensity",
        ))?;

        let positive_noise_grid =
            BISTABLE_STAGE0_V2_GRID_FACTORS.map(|factor| factor * predicted_noise_intensity);

        Ok(Self {
            model,
            run,
            seeds: BISTABLE_STAGE0_V2_SEEDS,
            predicted_noise_intensity,
            positive_noise_grid,
        })
    }

    /// Recompute the analytic Kramers control for the frozen falsification frequency.
    ///
    /// This is deliberately a preflight-only calculation. `None` means that the
    /// half-period matching equation has no positive finite solution under the
    /// model's weak-noise prefactor, so a grid derived from such a control must
    /// not be fabricated or executed.
    pub fn falsification_kramers_control() -> Result<Option<f64>, LangevinError> {
        let model = DoubleWellLangevin::new(
            1.2,
            1.0,
            0.15,
            BISTABLE_STAGE0_V2_FALSIFICATION_FREQUENCY_HZ,
            0.0,
        )?;
        matched_kramers_noise_intensity(
            model.barrier_height(),
            model.kramers_prefactor_rate(),
            model.forcing_frequency_hz,
        )
        .map_err(LangevinError::Resonance)
    }
}

/// Outcome-blind executable form of the frozen FHN coherence calibration.
#[derive(Debug, Clone, PartialEq)]
pub struct FhnCoherenceStage0 {
    pub model: FitzHughNagumo,
    pub run: FhnRun,
    pub noise_amplitudes: [f64; 8],
    pub seeds: [u64; 8],
    pub uncertainty_weight: f64,
    pub acceptance_interval: [f64; 2],
}

impl FhnCoherenceStage0 {
    /// Materialize the preregistered calibration inputs without running a sweep.
    pub fn materialize() -> Result<Self, FhnError> {
        let model = FitzHughNagumo::new(0.01, 1.05)?;
        let run = FhnRun::new(0.001, 80_000, 10_000, 0.0, 5)?;

        debug_assert!(FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
            .windows(2)
            .all(|pair| pair[0] < pair[1]));
        debug_assert!(
            FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL[0]
                <= FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL[1]
        );

        Ok(Self {
            model,
            run,
            noise_amplitudes: FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES,
            seeds: FHN_COHERENCE_STAGE0_SEEDS,
            uncertainty_weight: FHN_COHERENCE_STAGE0_UNCERTAINTY_WEIGHT,
            acceptance_interval: FHN_COHERENCE_STAGE0_ACCEPTANCE_INTERVAL,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage0_v2_materialization_matches_frozen_protocol() {
        let protocol = BistableStage0V2::primary().unwrap();
        assert_eq!(protocol.model.a, 1.2);
        assert_eq!(protocol.model.b, 1.0);
        assert_eq!(protocol.model.forcing_amplitude, 0.15);
        assert_eq!(protocol.model.forcing_frequency_hz, 0.04);
        assert_eq!(protocol.model.phase_rad, 0.0);
        assert!(protocol.model.forcing_is_subthreshold());
        assert_eq!(protocol.run.dt, 0.02);
        assert_eq!(protocol.run.steps, 100_000);
        assert_eq!(protocol.run.burn_in_steps, 25_000);
        assert_eq!(protocol.seeds, BISTABLE_STAGE0_V2_SEEDS);
    }

    #[test]
    fn stage0_v2_kramers_control_is_analytic_and_grid_is_nonadaptive() {
        let protocol = BistableStage0V2::primary().unwrap();
        let expected = 0.2958709422522067;
        assert!((protocol.predicted_noise_intensity - expected).abs() < 1e-12);
        for (actual, factor) in protocol
            .positive_noise_grid
            .iter()
            .zip(BISTABLE_STAGE0_V2_GRID_FACTORS)
        {
            assert!((*actual - factor * expected).abs() < 1e-12);
        }
    }

    #[test]
    fn stage0_v2_falsification_regime_fails_closed_without_kramers_match() {
        let control = BistableStage0V2::falsification_kramers_control().unwrap();
        assert_eq!(control, None);
    }

    #[test]
    fn fhn_stage0_materialization_matches_frozen_protocol() {
        let protocol = FhnCoherenceStage0::materialize().unwrap();
        assert_eq!(protocol.model, FitzHughNagumo::new(0.01, 1.05).unwrap());
        assert_eq!(
            protocol.run,
            FhnRun::new(0.001, 80_000, 10_000, 0.0, 5).unwrap()
        );
        assert_eq!(
            protocol.noise_amplitudes,
            FHN_COHERENCE_STAGE0_NOISE_AMPLITUDES
        );
        assert_eq!(protocol.seeds, FHN_COHERENCE_STAGE0_SEEDS);
        assert_eq!(protocol.uncertainty_weight, 2.0);
        assert_eq!(protocol.acceptance_interval, [0.03, 0.20]);
    }

    #[test]
    fn fhn_stage0_scan_is_strictly_increasing_and_acceptance_is_interior() {
        let protocol = FhnCoherenceStage0::materialize().unwrap();
        assert!(protocol
            .noise_amplitudes
            .windows(2)
            .all(|pair| pair[0].is_finite() && pair[0] >= 0.0 && pair[0] < pair[1]));
        assert!(protocol.acceptance_interval[0] > protocol.noise_amplitudes[0]);
        assert!(protocol.acceptance_interval[1] < *protocol.noise_amplitudes.last().unwrap());
    }
}
