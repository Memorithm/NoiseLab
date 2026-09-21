//! Bistable Langevin / Kramers Stage 0 v2 report (executable, outcome-blind).
//!
//! Default / CI: non-scientific smoke with abbreviated periods
//! (`4` total / `1` burn-in) and the first `2` preregistered seeds. Scientific
//! load requires `NOISELAB_BISTABLE_FULL=1` and uses the preregistered budget
//! (`80` / `20` periods, all `16` seeds). Flags accept only 1/0 or true/false
//! (case-insensitive). Enabling FULL and SMOKE together is an error.
//!
//! Smoke output is not scientific evidence. A full run still requires a retained
//! producing commit SHA and review before any Kramers-compatibility or novelty
//! claim. No results.md is emitted by this runner.
//! Protocol: `docs/research/bistable-kramers-stage0-v2.md`.

use noiselab::report_mode::ReportMode;
use noiselab::{
    run_bistable_kramers_stage0, BistableStage0Config, BistableStage0Decision,
    BistableStage0DeterministicControl, BistableStage0FalsificationStatus, BistableStage0Mode,
    BistableStage0Summary,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mode = ReportMode::from_environment("NOISELAB_BISTABLE_FULL", "NOISELAB_BISTABLE_SMOKE")?;
    let config = if mode.is_scientific() {
        BistableStage0Config::scientific()
    } else {
        BistableStage0Config::non_scientific_smoke()
    };

    match config.mode {
        BistableStage0Mode::NonScientificSmoke => {
            eprintln!(
                "NoiseLab bistable Kramers Stage 0 v2: NON-SCIENTIFIC smoke (abbreviated periods/seeds). Not evidence."
            );
        }
        BistableStage0Mode::Scientific => {
            eprintln!(
                "NoiseLab bistable Kramers Stage 0 v2: scientific panel (preregistered budget). Still exploratory until a results record with commit SHA."
            );
        }
    }

    let result = run_bistable_kramers_stage0(&config)?;
    println!("# report_schema\t1");
    println!("# report_kind\tBISTABLE_KRAMERS_STAGE0");
    println!("# mode\t{:?}", result.mode);
    println!(
        "# scientific_claim_permitted\t{}",
        result.scientific_claim_permitted
    );
    println!("# summary\t{}", format_summary(&result.summary));
    println!(
        "# primary_decision\t{}",
        format_decision(&result.primary_decision)
    );
    println!("# total_periods\t{}", result.total_periods);
    println!("# burn_in_periods\t{}", result.burn_in_periods);
    println!("# seed_count\t{}", result.seed_count);
    println!(
        "# predicted_noise_intensity\t{:.12}",
        result.analytic.predicted_noise_intensity
    );
    println!("# barrier_height\t{:.12}", result.analytic.barrier_height);
    println!(
        "# kramers_prefactor_rate\t{:.12}",
        result.analytic.kramers_prefactor_rate
    );
    println!("# protocol_blob_sha\t{}", result.analytic.protocol_blob_sha);
    println!(
        "# forcing_is_subthreshold\t{}",
        result.analytic.forcing_is_subthreshold
    );
    println!(
        "# deterministic_control\t{}",
        format_deterministic(&result.deterministic_control)
    );
    println!(
        "# falsification\t{}",
        format_falsification(&result.falsification)
    );
    println!(
        "regime\tnoise_intensity\tmean_amplitude\tsample_stddev\tstandard_error\treplicates\tper_seed_amplitudes"
    );
    for point in &result.primary_grid {
        println!(
            "primary\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{}\t{}",
            point.noise_intensity,
            point.mean_coherent_amplitude,
            point.sample_stddev,
            point.standard_error,
            point.replicates,
            join_seeds(&point.per_seed_amplitudes),
        );
    }
    for point in &result.zero_forcing_grid {
        println!(
            "zero_forcing\t{:.12}\t{:.12}\t{:.12}\t{:.12}\t{}\t{}",
            point.noise_intensity,
            point.mean_coherent_amplitude,
            point.sample_stddev,
            point.standard_error,
            point.replicates,
            join_seeds(&point.per_seed_amplitudes),
        );
    }

    Ok(())
}

fn join_seeds(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:.12}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn format_summary(summary: &BistableStage0Summary) -> String {
    match summary {
        BistableStage0Summary::Classified(decision) => {
            format!("Classified({})", format_decision(decision))
        }
        other => format!("{other:?}"),
    }
}

fn format_decision(decision: &BistableStage0Decision) -> String {
    match decision {
        BistableStage0Decision::BestSampledPoint {
            noise_intensity,
            mean_coherent_amplitude,
            grid_index,
        } => format!(
            "BestSampledPoint(D={noise_intensity:.12},amp={mean_coherent_amplitude:.12},idx={grid_index})"
        ),
        BistableStage0Decision::H0NotRejected { best_sampled } => match best_sampled {
            Some(best) => format!(
                "H0NotRejected(best_D={:.12},best_amp={:.12},idx={})",
                best.noise_intensity, best.mean_coherent_amplitude, best.grid_index
            ),
            None => "H0NotRejected".to_string(),
        },
        BistableStage0Decision::InteriorPeakWithoutMechanismCompatibility {
            peak,
            predicted_noise_intensity,
            best_sampled,
        } => format!(
            "InteriorPeakWithoutMechanismCompatibility(peak_D={:.12},pred_D={:.12},best_idx={})",
            peak.coordinate, predicted_noise_intensity, best_sampled.grid_index
        ),
        BistableStage0Decision::MechanismCompatibleInteriorPeak {
            peak,
            predicted_noise_intensity,
            best_sampled,
        } => format!(
            "MechanismCompatibleInteriorPeak(peak_D={:.12},pred_D={:.12},best_idx={})",
            peak.coordinate, predicted_noise_intensity, best_sampled.grid_index
        ),
        BistableStage0Decision::ControlFailure { detail } => format!(
            "ControlFailure({})",
            detail.replace(['\t', '\n', '\r'], " ")
        ),
        BistableStage0Decision::ProtocolMismatch { detail } => format!(
            "ProtocolMismatch({})",
            detail.replace(['\t', '\n', '\r'], " ")
        ),
    }
}

fn format_deterministic(control: &BistableStage0DeterministicControl) -> String {
    match control {
        BistableStage0DeterministicControl::SeedIndependent { amplitude } => {
            format!("SeedIndependent(amp={amplitude:.12})")
        }
        BistableStage0DeterministicControl::SeedDependent => "SeedDependent".to_string(),
        BistableStage0DeterministicControl::Failed { detail } => {
            format!("Failed({})", detail.replace(['\t', '\n', '\r'], " "))
        }
    }
}

fn format_falsification(status: &BistableStage0FalsificationStatus) -> String {
    match status {
        BistableStage0FalsificationStatus::KramersControlInapplicable => {
            "KramersControlInapplicable".to_string()
        }
        BistableStage0FalsificationStatus::Executed {
            predicted_noise_intensity,
            decision,
        } => format!(
            "Executed(pred_D={:.12},decision={})",
            predicted_noise_intensity,
            format_decision(decision)
        ),
        BistableStage0FalsificationStatus::Failed { detail } => {
            format!("Failed({})", detail.replace(['\t', '\n', '\r'], " "))
        }
    }
}
