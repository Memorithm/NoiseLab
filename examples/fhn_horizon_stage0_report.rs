//! FHN observation-horizon Stage 0 report (executable, outcome-blind runner).
//!
//! Default / CI: non-scientific smoke with abbreviated horizons
//! `[2_000, 4_000, 8_000]` (burn-in `500`). Scientific load requires
//! `NOISELAB_FHN_HORIZON_FULL=1` and uses the preregistered horizons
//! `[40_000, 80_000, 160_000]` (burn-in `10_000`). Flags accept only 1/0 or
//! true/false (case-insensitive). Enabling FULL and SMOKE together is an error.
//!
//! Smoke output is not scientific evidence. A full run still requires a retained
//! producing commit SHA and review before any robustness or novelty claim.
//! Protocol: `docs/research/fhn-observation-horizon-stage0.md`.

use noiselab::report_mode::ReportMode;
use noiselab::{
    run_fhn_horizon_stage0, FhnHorizonStage0Config, FhnHorizonStage0HorizonDecision,
    FhnHorizonStage0Mode, FhnHorizonStage0Summary,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mode =
        ReportMode::from_environment("NOISELAB_FHN_HORIZON_FULL", "NOISELAB_FHN_HORIZON_SMOKE")?;
    let config = if mode.is_scientific() {
        FhnHorizonStage0Config::scientific()
    } else {
        FhnHorizonStage0Config::non_scientific_smoke()
    };

    match config.mode {
        FhnHorizonStage0Mode::NonScientificSmoke => {
            eprintln!(
                "NoiseLab FHN horizon Stage 0: NON-SCIENTIFIC smoke (abbreviated horizons). Not evidence."
            );
        }
        FhnHorizonStage0Mode::Scientific => {
            eprintln!(
                "NoiseLab FHN horizon Stage 0: scientific panel (preregistered horizons). Still exploratory until a results record with commit SHA."
            );
        }
    }

    let result = run_fhn_horizon_stage0(&config)?;
    println!("# report_schema\t1");
    println!("# report_kind\tFHN_HORIZON_STAGE0");
    println!("# mode\t{:?}", result.mode);
    println!(
        "# scientific_claim_permitted\t{}",
        result.scientific_claim_permitted
    );
    println!("# summary\t{}", format_summary(&result.summary));
    println!("# burn_in_steps\t{}", config.mode.burn_in_steps());
    println!("steps\tdecision\tnoise_amplitude\tseed\tobserved_spikes\trequired_spikes\tdetail");
    for observation in &result.observations {
        let row = format_decision(&observation.decision);
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            observation.steps,
            row.decision,
            row.noise_amplitude,
            row.seed,
            row.observed,
            row.required,
            row.detail,
        );
    }

    Ok(())
}

struct DecisionRow {
    decision: String,
    noise_amplitude: String,
    seed: String,
    observed: String,
    required: String,
    detail: String,
}

fn format_summary(summary: &FhnHorizonStage0Summary) -> String {
    match summary {
        FhnHorizonStage0Summary::Classified(robustness) => format!("Classified({robustness:?})"),
        other => format!("{other:?}"),
    }
}

fn format_decision(decision: &FhnHorizonStage0HorizonDecision) -> DecisionRow {
    match decision {
        FhnHorizonStage0HorizonDecision::Stage0(stage0) => DecisionRow {
            decision: format!("Stage0({stage0:?})"),
            noise_amplitude: match stage0 {
                noiselab::FhnStage0Decision::Accepted { noise_amplitude }
                | noiselab::FhnStage0Decision::OutsideAcceptance { noise_amplitude } => {
                    format!("{noise_amplitude:.12}")
                }
                _ => String::new(),
            },
            seed: String::new(),
            observed: String::new(),
            required: String::new(),
            detail: String::new(),
        },
        FhnHorizonStage0HorizonDecision::InsufficientSpikes {
            noise_amplitude,
            seed,
            observed,
            required,
        } => DecisionRow {
            decision: "InsufficientSpikes".to_string(),
            noise_amplitude: format!("{noise_amplitude:.12}"),
            seed: seed.to_string(),
            observed: observed.to_string(),
            required: required.to_string(),
            detail: String::new(),
        },
        FhnHorizonStage0HorizonDecision::ProtocolFailure { detail } => DecisionRow {
            decision: "ProtocolFailure".to_string(),
            noise_amplitude: String::new(),
            seed: String::new(),
            observed: String::new(),
            required: String::new(),
            detail: detail.replace(['\t', '\n', '\r'], " "),
        },
    }
}
