//! FHN observation-horizon Stage 0 report (executable, outcome-blind runner).
//!
//! ## Smoke vs full
//!
//! - Default / CI: **non-scientific smoke** with abbreviated horizons
//!   `[2_000, 4_000, 8_000]` (burn-in `500`). Set `NOISELAB_FHN_HORIZON_SMOKE=1`
//!   explicitly if desired; smoke is already the default when `FULL` is unset.
//! - Scientific load (preregistered `[40_000, 80_000, 160_000]`, burn-in
//!   `10_000`, frozen noise grid / seeds / acceptance): set
//!   `NOISELAB_FHN_HORIZON_FULL=1`.
//!
//! Smoke output is **not** scientific evidence and must not be written into a
//! results markdown file. A full run still requires documenting the producing
//! commit SHA before any robustness claim; this example alone does not authorize
//! a novelty claim.
//!
//! Protocol: `docs/research/fhn-observation-horizon-stage0.md`.

use noiselab::{
    run_fhn_horizon_stage0, FhnHorizonStage0Config, FhnHorizonStage0HorizonDecision,
    FhnHorizonStage0Mode, FhnHorizonStage0Summary,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let full = std::env::var("NOISELAB_FHN_HORIZON_FULL")
        .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    let smoke_env = std::env::var("NOISELAB_FHN_HORIZON_SMOKE")
        .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
        .unwrap_or(false);

    let config = if full {
        FhnHorizonStage0Config::scientific()
    } else {
        // Default is smoke; NOISELAB_FHN_HORIZON_SMOKE=1 documents CI intent.
        let _ = smoke_env;
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
    println!(
        "mode\tscientific_claim_permitted\tsummary\t{}",
        format_summary(&result.summary)
    );
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
            detail: detail.replace(['\t', '\n'], " "),
        },
    }
}
