//! Stage U2 panel report (executable, outcome-blind runner).
//!
//! Default / CI: non-scientific smoke with 19 surrogates per null. Scientific
//! load (199 surrogates) requires `NOISELAB_U2_FULL=1`. Flags accept only 1/0 or
//! true/false (case-insensitive). Enabling FULL and SMOKE together is an error.
//!
//! Smoke output is not scientific evidence and must not be recorded as U2
//! results. The full-load selector does not itself validate a scientific claim.

use noiselab::report_mode::ReportMode;
use noiselab::{run_stage_u2_panel, StageU2PanelConfig, StageU2PanelMode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mode = ReportMode::from_environment("NOISELAB_U2_FULL", "NOISELAB_U2_SMOKE")?;
    let config = if mode.is_scientific() {
        StageU2PanelConfig::scientific()
    } else {
        StageU2PanelConfig::non_scientific_smoke()
    };

    match config.mode {
        StageU2PanelMode::NonScientificSmoke => {
            eprintln!(
                "NoiseLab Stage U2: NON-SCIENTIFIC smoke (19 surrogates/null). Not evidence."
            );
        }
        StageU2PanelMode::Scientific => {
            eprintln!(
                "NoiseLab Stage U2: scientific panel (199 surrogates/null). Still exploratory."
            );
        }
    }

    let result = run_stage_u2_panel(&config)?;
    println!("# report_schema\t1");
    println!("# report_kind\tU2");
    println!("# mode\t{:?}", result.mode);
    println!(
        "# scientific_claim_permitted\t{}",
        result.scientific_claim_permitted
    );
    println!(
        "# surrogates_per_null\t{}",
        config.mode.surrogates_per_null()
    );
    println!("# data_seed\t{}", config.data_seed);
    println!("# surrogate_seed_root\t{}", config.surrogate_seed_root);
    println!(
        "left\tright\tfine_distance\tterminal_distance\tconvergence\tp_shuffle\tp_phase\tdecision\tprotocol_error"
    );
    for pair in result.pairs {
        println!(
            "{:?}\t{:?}\t{:.12}\t{:.12}\t{:.12}\t{}\t{}\t{}\t{}",
            pair.pair.left,
            pair.pair.right,
            pair.fine_scale_distance,
            pair.terminal_scale_distance,
            pair.convergence_score,
            option_f64(pair.p_shuffle),
            option_f64(pair.p_phase),
            option_debug(pair.decision),
            pair.protocol_error
                .as_deref()
                .unwrap_or("")
                .replace(['\t', '\n', '\r'], " "),
        );
    }

    Ok(())
}

fn option_f64(value: Option<f64>) -> String {
    match value {
        Some(v) => format!("{v:.12}"),
        None => String::new(),
    }
}

fn option_debug<T: std::fmt::Debug>(value: Option<T>) -> String {
    match value {
        Some(v) => format!("{v:?}"),
        None => String::new(),
    }
}
