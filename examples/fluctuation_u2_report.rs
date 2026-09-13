//! Stage U2 panel report (executable, outcome-blind runner).
//!
//! ## Smoke vs full
//!
//! - Default / CI: **non-scientific smoke** with 19 surrogates per null.
//!   Set `NOISELAB_U2_SMOKE=1` explicitly (CI does this).
//! - Scientific load (199 surrogates): set `NOISELAB_U2_FULL=1`.
//!
//! Smoke output is **not** scientific evidence and must not be written into
//! `docs/research/fluctuation-universality-u2-results.md`.

use noiselab::{run_stage_u2_panel, StageU2PanelConfig, StageU2PanelMode};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let full = std::env::var("NOISELAB_U2_FULL")
        .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    let smoke_env = std::env::var("NOISELAB_U2_SMOKE")
        .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
        .unwrap_or(false);

    let config = if full {
        StageU2PanelConfig::scientific()
    } else {
        // Default is smoke; NOISELAB_U2_SMOKE=1 documents the CI intent.
        let _ = smoke_env;
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
    println!(
        "mode\tscientific_claim_permitted\t{}",
        result.scientific_claim_permitted
    );
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
                .replace(['\t', '\n'], " "),
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
