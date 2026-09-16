//! Exact Stage U2 pair-result capture bound to retained score replay.
//!
//! This support module writes the post-analysis pair fields as exact binary64 bit
//! patterns only after re-running the frozen score-only replay and confirming that
//! the retained p-values and decision match the analysis object exactly. It is an
//! integrity/reproducibility surface, not a scientific acceptance path.

use noiselab::u2_plan::{U2ExecutionPlan, U2_FROZEN_PAIRS};
use noiselab::{replay_u2_statistics_from_scores, StageU2PairAnalysis, StageU2PanelConfig};
use std::fs::{File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Persist exact post-analysis U2 pair fields after verifying score replay.
///
/// `SURROGATE_SCORES_COMPLETE` must already exist in `destination`. Successful
/// pair rows are replayed from the in-memory score evidence using the frozen
/// `+1` p-value rule and decision function. The capture rejects any bit-level
/// p-value drift or decision mismatch before publishing `PAIR_RESULTS_COMPLETE`.
/// Protocol-failure rows are preserved without fabricating p-values or decisions;
/// their error text must be non-empty so a retained row cannot become ambiguous.
///
/// # Errors
///
/// Returns `InvalidData` when pair identity/order drifts, score evidence has not
/// been retained first, a successful row contains non-finite or incomplete
/// values, replay disagrees with the analysis fields, or a protocol-failure row
/// has empty error text or contains a partial statistical result.
pub fn capture_pair_results(
    destination: &Path,
    config: &StageU2PanelConfig,
    pairs: &[StageU2PairAnalysis],
) -> io::Result<()> {
    if !destination.join("SURROGATE_SCORES_COMPLETE").is_file() {
        return Err(invalid_data(
            "U2 pair-result capture requires completed surrogate-score retention",
        ));
    }
    if pairs.len() != U2_FROZEN_PAIRS.len() {
        return Err(invalid_data(
            "U2 pair-result capture does not contain all frozen pairs",
        ));
    }

    let plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };

    let mut rows = new_file(destination, "pair_results.tsv")?;
    writeln!(
        rows,
        "pair_index\tleft\tright\tfine_distance_bits\tterminal_distance_bits\tconvergence_score_bits\tp_shuffle_bits\tp_phase_bits\tdecision\tprotocol_error"
    )?;

    let mut successful = 0usize;
    let mut failed = 0usize;
    for (expected_index, pair) in pairs.iter().enumerate() {
        if pair.pair_index != expected_index || pair.pair != U2_FROZEN_PAIRS[expected_index] {
            return Err(invalid_data("U2 pair-result identity mismatch"));
        }

        let (p_shuffle_bits, p_phase_bits, decision) = if let Some(error) = &pair.protocol_error {
            if error.is_empty() {
                return Err(invalid_data("protocol-failed U2 pair has empty error text"));
            }
            if pair.p_shuffle.is_some()
                || pair.p_phase.is_some()
                || pair.decision.is_some()
                || !pair.surrogate_scores.is_empty()
            {
                return Err(invalid_data(
                    "protocol-failed U2 pair contains partial statistical evidence",
                ));
            }
            failed = failed
                .checked_add(1)
                .ok_or(invalid_data("U2 pair-result failure count overflow"))?;
            (String::new(), String::new(), String::new())
        } else {
            if !pair.fine_scale_distance.is_finite()
                || !pair.terminal_scale_distance.is_finite()
                || !pair.convergence_score.is_finite()
            {
                return Err(invalid_data(
                    "successful U2 pair contains non-finite observed statistics",
                ));
            }
            let p_shuffle = pair
                .p_shuffle
                .ok_or(invalid_data("successful U2 pair is missing p_shuffle"))?;
            let p_phase = pair
                .p_phase
                .ok_or(invalid_data("successful U2 pair is missing p_phase"))?;
            let pair_decision = pair
                .decision
                .ok_or(invalid_data("successful U2 pair is missing a decision"))?;
            if !p_shuffle.is_finite() || !p_phase.is_finite() {
                return Err(invalid_data(
                    "successful U2 pair contains non-finite p-values",
                ));
            }

            let replay = replay_u2_statistics_from_scores(
                plan,
                expected_index,
                pair.convergence_score,
                config.alpha,
                &pair.surrogate_scores,
            )
            .map_err(|_| invalid_data("U2 pair-result score replay failed"))?;
            if replay.p_shuffle.to_bits() != p_shuffle.to_bits()
                || replay.p_phase.to_bits() != p_phase.to_bits()
                || replay.decision != pair_decision
            {
                return Err(invalid_data(
                    "U2 pair-result statistics do not match frozen score replay",
                ));
            }

            successful = successful
                .checked_add(1)
                .ok_or(invalid_data("U2 pair-result success count overflow"))?;
            (
                format!("{:016x}", p_shuffle.to_bits()),
                format!("{:016x}", p_phase.to_bits()),
                format!("{pair_decision:?}"),
            )
        };

        writeln!(
            rows,
            "{}\t{:?}\t{:?}\t{:016x}\t{:016x}\t{:016x}\t{}\t{}\t{}\t{}",
            pair.pair_index,
            pair.pair.left,
            pair.pair.right,
            pair.fine_scale_distance.to_bits(),
            pair.terminal_scale_distance.to_bits(),
            pair.convergence_score.to_bits(),
            p_shuffle_bits,
            p_phase_bits,
            decision,
            clean(pair.protocol_error.as_deref().unwrap_or("")),
        )?;
    }
    finish(rows)?;

    let mut marker = new_file(destination, "PAIR_RESULTS_COMPLETE")?;
    writeln!(marker, "pair_result_export_complete=true")?;
    writeln!(marker, "rows={}", pairs.len())?;
    writeln!(marker, "successful_rows={successful}")?;
    writeln!(marker, "protocol_error_rows={failed}")?;
    writeln!(marker, "alpha_bits={:016x}", config.alpha.to_bits())?;
    writeln!(marker, "surrogate_scores_persisted=true")?;
    writeln!(marker, "scientific_evidence=false")?;
    finish(marker)
}

fn new_file(directory: &Path, name: &str) -> io::Result<BufWriter<File>> {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(directory.join(name))
        .map(BufWriter::new)
}

fn finish(mut writer: BufWriter<File>) -> io::Result<()> {
    writer.flush()?;
    writer.get_ref().sync_all()
}

fn clean(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}
