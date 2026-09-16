//! Typed replay of retained Stage U2 pair-result artifacts.
//!
//! This verifier parses `pair_results.tsv` and its completion marker, validates
//! the frozen six-pair identity/order and exact binary64 encodings, and reapplies
//! the frozen decision function. A stronger replay path also re-reads the
//! retained `surrogate_scores.tsv` archive and requires every successful pair
//! result to reproduce the frozen p-values and decision directly from those
//! on-disk score rows.

use noiselab::u2_decision::{classify_stage_u2, StageU2Decision};
use noiselab::u2_plan::{U2ExecutionPlan, U2NullFamily, U2_FROZEN_PAIRS};
use noiselab::{replay_u2_statistics_from_scores, StageU2PanelConfig, U2SurrogateScore};
use std::fs;
use std::io;
use std::path::Path;

/// Typed representation of one retained Stage U2 pair result.
#[derive(Debug, Clone, PartialEq)]
pub struct ArchivedU2PairResult {
    pub pair_index: usize,
    pub fine_scale_distance: f64,
    pub terminal_scale_distance: f64,
    pub convergence_score: f64,
    pub p_shuffle: Option<f64>,
    pub p_phase: Option<f64>,
    pub decision: Option<StageU2Decision>,
    pub protocol_error: Option<String>,
}

/// Parse and validate the retained six-pair Stage U2 result artifact.
///
/// Successful rows must contain finite exact-binary64 statistics and a decision
/// that is reproduced by the frozen classifier under the supplied alpha.
/// Protocol-failure rows must contain no p-values or decision. The completion
/// marker must agree with the six parsed rows and the exact configured alpha.
///
/// This is an integrity/reproducibility check only. It does not authenticate the
/// bundle, re-read surrogate-score evidence, open a protected holdout or grant a
/// scientific verdict.
///
/// # Errors
///
/// Returns `InvalidData` for malformed schema/rows, pair-order drift, non-finite
/// successful values, decision drift, completion-marker mismatch or unexpected
/// extra/missing rows.
pub fn verify_archived_pair_results(
    destination: &Path,
    config: &StageU2PanelConfig,
) -> io::Result<Vec<ArchivedU2PairResult>> {
    let marker = fs::read_to_string(destination.join("PAIR_RESULTS_COMPLETE"))?;
    validate_marker(&marker, config)?;

    let contents = fs::read_to_string(destination.join("pair_results.tsv"))?;
    let mut lines = contents.lines();
    if lines.next()
        != Some(
            "pair_index\tleft\tright\tfine_distance_bits\tterminal_distance_bits\tconvergence_score_bits\tp_shuffle_bits\tp_phase_bits\tdecision\tprotocol_error",
        )
    {
        return Err(invalid_data("invalid U2 pair-result header"));
    }

    let mut rows = Vec::new();
    rows.try_reserve_exact(U2_FROZEN_PAIRS.len())
        .map_err(|_| invalid_data("unable to allocate U2 pair-result replay rows"))?;
    for (pair_index, expected_pair) in U2_FROZEN_PAIRS.iter().copied().enumerate() {
        let line = lines
            .next()
            .ok_or(invalid_data("incomplete U2 pair-result artifact"))?;
        rows.push(parse_row(line, pair_index, expected_pair, config.alpha)?);
    }
    if lines.next().is_some() {
        return Err(invalid_data("unexpected extra U2 pair-result row"));
    }

    let successful = rows
        .iter()
        .filter(|row| row.protocol_error.is_none())
        .count();
    let failed = rows.len() - successful;
    if marker_usize(&marker, "successful_rows")? != successful
        || marker_usize(&marker, "protocol_error_rows")? != failed
    {
        return Err(invalid_data("U2 pair-result completion counts mismatch"));
    }
    Ok(rows)
}

/// Re-read retained pair results and bind every successful row to the retained
/// surrogate-score archive on disk.
///
/// This first performs [`verify_archived_pair_results`], then parses
/// `surrogate_scores.tsv` under the frozen execution plan. Successful pairs must
/// have the complete canonical dual-null score sequence; protocol-failure pairs
/// must have no retained scores. The frozen score-only replay must reproduce the
/// archived pair-result p-values bit-for-bit and the archived decision exactly.
///
/// This is stronger disk-to-disk integrity evidence than capture-time binding,
/// but it still does not authenticate wholesale bundle replacement, revalidate
/// how surrogate scores were produced, open a protected holdout or grant a
/// scientific verdict.
///
/// # Errors
///
/// Returns `InvalidData` for malformed or reordered score archives, incomplete
/// successful-pair score sets, retained scores for a protocol-failure pair, or
/// any p-value/decision disagreement between the score archive and pair results.
pub fn verify_archived_pair_results_against_scores(
    destination: &Path,
    config: &StageU2PanelConfig,
) -> io::Result<Vec<ArchivedU2PairResult>> {
    let rows = verify_archived_pair_results(destination, config)?;
    let plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    let scores = read_archived_scores(destination, plan)?;

    for row in &rows {
        let pair_scores = &scores[row.pair_index];
        if row.protocol_error.is_some() {
            if !pair_scores.is_empty() {
                return Err(invalid_data(
                    "protocol-failed U2 pair retains surrogate scores",
                ));
            }
            continue;
        }

        let replay = replay_u2_statistics_from_scores(
            plan,
            row.pair_index,
            row.convergence_score,
            config.alpha,
            pair_scores,
        )
        .map_err(|_| invalid_data("archived U2 score-to-pair replay failed"))?;

        if row.p_shuffle.map(f64::to_bits) != Some(replay.p_shuffle.to_bits())
            || row.p_phase.map(f64::to_bits) != Some(replay.p_phase.to_bits())
            || row.decision != Some(replay.decision)
        {
            return Err(invalid_data(
                "archived U2 pair result does not match retained scores",
            ));
        }
    }

    Ok(rows)
}

fn read_archived_scores(
    destination: &Path,
    plan: U2ExecutionPlan,
) -> io::Result<Vec<Vec<U2SurrogateScore>>> {
    let marker = fs::read_to_string(destination.join("SURROGATE_SCORES_COMPLETE"))?;
    if marker_value(&marker, "surrogate_score_export_complete") != Some("true")
        || marker_value(&marker, "scientific_evidence") != Some("false")
    {
        return Err(invalid_data("invalid U2 surrogate-score completion marker"));
    }

    let contents = fs::read_to_string(destination.join("surrogate_scores.tsv"))?;
    let mut lines = contents.lines();
    if lines.next()
        != Some("pair_index\tleft\tright\tnull_family\trepetition\tseed\tconvergence_score_bits")
    {
        return Err(invalid_data("invalid U2 surrogate-score header"));
    }

    let mut grouped = (0..U2_FROZEN_PAIRS.len())
        .map(|_| Vec::new())
        .collect::<Vec<Vec<U2SurrogateScore>>>();
    let mut total_rows = 0usize;
    let mut previous_pair_index = None;
    for line in lines {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 7 {
            return Err(invalid_data("invalid U2 surrogate-score row"));
        }

        let pair_index = parse_usize(fields[0], "invalid U2 score pair index")?;
        if previous_pair_index.is_some_and(|previous| pair_index < previous) {
            return Err(invalid_data("reordered U2 surrogate-score pair groups"));
        }
        previous_pair_index = Some(pair_index);

        let null_family = parse_null_family(fields[3])?;
        let repetition = parse_usize(fields[4], "invalid U2 score repetition")?;
        let seed = parse_u64(fields[5], "invalid U2 score seed")?;
        let score_bits = u64::from_str_radix(fields[6], 16)
            .map_err(|_| invalid_data("invalid U2 convergence-score bits"))?;
        let convergence_score = f64::from_bits(score_bits);
        if !convergence_score.is_finite() {
            return Err(invalid_data("non-finite archived U2 convergence score"));
        }

        let expected_job = plan
            .surrogate_job(pair_index, null_family, repetition)
            .ok_or(invalid_data("U2 archived score references an invalid job"))?;
        if fields[1] != format!("{:?}", expected_job.pair.left)
            || fields[2] != format!("{:?}", expected_job.pair.right)
            || seed != expected_job.seed
        {
            return Err(invalid_data("U2 archived score identity mismatch"));
        }

        grouped[pair_index].push(U2SurrogateScore {
            job: expected_job,
            convergence_score,
        });
        total_rows = total_rows
            .checked_add(1)
            .ok_or(invalid_data("U2 surrogate-score row-count overflow"))?;
    }

    if marker_usize(&marker, "rows")? != total_rows {
        return Err(invalid_data("U2 surrogate-score completion count mismatch"));
    }

    let expected_per_pair = plan
        .surrogates_per_null
        .checked_mul(2)
        .ok_or(invalid_data("U2 score replay count overflow"))?;
    if grouped
        .iter()
        .any(|pair_scores| !pair_scores.is_empty() && pair_scores.len() != expected_per_pair)
    {
        return Err(invalid_data("partial successful-pair U2 score archive"));
    }

    Ok(grouped)
}

fn validate_marker(marker: &str, config: &StageU2PanelConfig) -> io::Result<()> {
    if marker_value(marker, "pair_result_export_complete") != Some("true")
        || marker_usize(marker, "rows")? != U2_FROZEN_PAIRS.len()
        || marker_value(marker, "surrogate_scores_persisted") != Some("true")
        || marker_value(marker, "scientific_evidence") != Some("false")
    {
        return Err(invalid_data("invalid U2 pair-result completion marker"));
    }
    let alpha_bits = marker_value(marker, "alpha_bits")
        .ok_or(invalid_data("missing U2 pair-result alpha bits"))?;
    let parsed_alpha = u64::from_str_radix(alpha_bits, 16)
        .map_err(|_| invalid_data("invalid U2 pair-result alpha bits"))?;
    if parsed_alpha != config.alpha.to_bits() {
        return Err(invalid_data("U2 pair-result alpha mismatch"));
    }
    Ok(())
}

fn parse_row(
    line: &str,
    expected_index: usize,
    expected_pair: noiselab::u2_plan::U2Pair,
    alpha: f64,
) -> io::Result<ArchivedU2PairResult> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 10 {
        return Err(invalid_data("invalid U2 pair-result row"));
    }
    let pair_index = parse_usize(fields[0], "invalid U2 pair-result pair index")?;
    if pair_index != expected_index
        || fields[1] != format!("{:?}", expected_pair.left)
        || fields[2] != format!("{:?}", expected_pair.right)
    {
        return Err(invalid_data("U2 pair-result identity mismatch"));
    }

    let fine_scale_distance = parse_f64_bits(fields[3], "invalid U2 fine-distance bits")?;
    let terminal_scale_distance = parse_f64_bits(fields[4], "invalid U2 terminal-distance bits")?;
    let convergence_score = parse_f64_bits(fields[5], "invalid U2 convergence-score bits")?;
    let protocol_error = (!fields[9].is_empty()).then(|| fields[9].to_owned());

    if protocol_error.is_some() {
        if !fields[6].is_empty() || !fields[7].is_empty() || !fields[8].is_empty() {
            return Err(invalid_data(
                "protocol-failed U2 pair-result row contains statistical fields",
            ));
        }
        return Ok(ArchivedU2PairResult {
            pair_index,
            fine_scale_distance,
            terminal_scale_distance,
            convergence_score,
            p_shuffle: None,
            p_phase: None,
            decision: None,
            protocol_error,
        });
    }

    if !fine_scale_distance.is_finite()
        || !terminal_scale_distance.is_finite()
        || !convergence_score.is_finite()
    {
        return Err(invalid_data(
            "successful U2 pair-result row contains non-finite observed statistics",
        ));
    }
    let p_shuffle = parse_finite_f64_bits(fields[6], "invalid U2 p_shuffle bits")?;
    let p_phase = parse_finite_f64_bits(fields[7], "invalid U2 p_phase bits")?;
    let decision = parse_decision(fields[8])?;
    let replayed = classify_stage_u2(convergence_score, p_shuffle, p_phase, alpha)
        .map_err(|_| invalid_data("invalid U2 pair-result decision inputs"))?;
    if replayed != decision {
        return Err(invalid_data("U2 pair-result decision mismatch"));
    }

    Ok(ArchivedU2PairResult {
        pair_index,
        fine_scale_distance,
        terminal_scale_distance,
        convergence_score,
        p_shuffle: Some(p_shuffle),
        p_phase: Some(p_phase),
        decision: Some(decision),
        protocol_error: None,
    })
}

fn parse_null_family(value: &str) -> io::Result<U2NullFamily> {
    match value {
        "ShuffledMarginal" => Ok(U2NullFamily::ShuffledMarginal),
        "PhaseRandomizedSpectrum" => Ok(U2NullFamily::PhaseRandomizedSpectrum),
        _ => Err(invalid_data("invalid U2 archived null family")),
    }
}

fn parse_decision(value: &str) -> io::Result<StageU2Decision> {
    match value {
        "NoObservedConvergence" => Ok(StageU2Decision::NoObservedConvergence),
        "CompatibleWithMarginalNull" => Ok(StageU2Decision::CompatibleWithMarginalNull),
        "SpectrumExplainedCandidate" => Ok(StageU2Decision::SpectrumExplainedCandidate),
        "CrossMechanismCandidate" => Ok(StageU2Decision::CrossMechanismCandidate),
        _ => Err(invalid_data("invalid U2 archived decision")),
    }
}

fn parse_finite_f64_bits(value: &str, message: &'static str) -> io::Result<f64> {
    let parsed = parse_f64_bits(value, message)?;
    if !parsed.is_finite() {
        return Err(invalid_data(message));
    }
    Ok(parsed)
}

fn parse_f64_bits(value: &str, message: &'static str) -> io::Result<f64> {
    let bits = u64::from_str_radix(value, 16).map_err(|_| invalid_data(message))?;
    Ok(f64::from_bits(bits))
}

fn marker_value<'a>(marker: &'a str, key: &str) -> Option<&'a str> {
    marker.lines().find_map(|line| {
        let (candidate, value) = line.split_once('=')?;
        (candidate == key).then_some(value)
    })
}

fn marker_usize(marker: &str, key: &str) -> io::Result<usize> {
    marker_value(marker, key)
        .ok_or(invalid_data("missing U2 pair-result marker field"))?
        .parse::<usize>()
        .map_err(|_| invalid_data("invalid U2 pair-result marker integer"))
}

fn parse_usize(value: &str, message: &'static str) -> io::Result<usize> {
    value.parse::<usize>().map_err(|_| invalid_data(message))
}

fn parse_u64(value: &str, message: &'static str) -> io::Result<u64> {
    value.parse::<u64>().map_err(|_| invalid_data(message))
}

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}
