//! Typed replay of retained Stage U2 pair-result artifacts.
//!
//! This verifier parses `pair_results.tsv` and its completion marker, validates
//! the frozen six-pair identity/order and exact binary64 encodings, and reapplies
//! the frozen decision function. It does not re-read surrogate-score evidence;
//! score-to-pair binding remains the responsibility of capture-time replay and
//! byte-integrity sealing.

use noiselab::u2_decision::{classify_stage_u2, StageU2Decision};
use noiselab::u2_plan::U2_FROZEN_PAIRS;
use noiselab::StageU2PanelConfig;
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

    let successful = rows.iter().filter(|row| row.protocol_error.is_none()).count();
    let failed = rows.len() - successful;
    if marker_usize(&marker, "successful_rows")? != successful
        || marker_usize(&marker, "protocol_error_rows")? != failed
    {
        return Err(invalid_data("U2 pair-result completion counts mismatch"));
    }
    Ok(rows)
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
    let terminal_scale_distance =
        parse_f64_bits(fields[4], "invalid U2 terminal-distance bits")?;
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

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}
