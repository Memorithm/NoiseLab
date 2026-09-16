//! Archived pair-result replay regressions. Synthetic fixtures are not scientific evidence.

#[path = "../examples/support/u2_pair_result_replay.rs"]
mod u2_pair_result_replay;
#[path = "../examples/support/u2_pair_results.rs"]
mod u2_pair_results;

use noiselab::u2_plan::{U2ExecutionPlan, U2NullFamily, U2_FROZEN_PAIRS};
use noiselab::{
    replay_u2_statistics_from_scores, StageU2PairAnalysis, StageU2PanelConfig, U2SurrogateScore,
};
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "noiselab-u2-pair-result-replay-{}-{}",
            std::process::id(),
            NEXT_ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for TemporaryDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn retained_score_marker(path: &std::path::Path) {
    fs::write(
        path.join("SURROGATE_SCORES_COMPLETE"),
        "surrogate_score_export_complete=true\nrows=0\nscientific_evidence=false\n",
    )
    .unwrap();
}

fn synthetic_pairs(config: &StageU2PanelConfig) -> Vec<StageU2PairAnalysis> {
    let plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    U2_FROZEN_PAIRS
        .iter()
        .copied()
        .enumerate()
        .map(|(pair_index, pair)| {
            let observed = 0.5 + pair_index as f64 * 0.001;
            let mut surrogate_scores = Vec::new();
            for null_family in [
                U2NullFamily::ShuffledMarginal,
                U2NullFamily::PhaseRandomizedSpectrum,
            ] {
                for repetition in 0..plan.surrogates_per_null {
                    let convergence_score = match null_family {
                        U2NullFamily::ShuffledMarginal => {
                            if repetition % 2 == 0 {
                                0.75
                            } else {
                                0.25
                            }
                        }
                        U2NullFamily::PhaseRandomizedSpectrum => 0.125,
                    };
                    surrogate_scores.push(U2SurrogateScore {
                        job: plan
                            .surrogate_job(pair_index, null_family, repetition)
                            .unwrap(),
                        convergence_score,
                    });
                }
            }
            let replay = replay_u2_statistics_from_scores(
                plan,
                pair_index,
                observed,
                config.alpha,
                &surrogate_scores,
            )
            .unwrap();
            StageU2PairAnalysis {
                pair_index,
                pair,
                fine_scale_distance: 1.0 + pair_index as f64,
                terminal_scale_distance: 0.25 + pair_index as f64 * 0.01,
                convergence_score: observed,
                p_shuffle: Some(replay.p_shuffle),
                p_phase: Some(replay.p_phase),
                decision: Some(replay.decision),
                surrogate_scores,
                protocol_error: None,
            }
        })
        .collect()
}

fn write_score_archive(path: &std::path::Path, pairs: &[StageU2PairAnalysis]) {
    let mut contents = String::from(
        "pair_index\tleft\tright\tnull_family\trepetition\tseed\tconvergence_score_bits\n",
    );
    let mut rows = 0usize;
    for pair in pairs {
        for score in &pair.surrogate_scores {
            contents.push_str(&format!(
                "{}\t{:?}\t{:?}\t{:?}\t{}\t{}\t{:016x}\n",
                score.job.pair_index,
                score.job.pair.left,
                score.job.pair.right,
                score.job.null_family,
                score.job.repetition,
                score.job.seed,
                score.convergence_score.to_bits()
            ));
            rows += 1;
        }
    }
    fs::write(path.join("surrogate_scores.tsv"), contents).unwrap();
    fs::write(
        path.join("SURROGATE_SCORES_COMPLETE"),
        format!(
            "surrogate_score_export_complete=true\nrows={rows}\nscientific_evidence=false\n"
        ),
    )
    .unwrap();
}

fn fixture(config: &StageU2PanelConfig) -> TemporaryDirectory {
    let directory = TemporaryDirectory::new();
    let pairs = synthetic_pairs(config);
    write_score_archive(&directory.0, &pairs);
    u2_pair_results::capture_pair_results(&directory.0, config, &pairs).unwrap();
    directory
}

#[test]
fn archived_pair_results_replay_to_typed_rows() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = fixture(&config);

    let rows = u2_pair_result_replay::verify_archived_pair_results(&directory.0, &config).unwrap();

    assert_eq!(rows.len(), U2_FROZEN_PAIRS.len());
    assert!(rows.iter().all(|row| row.protocol_error.is_none()));
    assert!(rows.iter().all(|row| row.p_shuffle.is_some()));
    assert!(rows.iter().all(|row| row.p_phase.is_some()));
    assert!(rows.iter().all(|row| row.decision.is_some()));
}

#[test]
fn archived_pair_results_bind_to_archived_scores() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = fixture(&config);

    let rows = u2_pair_result_replay::verify_archived_pair_results_against_scores(
        &directory.0,
        &config,
    )
    .unwrap();

    assert_eq!(rows.len(), U2_FROZEN_PAIRS.len());
    assert!(rows.iter().all(|row| row.protocol_error.is_none()));
}

#[test]
fn archived_score_drift_breaks_pair_result_binding() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = fixture(&config);
    let path = directory.0.join("surrogate_scores.tsv");
    let contents = fs::read_to_string(&path).unwrap();
    let mut lines = contents.lines().map(str::to_owned).collect::<Vec<_>>();
    let mut fields = lines[1].split('\t').map(str::to_owned).collect::<Vec<_>>();
    fields[6] = format!("{:016x}", 0.25f64.to_bits());
    lines[1] = fields.join("\t");
    fs::write(&path, format!("{}\n", lines.join("\n"))).unwrap();

    let error = u2_pair_result_replay::verify_archived_pair_results_against_scores(
        &directory.0,
        &config,
    )
    .unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn archived_decision_drift_is_rejected() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = fixture(&config);
    let path = directory.0.join("pair_results.tsv");
    let contents = fs::read_to_string(&path).unwrap();
    let mut lines = contents.lines().map(str::to_owned).collect::<Vec<_>>();
    let mut fields = lines[1].split('\t').map(str::to_owned).collect::<Vec<_>>();
    fields[8] = "NoObservedConvergence".to_string();
    lines[1] = fields.join("\t");
    fs::write(&path, format!("{}\n", lines.join("\n"))).unwrap();

    let error =
        u2_pair_result_replay::verify_archived_pair_results(&directory.0, &config).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn completion_marker_alpha_drift_is_rejected() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = fixture(&config);
    let path = directory.0.join("PAIR_RESULTS_COMPLETE");
    let contents = fs::read_to_string(&path).unwrap();
    let changed = contents.replace(
        &format!("alpha_bits={:016x}", config.alpha.to_bits()),
        &format!("alpha_bits={:016x}", 0.025f64.to_bits()),
    );
    fs::write(&path, changed).unwrap();

    let error =
        u2_pair_result_replay::verify_archived_pair_results(&directory.0, &config).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn protocol_failure_row_replays_without_fabricated_statistics() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = TemporaryDirectory::new();
    let mut pairs = synthetic_pairs(&config);
    pairs[0].fine_scale_distance = f64::NAN;
    pairs[0].terminal_scale_distance = f64::NAN;
    pairs[0].convergence_score = f64::NAN;
    pairs[0].p_shuffle = None;
    pairs[0].p_phase = None;
    pairs[0].decision = None;
    pairs[0].surrogate_scores.clear();
    pairs[0].protocol_error = Some("synthetic protocol failure".to_string());
    write_score_archive(&directory.0, &pairs);
    u2_pair_results::capture_pair_results(&directory.0, &config, &pairs).unwrap();

    let rows = u2_pair_result_replay::verify_archived_pair_results_against_scores(
        &directory.0,
        &config,
    )
    .unwrap();
    assert_eq!(
        rows[0].protocol_error.as_deref(),
        Some("synthetic protocol failure")
    );
    assert!(rows[0].p_shuffle.is_none());
    assert!(rows[0].p_phase.is_none());
    assert!(rows[0].decision.is_none());
}

#[test]
fn protocol_failure_with_retained_scores_is_rejected_by_disk_binding() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = TemporaryDirectory::new();
    let original_pairs = synthetic_pairs(&config);
    let mut captured_pairs = original_pairs.clone();
    captured_pairs[0].fine_scale_distance = f64::NAN;
    captured_pairs[0].terminal_scale_distance = f64::NAN;
    captured_pairs[0].convergence_score = f64::NAN;
    captured_pairs[0].p_shuffle = None;
    captured_pairs[0].p_phase = None;
    captured_pairs[0].decision = None;
    captured_pairs[0].surrogate_scores.clear();
    captured_pairs[0].protocol_error = Some("synthetic protocol failure".to_string());
    write_score_archive(&directory.0, &original_pairs);
    u2_pair_results::capture_pair_results(&directory.0, &config, &captured_pairs).unwrap();

    let error = u2_pair_result_replay::verify_archived_pair_results_against_scores(
        &directory.0,
        &config,
    )
    .unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn empty_protocol_failure_text_is_rejected_before_capture_completion() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let directory = TemporaryDirectory::new();
    retained_score_marker(&directory.0);
    let mut pairs = synthetic_pairs(&config);
    pairs[0].fine_scale_distance = f64::NAN;
    pairs[0].terminal_scale_distance = f64::NAN;
    pairs[0].convergence_score = f64::NAN;
    pairs[0].p_shuffle = None;
    pairs[0].p_phase = None;
    pairs[0].decision = None;
    pairs[0].surrogate_scores.clear();
    pairs[0].protocol_error = Some(String::new());

    let error = u2_pair_results::capture_pair_results(&directory.0, &config, &pairs).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert!(!directory.0.join("PAIR_RESULTS_COMPLETE").exists());
}
