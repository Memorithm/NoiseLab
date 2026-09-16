//! Pair-result capture regressions. Synthetic fixtures are not scientific evidence.

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
            "noiselab-u2-pair-results-{}-{}",
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
        "surrogate_score_export_complete=true\nscientific_evidence=false\n",
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

#[test]
fn exact_pair_results_are_retained_only_after_score_replay_agrees() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let pairs = synthetic_pairs(&config);
    let directory = TemporaryDirectory::new();
    retained_score_marker(&directory.0);

    u2_pair_results::capture_pair_results(&directory.0, &config, &pairs).unwrap();

    let rows = fs::read_to_string(directory.0.join("pair_results.tsv")).unwrap();
    assert_eq!(rows.lines().count(), U2_FROZEN_PAIRS.len() + 1);
    let first = rows.lines().nth(1).unwrap();
    assert!(first.contains(&format!("{:016x}", pairs[0].convergence_score.to_bits())));
    assert!(first.contains(&format!("{:016x}", pairs[0].p_shuffle.unwrap().to_bits())));

    let marker = fs::read_to_string(directory.0.join("PAIR_RESULTS_COMPLETE")).unwrap();
    assert!(marker.contains(&format!("rows={}", U2_FROZEN_PAIRS.len())));
    assert!(marker.contains(&format!("successful_rows={}", U2_FROZEN_PAIRS.len())));
    assert!(marker.contains("protocol_error_rows=0"));
    assert!(marker.contains("surrogate_scores_persisted=true"));
    assert!(marker.contains("scientific_evidence=false"));
}

#[test]
fn one_bit_p_value_drift_is_rejected_before_completion_marker() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let mut pairs = synthetic_pairs(&config);
    let original = pairs[0].p_shuffle.unwrap();
    pairs[0].p_shuffle = Some(f64::from_bits(original.to_bits() ^ 1));

    let directory = TemporaryDirectory::new();
    retained_score_marker(&directory.0);
    let error = u2_pair_results::capture_pair_results(&directory.0, &config, &pairs).unwrap_err();

    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert!(!directory.0.join("PAIR_RESULTS_COMPLETE").exists());
}

#[test]
fn protocol_failure_rows_remain_explicit_and_cannot_carry_partial_decisions() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let mut pairs = synthetic_pairs(&config);
    pairs[0].fine_scale_distance = f64::NAN;
    pairs[0].terminal_scale_distance = f64::NAN;
    pairs[0].convergence_score = f64::NAN;
    pairs[0].p_shuffle = None;
    pairs[0].p_phase = None;
    pairs[0].decision = None;
    pairs[0].surrogate_scores.clear();
    pairs[0].protocol_error = Some("synthetic protocol failure".to_string());

    let directory = TemporaryDirectory::new();
    retained_score_marker(&directory.0);
    u2_pair_results::capture_pair_results(&directory.0, &config, &pairs).unwrap();
    let marker = fs::read_to_string(directory.0.join("PAIR_RESULTS_COMPLETE")).unwrap();
    assert!(marker.contains("successful_rows=5"));
    assert!(marker.contains("protocol_error_rows=1"));

    let mut invalid = pairs;
    invalid[0].p_shuffle = Some(0.5);
    let invalid_directory = TemporaryDirectory::new();
    retained_score_marker(&invalid_directory.0);
    let error =
        u2_pair_results::capture_pair_results(&invalid_directory.0, &config, &invalid).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert!(!invalid_directory.0.join("PAIR_RESULTS_COMPLETE").exists());
}
