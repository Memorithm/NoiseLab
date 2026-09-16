//! Archived-surrogate replay regressions. Synthetic fixtures are not scientific evidence.

#[path = "../examples/support/u2_capture.rs"]
mod u2_capture;
#[path = "../examples/support/u2_surrogate_arrays.rs"]
mod u2_surrogate_arrays;

use noiselab::u2_manifest::materialize_u2_manifest;
use noiselab::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_PAIRS, U2_FROZEN_SOURCES};
use noiselab::{
    analyze_u2_pair, residual_for_family, StageU2PairAnalysis, StageU2PairRequest,
    StageU2PanelConfig, U2ResidualProvenance, U2ResidualSeries, U2_SCIRUST_REVISION,
};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "noiselab-u2-archive-replay-{}-{}",
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

fn fixtures() -> [U2ResidualSeries; 4] {
    std::array::from_fn(|index| {
        let family = U2_FROZEN_SOURCES[index];
        let residual = (0..128)
            .map(|sample| {
                let x = sample as f64;
                (((index + 1) as f64) * 0.03125 * x).sin() + (index as f64) * 0.001 * x
            })
            .collect::<Vec<_>>();
        U2ResidualSeries {
            family,
            provenance: U2ResidualProvenance {
                family,
                scirust_revision: U2_SCIRUST_REVISION,
                data_seed: index as u64,
                samples: residual.len(),
                burn_in_samples: 0,
                extraction_rule: "synthetic archive-replay fixture only",
                parameter_summary: "synthetic serialization test; not a U2 observation",
            },
            residual,
        }
    })
}

fn fixture_bundle(
    root: &Path,
) -> (
    StageU2PanelConfig,
    Vec<U2SurrogateJob>,
    Vec<StageU2PairAnalysis>,
) {
    let mut config = StageU2PanelConfig::non_scientific_smoke();
    config.scales = vec![1, 2, 4, 8];
    let residuals = fixtures();
    let plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    let jobs = materialize_u2_manifest(plan);
    let pairs = U2_FROZEN_PAIRS
        .iter()
        .copied()
        .enumerate()
        .map(|(pair_index, pair)| {
            analyze_u2_pair(StageU2PairRequest {
                pair_index,
                pair,
                series_a: &residual_for_family(&residuals, pair.left).residual,
                series_b: &residual_for_family(&residuals, pair.right).residual,
                scales: &config.scales,
                alpha: config.alpha,
                jobs: &jobs,
                retain_protocol_failures: true,
            })
            .unwrap()
        })
        .collect::<Vec<_>>();

    fs::create_dir(root).unwrap();
    u2_surrogate_arrays::capture_surrogate_arrays(root, &config, &residuals, &jobs).unwrap();
    u2_capture::capture_surrogate_scores(root, &config, &pairs).unwrap();
    (config, jobs, pairs)
}

#[test]
fn archived_arrays_reproduce_every_retained_score_bit() {
    let directory = TemporaryDirectory::new();
    let bundle = directory.0.join("valid");
    let (config, jobs, _) = fixture_bundle(&bundle);

    let verified = u2_surrogate_arrays::verify_archived_surrogate_scores(&bundle, &config).unwrap();
    assert_eq!(verified, jobs.len());
}

#[test]
fn replay_rejects_a_tampered_stored_score_bit() {
    let directory = TemporaryDirectory::new();
    let bundle = directory.0.join("score-tampered");
    let (config, _, _) = fixture_bundle(&bundle);

    let score_path = bundle.join("surrogate_scores.tsv");
    let contents = fs::read_to_string(&score_path).unwrap();
    let mut lines = contents.lines().map(str::to_owned).collect::<Vec<_>>();
    let fields = lines[1].split('\t').collect::<Vec<_>>();
    let bits = u64::from_str_radix(fields[6], 16).unwrap();
    let mut replacement = fields[..6].join("\t");
    replacement.push('\t');
    replacement.push_str(&format!("{:016x}", bits ^ 1));
    lines[1] = replacement;
    fs::write(&score_path, format!("{}\n", lines.join("\n"))).unwrap();

    let error =
        u2_surrogate_arrays::verify_archived_surrogate_scores(&bundle, &config).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}

#[test]
fn replay_rejects_non_finite_archived_surrogate_bytes() {
    let directory = TemporaryDirectory::new();
    let bundle = directory.0.join("array-tampered");
    let (config, jobs, _) = fixture_bundle(&bundle);

    let first_job = jobs[0];
    let first_file = format!(
        "pair{:02}-{:?}-{:03}.f64le",
        first_job.pair_index, first_job.null_family, first_job.repetition
    );
    let path = bundle.join("surrogate_arrays").join(first_file);
    let mut bytes = fs::read(&path).unwrap();
    bytes[..8].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    fs::write(&path, bytes).unwrap();

    let error =
        u2_surrogate_arrays::verify_archived_surrogate_scores(&bundle, &config).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
}
