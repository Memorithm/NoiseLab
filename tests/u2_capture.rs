//! Input-capture regressions. Synthetic fixtures are not scientific evidence.

#[path = "../examples/support/u2_capture.rs"]
mod u2_capture;
#[path = "../examples/support/u2_surrogate_arrays.rs"]
mod u2_surrogate_arrays;

use noiselab::u2_manifest::materialize_u2_manifest;
use noiselab::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_PAIRS, U2_FROZEN_SOURCES};
use noiselab::universality_u2_panel::{run_stage_u2_panel_with_capture, StageU2PanelError};
use noiselab::{
    analyze_u2_pair, materialize_u2_surrogate_pair, residual_for_family, StageU2PairRequest,
    StageU2PanelConfig, U2ResidualProvenance, U2ResidualSeries, U2_SCIRUST_REVISION,
};
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

struct TemporaryDirectory(PathBuf);

impl TemporaryDirectory {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "noiselab-capture-{}-{}",
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
        let residual = vec![
            -0.0,
            0.0,
            f64::from_bits(1),
            -1.25,
            2.5,
            -3.0,
            4.0,
            -5.5,
        ];
        U2ResidualSeries {
            family,
            provenance: U2ResidualProvenance {
                family,
                scirust_revision: U2_SCIRUST_REVISION,
                data_seed: index as u64,
                samples: residual.len(),
                burn_in_samples: 0,
                extraction_rule: "synthetic fixture only",
                parameter_summary: "synthetic serialization test; not a U2 observation",
            },
            residual,
        }
    })
}

#[test]
fn invalid_config_never_reaches_capture() {
    let mut config = StageU2PanelConfig::scientific();
    config.data_seed ^= 1;
    let result = run_stage_u2_panel_with_capture(&config, |_, _| {
        panic!("capture must not run before successful validation")
    });
    assert!(matches!(
        result,
        Err(StageU2PanelError::PreregistrationDrift(_))
    ));
}

#[test]
fn failed_sink_stops_before_analysis() {
    let result =
        run_stage_u2_panel_with_capture(&StageU2PanelConfig::non_scientific_smoke(), |_, _| {
            Err("simulated full disk".to_string())
        });
    assert_eq!(
        result,
        Err(StageU2PanelError::Capture(
            "simulated full disk".to_string()
        ))
    );
}

#[test]
fn captured_arrays_and_jobs_reproduce_the_actual_pair_analysis() {
    let config = StageU2PanelConfig::non_scientific_smoke();
    let mut saved: Option<([U2ResidualSeries; 4], Vec<U2SurrogateJob>)> = None;
    let result = run_stage_u2_panel_with_capture(&config, |residuals, jobs| {
        assert_eq!(jobs.len(), 6 * 2 * 19);
        assert!(saved.is_none());
        saved = Some(((*residuals).clone(), jobs.to_vec()));
        Ok(())
    })
    .unwrap();
    let (residuals, jobs) = saved.unwrap();
    assert_eq!(
        result.residual_provenance,
        residuals
            .iter()
            .map(|series| series.provenance.clone())
            .collect::<Vec<_>>()
    );
    for original in &result.pairs {
        let replay = analyze_u2_pair(StageU2PairRequest {
            pair_index: original.pair_index,
            pair: original.pair,
            series_a: &residual_for_family(&residuals, original.pair.left).residual,
            series_b: &residual_for_family(&residuals, original.pair.right).residual,
            scales: &config.scales,
            alpha: config.alpha,
            jobs: &jobs,
            retain_protocol_failures: true,
        })
        .unwrap();
        // Compare bits explicitly: failure rows deliberately contain NaNs.
        assert_eq!(
            original.convergence_score.to_bits(),
            replay.convergence_score.to_bits()
        );
        assert_eq!(
            original.fine_scale_distance.to_bits(),
            replay.fine_scale_distance.to_bits()
        );
        assert_eq!(
            original.terminal_scale_distance.to_bits(),
            replay.terminal_scale_distance.to_bits()
        );
        assert_eq!(original.p_shuffle, replay.p_shuffle);
        assert_eq!(original.p_phase, replay.p_phase);
        assert_eq!(original.decision, replay.decision);
        assert_eq!(original.protocol_error, replay.protocol_error);
    }

    let directory = TemporaryDirectory::new();
    let valid_path = directory.0.join("scores-valid");
    fs::create_dir(&valid_path).unwrap();
    u2_capture::capture_surrogate_scores(&valid_path, &config, &result.pairs).unwrap();
    assert!(valid_path.join("surrogate_scores.tsv").is_file());
    assert!(valid_path.join("SURROGATE_SCORES_COMPLETE").is_file());

    let mut tampered_pairs = result.pairs.clone();
    let tampered_score = tampered_pairs
        .iter_mut()
        .find_map(|pair| pair.surrogate_scores.first_mut())
        .expect("smoke panel must realize at least one surrogate score");
    tampered_score.job.seed ^= 1;
    let tampered_path = directory.0.join("scores-tampered");
    fs::create_dir(&tampered_path).unwrap();
    let error =
        u2_capture::capture_surrogate_scores(&tampered_path, &config, &tampered_pairs).unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert!(!tampered_path.join("SURROGATE_SCORES_COMPLETE").exists());
}

#[test]
fn snapshot_preserves_bits_including_signed_zero_and_subnormals() {
    let directory = TemporaryDirectory::new();
    let path = directory.0.join("inputs");
    let residuals = fixtures();
    u2_capture::capture_inputs(
        &path,
        &StageU2PanelConfig::non_scientific_smoke(),
        &residuals,
        &[],
    )
    .unwrap();
    for series in &residuals {
        let bytes = fs::read(path.join(format!("{:?}.f64le", series.family))).unwrap();
        let expected = series
            .residual
            .iter()
            .flat_map(|value| value.to_bits().to_le_bytes())
            .collect::<Vec<_>>();
        assert_eq!(bytes, expected);
    }
    assert!(path.join("INPUTS_COMPLETE").is_file());
    // These fixtures are too short for the descriptor; rejection is retained.
    assert_eq!(
        fs::read_to_string(path.join("descriptor_errors.tsv"))
            .unwrap()
            .lines()
            .count(),
        5
    );
}

#[test]
fn surrogate_array_capture_is_complete_bit_exact_and_manifest_bound() {
    let directory = TemporaryDirectory::new();
    let path = directory.0.join("inputs");
    let residuals = fixtures();
    let config = StageU2PanelConfig::non_scientific_smoke();
    let plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    let jobs = materialize_u2_manifest(plan);

    u2_capture::capture_inputs(&path, &config, &residuals, &jobs).unwrap();
    u2_surrogate_arrays::capture_surrogate_arrays(&path, &config, &residuals, &jobs).unwrap();

    let index = fs::read_to_string(path.join("surrogate_arrays.tsv")).unwrap();
    assert_eq!(index.lines().count(), jobs.len() + 1);
    let marker = fs::read_to_string(path.join("SURROGATE_ARRAYS_COMPLETE")).unwrap();
    assert!(marker.contains(&format!("jobs={}", jobs.len())));
    assert!(marker.contains(&format!("bytes={}", jobs.len() * 8 * 16)));
    assert!(marker.contains("scientific_evidence=false"));

    let first_job = jobs[0];
    let first_left = residual_for_family(&residuals, first_job.pair.left);
    let first_right = residual_for_family(&residuals, first_job.pair.right);
    let expected =
        materialize_u2_surrogate_pair(first_job, &first_left.residual, &first_right.residual)
            .unwrap();
    let expected_bytes = expected
        .left
        .iter()
        .chain(&expected.right)
        .flat_map(|value| value.to_bits().to_le_bytes())
        .collect::<Vec<_>>();
    let first_file = format!(
        "pair{:02}-{:?}-{:03}.f64le",
        first_job.pair_index, first_job.null_family, first_job.repetition
    );
    assert_eq!(
        fs::read(path.join("surrogate_arrays").join(first_file)).unwrap(),
        expected_bytes
    );

    let tampered_path = directory.0.join("tampered-inputs");
    let mut tampered_jobs = jobs.clone();
    tampered_jobs[0].seed ^= 1;
    u2_capture::capture_inputs(&tampered_path, &config, &residuals, &tampered_jobs).unwrap();
    let error = u2_surrogate_arrays::capture_surrogate_arrays(
        &tampered_path,
        &config,
        &residuals,
        &tampered_jobs,
    )
    .unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::InvalidData);
    assert!(!tampered_path.join("surrogate_arrays").exists());
    assert!(!tampered_path.join("SURROGATE_ARRAYS_COMPLETE").exists());
}

#[test]
fn existing_snapshot_is_never_overwritten() {
    let directory = TemporaryDirectory::new();
    let path = directory.0.join("inputs");
    fs::create_dir(&path).unwrap();
    fs::write(path.join("sentinel"), "preserve").unwrap();
    let error = u2_capture::capture_inputs(
        &path,
        &StageU2PanelConfig::non_scientific_smoke(),
        &fixtures(),
        &[],
    )
    .unwrap_err();
    assert_eq!(error.kind(), std::io::ErrorKind::AlreadyExists);
    assert_eq!(
        fs::read_to_string(path.join("sentinel")).unwrap(),
        "preserve"
    );
    assert!(!path.join("INPUTS_COMPLETE").exists());
}

#[test]
fn invalid_inputs_are_rejected_before_creating_a_directory() {
    let directory = TemporaryDirectory::new();
    let path = directory.0.join("inputs");
    for variant in 0..3 {
        let mut residuals = fixtures();
        match variant {
            0 => residuals[0].residual[0] = f64::NAN,
            1 => residuals[0].provenance.samples += 1,
            _ => residuals[0].family = U2_FROZEN_SOURCES[1],
        }
        assert!(u2_capture::capture_inputs(
            &path,
            &StageU2PanelConfig::non_scientific_smoke(),
            &residuals,
            &[]
        )
        .is_err());
        assert!(!path.exists());
    }
}
