//! Exact realized-surrogate array capture for the Stage U2 runner.
//!
//! These artifacts support byte-level replay and auditing. A complete capture
//! does not by itself establish that panel analysis completed or that any U2
//! scientific hypothesis is supported.

use noiselab::u2_manifest::materialize_u2_manifest;
use noiselab::u2_plan::{U2ExecutionPlan, U2SurrogateJob, U2_FROZEN_PAIRS};
use noiselab::{
    materialize_u2_surrogate_pair, StageU2PanelConfig, U2ResidualSeries, U2SourceFamily,
};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Persist the exact left/right surrogate arrays for every preregistered U2 job.
///
/// The destination is the already-created U2 input-capture directory. This
/// function creates a new `surrogate_arrays` child plus an index and writes a
/// completion marker only after every file and the index have been synchronized.
/// Partial files are retained for diagnosis when an error occurs.
pub fn capture_surrogate_arrays(
    destination: &Path,
    config: &StageU2PanelConfig,
    residuals: &[U2ResidualSeries; 4],
    jobs: &[U2SurrogateJob],
) -> io::Result<()> {
    let expected_plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    let expected_jobs = materialize_u2_manifest(expected_plan);
    if jobs != expected_jobs.as_slice() {
        return Err(invalid_data("U2 surrogate job manifest mismatch"));
    }

    for series in residuals {
        if series.family != series.provenance.family
            || series.residual.len() != series.provenance.samples
            || series.residual.iter().any(|value| !value.is_finite())
        {
            return Err(invalid_data("invalid U2 residual for surrogate capture"));
        }
    }

    let arrays_dir = destination.join("surrogate_arrays");
    fs::create_dir(&arrays_dir)?;
    let mut index = new_file(destination, "surrogate_arrays.tsv")?;
    writeln!(
        index,
        "pair_index\tleft\tright\tnull_family\trepetition\tseed\tfile\tsamples_per_side\tbytes"
    )?;

    let mut total_bytes = 0u128;
    for job in jobs {
        let left = series_for_family(residuals, job.pair.left)?;
        let right = series_for_family(residuals, job.pair.right)?;
        let realized = materialize_u2_surrogate_pair(job.to_owned(), &left.residual, &right.residual)
            .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        let samples_per_side = realized.left.len();
        let bytes = samples_per_side
            .checked_mul(16)
            .ok_or(invalid_data("U2 surrogate byte-count overflow"))?;
        total_bytes = total_bytes
            .checked_add(bytes as u128)
            .ok_or(invalid_data("U2 surrogate total byte-count overflow"))?;

        let file_name = format!(
            "pair{:02}-{:?}-{:03}.f64le",
            job.pair_index, job.null_family, job.repetition
        );
        let mut values = new_file(&arrays_dir, &file_name)?;
        for value in realized.left.iter().chain(&realized.right) {
            values.write_all(&value.to_bits().to_le_bytes())?;
        }
        finish(values)?;

        writeln!(
            index,
            "{}\t{:?}\t{:?}\t{:?}\t{}\t{}\t{}\t{}\t{}",
            job.pair_index,
            job.pair.left,
            job.pair.right,
            job.null_family,
            job.repetition,
            job.seed,
            file_name,
            samples_per_side,
            bytes
        )?;
    }
    finish(index)?;

    let mut marker = new_file(destination, "SURROGATE_ARRAYS_COMPLETE")?;
    writeln!(marker, "surrogate_array_export_complete=true")?;
    writeln!(marker, "jobs={}", jobs.len())?;
    writeln!(marker, "files={}", jobs.len())?;
    writeln!(marker, "bytes={total_bytes}")?;
    writeln!(marker, "layout=left_then_right_ieee754_binary64_little_endian")?;
    writeln!(marker, "scientific_evidence=false")?;
    finish(marker)
}

fn series_for_family(
    residuals: &[U2ResidualSeries; 4],
    family: U2SourceFamily,
) -> io::Result<&U2ResidualSeries> {
    residuals
        .iter()
        .find(|series| series.family == family)
        .ok_or(invalid_data("missing U2 residual family for surrogate capture"))
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

fn invalid_data(message: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}
