//! Exact realized-surrogate array capture and replay for the Stage U2 runner.
//!
//! These artifacts support byte-level replay and auditing. A complete capture
//! or replay does not by itself establish that panel analysis completed or that
//! any U2 scientific hypothesis is supported.

use noiselab::u2_manifest::materialize_u2_manifest;
use noiselab::u2_plan::{
    U2ExecutionPlan, U2NullFamily, U2SurrogateJob, U2_FROZEN_PAIRS, U2_FROZEN_SOURCES,
};
use noiselab::{
    materialize_u2_surrogate_pair, observed_multiscale_comparison, StageU2PanelConfig,
    U2ResidualSeries, U2SourceFamily,
};
use std::collections::HashSet;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;

const MAX_SURROGATE_ARRAY_BYTES: u64 = 128 * 1024 * 1024;

#[derive(Debug)]
struct ArchivedSurrogateRow {
    job: U2SurrogateJob,
    file_name: String,
    samples_per_side: usize,
    bytes: usize,
}

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

    for (series, expected_family) in residuals.iter().zip(U2_FROZEN_SOURCES) {
        if series.family != expected_family
            || series.provenance.family != expected_family
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
        let realized =
            materialize_u2_surrogate_pair(job.to_owned(), &left.residual, &right.residual)
                .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
        let samples_per_side = realized.left.len();
        let bytes = samples_per_side
            .checked_mul(16)
            .ok_or(invalid_data("U2 surrogate byte-count overflow"))?;
        total_bytes = total_bytes
            .checked_add(bytes as u128)
            .ok_or(invalid_data("U2 surrogate total byte-count overflow"))?;

        let file_name = surrogate_file_name(job);
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
    writeln!(
        marker,
        "layout=left_then_right_ieee754_binary64_little_endian"
    )?;
    writeln!(marker, "scientific_evidence=false")?;
    finish(marker)
}

/// Replay archived U2 surrogate arrays and verify every retained score bit.
///
/// The verifier consumes only the on-disk array/index/score artifacts and the
/// frozen panel configuration. It never regenerates source trajectories or null
/// transformations. It validates the complete array index against the expected
/// preregistered manifest, decodes each referenced binary64 payload, recomputes
/// the multiscale convergence score, and requires exact IEEE-754 bit equality
/// with the corresponding row in `surrogate_scores.tsv`.
///
/// A successful replay is reproducibility/integrity evidence only. It does not
/// authorize a U2 decision, validate a p-value, or turn a stored bundle into a
/// positive scientific result.
pub fn verify_archived_surrogate_scores(
    destination: &Path,
    config: &StageU2PanelConfig,
) -> io::Result<usize> {
    let expected_plan = U2ExecutionPlan {
        pairs: &U2_FROZEN_PAIRS,
        surrogates_per_null: config.mode.surrogates_per_null(),
        seed_root: config.surrogate_seed_root,
    };
    let expected_jobs = materialize_u2_manifest(expected_plan);
    validate_completion_markers(destination, expected_jobs.len())?;
    let archive_rows = read_array_index(destination, &expected_jobs)?;
    replay_score_file(destination, config, expected_plan, &archive_rows)
}

fn validate_completion_markers(destination: &Path, expected_jobs: usize) -> io::Result<()> {
    let arrays_marker = fs::read_to_string(destination.join("SURROGATE_ARRAYS_COMPLETE"))?;
    if marker_value(&arrays_marker, "surrogate_array_export_complete") != Some("true")
        || marker_value(&arrays_marker, "scientific_evidence") != Some("false")
        || marker_usize(&arrays_marker, "jobs")? != expected_jobs
        || marker_usize(&arrays_marker, "files")? != expected_jobs
    {
        return Err(invalid_data("invalid U2 surrogate-array completion marker"));
    }

    let scores_marker = fs::read_to_string(destination.join("SURROGATE_SCORES_COMPLETE"))?;
    if marker_value(&scores_marker, "surrogate_score_export_complete") != Some("true")
        || marker_value(&scores_marker, "scientific_evidence") != Some("false")
        || marker_value(&scores_marker, "surrogate_arrays_persisted") != Some("true")
    {
        return Err(invalid_data("invalid U2 surrogate-score completion marker"));
    }
    Ok(())
}

fn read_array_index(
    destination: &Path,
    expected_jobs: &[U2SurrogateJob],
) -> io::Result<Vec<ArchivedSurrogateRow>> {
    let contents = fs::read_to_string(destination.join("surrogate_arrays.tsv"))?;
    let mut lines = contents.lines();
    if lines.next()
        != Some(
            "pair_index\tleft\tright\tnull_family\trepetition\tseed\tfile\tsamples_per_side\tbytes",
        )
    {
        return Err(invalid_data("invalid U2 surrogate-array index header"));
    }

    let mut rows = Vec::new();
    rows.try_reserve_exact(expected_jobs.len())
        .map_err(|_| invalid_data("unable to allocate U2 surrogate replay index"))?;
    for expected_job in expected_jobs {
        let line = lines
            .next()
            .ok_or(invalid_data("incomplete U2 surrogate-array index"))?;
        rows.push(parse_array_row(line, *expected_job)?);
    }
    if lines.next().is_some() {
        return Err(invalid_data(
            "unexpected extra U2 surrogate-array index row",
        ));
    }
    Ok(rows)
}

fn parse_array_row(line: &str, expected_job: U2SurrogateJob) -> io::Result<ArchivedSurrogateRow> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 9 {
        return Err(invalid_data("invalid U2 surrogate-array index row"));
    }
    let pair_index = parse_usize(fields[0], "invalid U2 array pair index")?;
    let repetition = parse_usize(fields[4], "invalid U2 array repetition")?;
    let seed = parse_u64(fields[5], "invalid U2 array seed")?;
    let samples_per_side = parse_usize(fields[7], "invalid U2 array sample count")?;
    let bytes = parse_usize(fields[8], "invalid U2 array byte count")?;
    let expected_file = surrogate_file_name(&expected_job);
    let expected_bytes = samples_per_side
        .checked_mul(16)
        .ok_or(invalid_data("U2 replay byte-count overflow"))?;

    if pair_index != expected_job.pair_index
        || fields[1] != format!("{:?}", expected_job.pair.left)
        || fields[2] != format!("{:?}", expected_job.pair.right)
        || fields[3] != format!("{:?}", expected_job.null_family)
        || repetition != expected_job.repetition
        || seed != expected_job.seed
        || fields[6] != expected_file
        || samples_per_side == 0
        || bytes != expected_bytes
    {
        return Err(invalid_data("U2 surrogate-array index identity mismatch"));
    }

    Ok(ArchivedSurrogateRow {
        job: expected_job,
        file_name: expected_file,
        samples_per_side,
        bytes,
    })
}

fn replay_score_file(
    destination: &Path,
    config: &StageU2PanelConfig,
    expected_plan: U2ExecutionPlan,
    archive_rows: &[ArchivedSurrogateRow],
) -> io::Result<usize> {
    let contents = fs::read_to_string(destination.join("surrogate_scores.tsv"))?;
    let mut lines = contents.lines();
    if lines.next()
        != Some("pair_index\tleft\tright\tnull_family\trepetition\tseed\tconvergence_score_bits")
    {
        return Err(invalid_data("invalid U2 surrogate-score header"));
    }

    let mut verified = 0usize;
    let mut pair_counts = vec![0usize; U2_FROZEN_PAIRS.len()];
    let mut seen = HashSet::new();
    for line in lines {
        replay_score_row(
            destination,
            config,
            expected_plan,
            archive_rows,
            line,
            &mut seen,
            &mut pair_counts,
        )?;
        verified = verified
            .checked_add(1)
            .ok_or(invalid_data("U2 replay score-count overflow"))?;
    }

    let expected_per_pair = config
        .mode
        .surrogates_per_null()
        .checked_mul(2)
        .ok_or(invalid_data("U2 replay pair-count overflow"))?;
    if pair_counts
        .iter()
        .any(|count| *count != 0 && *count != expected_per_pair)
    {
        return Err(invalid_data("partial successful-pair U2 score archive"));
    }

    let score_marker = fs::read_to_string(destination.join("SURROGATE_SCORES_COMPLETE"))?;
    if marker_usize(&score_marker, "rows")? != verified {
        return Err(invalid_data("U2 score completion row-count mismatch"));
    }
    Ok(verified)
}

fn replay_score_row(
    destination: &Path,
    config: &StageU2PanelConfig,
    expected_plan: U2ExecutionPlan,
    archive_rows: &[ArchivedSurrogateRow],
    line: &str,
    seen: &mut HashSet<(usize, u8, usize)>,
    pair_counts: &mut [usize],
) -> io::Result<()> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 7 {
        return Err(invalid_data("invalid U2 surrogate-score row"));
    }
    let pair_index = parse_usize(fields[0], "invalid U2 score pair index")?;
    let null_family = parse_null_family(fields[3])?;
    let repetition = parse_usize(fields[4], "invalid U2 score repetition")?;
    let seed = parse_u64(fields[5], "invalid U2 score seed")?;
    let stored_bits = u64::from_str_radix(fields[6], 16)
        .map_err(|_| invalid_data("invalid U2 convergence-score bits"))?;
    if !f64::from_bits(stored_bits).is_finite() {
        return Err(invalid_data("non-finite archived U2 convergence score"));
    }

    let expected_job = expected_plan
        .surrogate_job(pair_index, null_family, repetition)
        .ok_or(invalid_data("U2 archived score references an invalid job"))?;
    if fields[1] != format!("{:?}", expected_job.pair.left)
        || fields[2] != format!("{:?}", expected_job.pair.right)
        || seed != expected_job.seed
        || !seen.insert((pair_index, null_family_tag(null_family), repetition))
    {
        return Err(invalid_data("U2 archived score identity mismatch"));
    }

    let archived = archive_rows
        .iter()
        .find(|row| row.job == expected_job)
        .ok_or(invalid_data("missing U2 surrogate array for score row"))?;
    let (left, right) = decode_archived_pair(destination, archived)?;
    let replay = observed_multiscale_comparison(&left, &right, &config.scales)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error.to_string()))?;
    if replay.convergence_score.to_bits() != stored_bits {
        return Err(invalid_data("U2 archived surrogate score bit mismatch"));
    }

    let count = pair_counts
        .get_mut(pair_index)
        .ok_or(invalid_data("U2 score pair index out of range"))?;
    *count = count
        .checked_add(1)
        .ok_or(invalid_data("U2 replay pair score-count overflow"))?;
    Ok(())
}

fn decode_archived_pair(
    destination: &Path,
    archived: &ArchivedSurrogateRow,
) -> io::Result<(Vec<f64>, Vec<f64>)> {
    let path = destination
        .join("surrogate_arrays")
        .join(&archived.file_name);
    let metadata = fs::metadata(&path)?;
    let expected_len = u64::try_from(archived.bytes)
        .map_err(|_| invalid_data("U2 archive byte count does not fit u64"))?;
    if !metadata.is_file()
        || metadata.len() != expected_len
        || metadata.len() > MAX_SURROGATE_ARRAY_BYTES
    {
        return Err(invalid_data("invalid U2 surrogate-array file size"));
    }

    let bytes = fs::read(path)?;
    let side_bytes = archived
        .samples_per_side
        .checked_mul(8)
        .ok_or(invalid_data("U2 replay side byte-count overflow"))?;
    if bytes.len() != archived.bytes || side_bytes.checked_mul(2) != Some(bytes.len()) {
        return Err(invalid_data("invalid U2 surrogate-array payload layout"));
    }
    let (left_bytes, right_bytes) = bytes.split_at(side_bytes);
    Ok((
        decode_binary64(left_bytes, archived.samples_per_side)?,
        decode_binary64(right_bytes, archived.samples_per_side)?,
    ))
}

fn decode_binary64(bytes: &[u8], expected_values: usize) -> io::Result<Vec<f64>> {
    if bytes.len() != expected_values.saturating_mul(8) || !bytes.len().is_multiple_of(8) {
        return Err(invalid_data("invalid binary64 surrogate payload length"));
    }
    let mut values = Vec::new();
    values
        .try_reserve_exact(expected_values)
        .map_err(|_| invalid_data("unable to allocate U2 replay values"))?;
    for word in bytes.as_chunks::<8>().0 {
        let value = f64::from_bits(u64::from_le_bytes(*word));
        if !value.is_finite() {
            return Err(invalid_data("non-finite archived U2 surrogate value"));
        }
        values.push(value);
    }
    Ok(values)
}

fn marker_value<'a>(marker: &'a str, key: &str) -> Option<&'a str> {
    marker.lines().find_map(|line| {
        let (candidate, value) = line.split_once('=')?;
        (candidate == key).then_some(value)
    })
}

fn marker_usize(marker: &str, key: &str) -> io::Result<usize> {
    marker_value(marker, key)
        .ok_or(invalid_data("missing U2 completion-marker field"))?
        .parse::<usize>()
        .map_err(|_| invalid_data("invalid U2 completion-marker integer"))
}

fn parse_usize(value: &str, message: &'static str) -> io::Result<usize> {
    value.parse::<usize>().map_err(|_| invalid_data(message))
}

fn parse_u64(value: &str, message: &'static str) -> io::Result<u64> {
    value.parse::<u64>().map_err(|_| invalid_data(message))
}

fn parse_null_family(value: &str) -> io::Result<U2NullFamily> {
    match value {
        "ShuffledMarginal" => Ok(U2NullFamily::ShuffledMarginal),
        "PhaseRandomizedSpectrum" => Ok(U2NullFamily::PhaseRandomizedSpectrum),
        _ => Err(invalid_data("invalid U2 archived null family")),
    }
}

fn null_family_tag(value: U2NullFamily) -> u8 {
    match value {
        U2NullFamily::ShuffledMarginal => 0,
        U2NullFamily::PhaseRandomizedSpectrum => 1,
    }
}

fn surrogate_file_name(job: &U2SurrogateJob) -> String {
    format!(
        "pair{:02}-{:?}-{:03}.f64le",
        job.pair_index, job.null_family, job.repetition
    )
}

fn series_for_family(
    residuals: &[U2ResidualSeries; 4],
    family: U2SourceFamily,
) -> io::Result<&U2ResidualSeries> {
    residuals
        .iter()
        .find(|series| series.family == family)
        .ok_or(invalid_data(
            "missing U2 residual family for surrogate capture",
        ))
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
