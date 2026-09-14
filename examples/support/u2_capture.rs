//! Versioned input snapshots for the U2 report example (not scientific verdicts).

use noiselab::u2_plan::{U2SurrogateJob, U2_FROZEN_SOURCES};
use noiselab::{multiscale_trace, StageU2PanelConfig, U2ResidualSeries, SPECTRAL_RIGHT_SEED_TAG};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// Persist the exact inputs supplied to the panel's pre-analysis capture sink.
///
/// The destination must not exist. Partial writes are retained for diagnosis,
/// but only a fully written snapshot receives `INPUTS_COMPLETE`. This marker
/// certifies neither completion of the panel nor a positive scientific result.
pub fn capture_inputs(
    destination: &Path,
    config: &StageU2PanelConfig,
    residuals: &[U2ResidualSeries; 4],
    jobs: &[U2SurrogateJob],
) -> io::Result<()> {
    for (series, family) in residuals.iter().zip(U2_FROZEN_SOURCES) {
        if series.family != family
            || series.provenance.family != family
            || series.residual.len() != series.provenance.samples
            || series.residual.iter().any(|value| !value.is_finite())
        {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid input snapshot"));
        }
    }
    // Never reuse an old artifact directory or overwrite a previous run.
    fs::create_dir(destination)?;
    let mut metadata = new_file(destination, "capture.tsv")?;
    writeln!(metadata, "capture_schema\t1")?;
    writeln!(metadata, "mode\t{:?}", config.mode)?;
    writeln!(metadata, "scientific_evidence\tfalse")?;
    writeln!(metadata, "input_kind\tpost_burn_in_extracted_residual")?;
    writeln!(metadata, "encoding\tieee754-binary64-little-endian")?;
    writeln!(metadata, "scales\t1,2,4,8,16")?;
    writeln!(metadata, "alpha_bits\t{:016x}", config.alpha.to_bits())?;
    writeln!(metadata, "data_seed_root\t{}", config.data_seed)?;
    writeln!(metadata, "surrogate_seed_root\t{}", config.surrogate_seed_root)?;
    writeln!(metadata, "surrogates_per_null\t{}", config.mode.surrogates_per_null())?;
    writeln!(metadata, "spectral_right_seed_tag\t{SPECTRAL_RIGHT_SEED_TAG}")?;
    writeln!(metadata, "shuffle_rule\tone_rng_per_job_left_then_right")?;
    writeln!(metadata, "controls\tjob_manifest_not_surrogate_realizations")?;
    finish(metadata)?;

    let mut sources = new_file(destination, "sources.tsv")?;
    writeln!(sources, "family\tfile\tsamples\tburn_in_samples\tdata_seed\tscirust_revision\textraction_rule\tparameters")?;
    let mut descriptors = new_file(destination, "descriptors.tsv")?;
    writeln!(descriptors, "family\tblock_size\tcoarse_samples\tdropped_samples\tvariance_bits\tlag1_bits\tkurtosis_bits\tsign_change_bits\tscaling_exponent_bits")?;
    let mut errors = new_file(destination, "descriptor_errors.tsv")?;
    writeln!(errors, "family\terror")?;
    for series in residuals {
        let name = format!("{:?}.f64le", series.family);
        let mut values = new_file(destination, &name)?;
        for value in &series.residual {
            values.write_all(&value.to_bits().to_le_bytes())?;
        }
        finish(values)?;
        let provenance = &series.provenance;
        writeln!(sources, "{:?}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            series.family, name, provenance.samples, provenance.burn_in_samples,
            provenance.data_seed, provenance.scirust_revision,
            clean(provenance.extraction_rule), clean(provenance.parameter_summary))?;
        // A numerical descriptor rejection is preserved, not turned into a
        // successful trace and not used to omit a source or stop the panel.
        match multiscale_trace(&series.residual, &config.scales) {
            Ok(trace) => {
                for signature in trace.signatures {
                    writeln!(descriptors, "{:?}\t{}\t{}\t{}\t{:016x}\t{:016x}\t{:016x}\t{:016x}\t{:016x}",
                        series.family, signature.block_size, signature.coarse_samples,
                        signature.dropped_samples, signature.variance.to_bits(),
                        signature.lag1_autocorrelation.to_bits(), signature.excess_kurtosis.to_bits(),
                        signature.sign_change_rate.to_bits(), trace.variance_scaling_exponent.to_bits())?;
                }
            }
            Err(error) => writeln!(errors, "{:?}\t{}", series.family, clean(&error.to_string()))?,
        }
    }
    finish(sources)?;
    finish(descriptors)?;
    finish(errors)?;

    let mut manifest = new_file(destination, "surrogate_jobs.tsv")?;
    writeln!(manifest, "pair_index\tleft\tright\tnull_family\trepetition\tseed")?;
    for job in jobs {
        writeln!(manifest, "{}\t{:?}\t{:?}\t{:?}\t{}\t{}", job.pair_index,
            job.pair.left, job.pair.right, job.null_family, job.repetition, job.seed)?;
    }
    finish(manifest)?;
    let mut marker = new_file(destination, "INPUTS_COMPLETE")?;
    writeln!(marker, "input_snapshot_complete=true\nscientific_evidence=false")?;
    finish(marker)
}

fn new_file(directory: &Path, name: &str) -> io::Result<BufWriter<File>> {
    OpenOptions::new().write(true).create_new(true).open(directory.join(name)).map(BufWriter::new)
}

fn finish(mut writer: BufWriter<File>) -> io::Result<()> {
    writer.flush()?;
    writer.get_ref().sync_all()
}

fn clean(value: &str) -> String {
    value.replace(['\t', '\n', '\r'], " ")
}
