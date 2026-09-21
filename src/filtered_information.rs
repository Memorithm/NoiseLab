//! Information retention across preregistered frequency-selective filters.
//!
//! Stage 0.5 measures how a declared linear frequency-selective
//! transformation `T` changes the frozen equal-width histogram mutual
//! information between a noise-like observation `N` and a discrete
//! hidden/system state `Z`:
//!
//! ```text
//! I(N; Z)
//! I(T(N); Z)
//! delta_I = I(T(N); Z) - I(N; Z)
//! ```
//!
//! Filtering is treated as an intervention. Bin edges are frozen from the
//! raw observation before filtering so the comparison audits that intervention
//! rather than re-tuning the estimator after seeing the filtered series.
//! This module does not establish causality, mechanism, or a general claim
//! that noise is information.

use crate::information::{
    mutual_information_from_quantized, observation_range,
    quantize_equal_width_with_reference_range, InformationError,
};
use scirust_signal::{
    butter_highpass_sos, butter_lowpass_sos, fir_highpass, fir_lowpass, hanning, lfilter,
    sos_filter, SignalError,
};
use std::error::Error;
use std::fmt::{Display, Formatter};

/// Hard safety ceiling for a single FIR design request.
pub const MAX_FILTERED_INFORMATION_FIR_TAPS: usize = 4_096;

/// Hard safety ceiling for a single Butterworth order request.
pub const MAX_FILTERED_INFORMATION_BUTTER_ORDER: usize = 32;

/// Preregistered linear frequency-selective filter applied as intervention `T`.
///
/// Cutoffs are normalized to Nyquist exactly as in the pinned SciRust filter
/// primitives: `1.0` means Nyquist. FIR designs use a Hann window of length
/// `numtaps`. Band-pass designs are an explicit SciRust high-pass then
/// low-pass cascade; SciRust does not currently expose a dedicated band-pass
/// designer at the pinned revision.
#[derive(Debug, Clone, PartialEq)]
pub enum FrequencySelectiveFilterSpec {
    /// Exact copy; known-answer control for zero information change.
    Identity,
    /// Windowed-sinc FIR low-pass via SciRust `fir_lowpass` + `lfilter`.
    FirLowpass { numtaps: usize, cutoff_nyquist: f64 },
    /// Windowed-sinc FIR high-pass via SciRust `fir_highpass` + `lfilter`.
    ///
    /// SciRust requires an odd `numtaps`.
    FirHighpass { numtaps: usize, cutoff_nyquist: f64 },
    /// FIR band-pass as high-pass then low-pass with the same `numtaps`.
    FirBandpass {
        numtaps: usize,
        low_cutoff_nyquist: f64,
        high_cutoff_nyquist: f64,
    },
    /// Butterworth low-pass SOS via SciRust `butter_lowpass_sos` + `sos_filter`.
    ButterLowpass { order: usize, cutoff_nyquist: f64 },
    /// Butterworth high-pass SOS via SciRust `butter_highpass_sos` + `sos_filter`.
    ButterHighpass { order: usize, cutoff_nyquist: f64 },
    /// Butterworth band-pass as high-pass then low-pass SOS cascade.
    ButterBandpass {
        order: usize,
        low_cutoff_nyquist: f64,
        high_cutoff_nyquist: f64,
    },
}

impl FrequencySelectiveFilterSpec {
    /// Stable protocol label for provenance records.
    pub fn family_label(&self) -> &'static str {
        match self {
            Self::Identity => "identity",
            Self::FirLowpass { .. } => "fir_lowpass",
            Self::FirHighpass { .. } => "fir_highpass",
            Self::FirBandpass { .. } => "fir_bandpass",
            Self::ButterLowpass { .. } => "butter_lowpass",
            Self::ButterHighpass { .. } => "butter_highpass",
            Self::ButterBandpass { .. } => "butter_bandpass",
        }
    }

    /// Compact human-readable provenance string retaining all declared parameters.
    pub fn provenance_detail(&self) -> String {
        match self {
            Self::Identity => "identity".to_string(),
            Self::FirLowpass {
                numtaps,
                cutoff_nyquist,
            } => format!(
                "fir_lowpass;numtaps={numtaps};cutoff_nyquist={cutoff_nyquist};window=hanning"
            ),
            Self::FirHighpass {
                numtaps,
                cutoff_nyquist,
            } => format!(
                "fir_highpass;numtaps={numtaps};cutoff_nyquist={cutoff_nyquist};window=hanning"
            ),
            Self::FirBandpass {
                numtaps,
                low_cutoff_nyquist,
                high_cutoff_nyquist,
            } => format!(
                "fir_bandpass;numtaps={numtaps};low_cutoff_nyquist={low_cutoff_nyquist};high_cutoff_nyquist={high_cutoff_nyquist};window=hanning;cascade=highpass_then_lowpass"
            ),
            Self::ButterLowpass {
                order,
                cutoff_nyquist,
            } => format!("butter_lowpass;order={order};cutoff_nyquist={cutoff_nyquist}"),
            Self::ButterHighpass {
                order,
                cutoff_nyquist,
            } => format!("butter_highpass;order={order};cutoff_nyquist={cutoff_nyquist}"),
            Self::ButterBandpass {
                order,
                low_cutoff_nyquist,
                high_cutoff_nyquist,
            } => format!(
                "butter_bandpass;order={order};low_cutoff_nyquist={low_cutoff_nyquist};high_cutoff_nyquist={high_cutoff_nyquist};cascade=highpass_then_lowpass"
            ),
        }
    }
}

/// Errors returned by the Stage 0.5 filtered-information diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub enum FilteredInformationError {
    /// Canonical observation/hidden-state validation or quantization failed.
    Information(InformationError),
    /// A FIR tap count exceeded the research-bench safety ceiling.
    TooManyFirTaps { requested: usize, maximum: usize },
    /// A Butterworth order exceeded the research-bench safety ceiling.
    TooManyButterOrder { requested: usize, maximum: usize },
    /// A declared cutoff is outside the SciRust open interval `(0, 1)`.
    InvalidCutoff { name: &'static str, value: f64 },
    /// Band-pass edges must satisfy `low < high` strictly inside `(0, 1)`.
    InvalidBandEdges { low: f64, high: f64 },
    /// SciRust rejected the filter design or application.
    Upstream(String),
}

impl Display for FilteredInformationError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Information(error) => write!(f, "filtered-information input error: {error}"),
            Self::TooManyFirTaps { requested, maximum } => write!(
                f,
                "requested {requested} FIR taps exceeds safety maximum {maximum}"
            ),
            Self::TooManyButterOrder { requested, maximum } => write!(
                f,
                "requested Butterworth order {requested} exceeds safety maximum {maximum}"
            ),
            Self::InvalidCutoff { name, value } => write!(
                f,
                "{name} must be in (0, 1) normalized to Nyquist, got {value}"
            ),
            Self::InvalidBandEdges { low, high } => write!(
                f,
                "band-pass requires 0 < low < high < 1 (Nyquist-normalized), got low={low}, high={high}"
            ),
            Self::Upstream(message) => write!(f, "SciRust rejected filter request: {message}"),
        }
    }
}

impl Error for FilteredInformationError {}

impl From<InformationError> for FilteredInformationError {
    fn from(value: InformationError) -> Self {
        Self::Information(value)
    }
}

impl From<SignalError> for FilteredInformationError {
    fn from(value: SignalError) -> Self {
        Self::Upstream(value.to_string())
    }
}

/// Frozen binning provenance for the Stage 0.5 intervention audit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FrozenEqualWidthBinning {
    pub bins: usize,
    pub reference_min: f64,
    pub reference_max: f64,
}

impl FrozenEqualWidthBinning {
    /// Protocol label for the Stage 0.5 bin-edge freeze policy.
    pub const POLICY: &'static str = "frozen_equal_width_from_raw";
}

/// Stage 0.5 audit of information change under a frequency-selective filter.
#[derive(Debug, Clone, PartialEq)]
pub struct FilteredMutualInformationAudit {
    pub raw_mutual_information_bits: f64,
    pub filtered_mutual_information_bits: f64,
    /// Signed change: filtered minus raw mutual information.
    pub information_delta_bits: f64,
    pub binning: FrozenEqualWidthBinning,
    pub binning_policy: &'static str,
    pub filter_family: &'static str,
    pub filter_provenance: String,
    pub filter: FrequencySelectiveFilterSpec,
}

/// Apply a preregistered frequency-selective filter fail-closed.
///
/// Rejects non-finite samples, invalid cutoffs, oversized designs, and any
/// SciRust design/application failure. Output length matches the input.
pub fn apply_frequency_selective_filter(
    observed: &[f64],
    filter: &FrequencySelectiveFilterSpec,
) -> Result<Vec<f64>, FilteredInformationError> {
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(InformationError::NonFiniteObservation { index, value }.into());
        }
    }
    match filter {
        FrequencySelectiveFilterSpec::Identity => Ok(observed.to_vec()),
        FrequencySelectiveFilterSpec::FirLowpass {
            numtaps,
            cutoff_nyquist,
        } => {
            validate_fir_taps(*numtaps)?;
            validate_cutoff("cutoff_nyquist", *cutoff_nyquist)?;
            let window = hanning(*numtaps);
            let taps = fir_lowpass(*numtaps, *cutoff_nyquist, &window)?;
            Ok(lfilter(&taps, &[1.0], observed)?)
        }
        FrequencySelectiveFilterSpec::FirHighpass {
            numtaps,
            cutoff_nyquist,
        } => {
            validate_fir_taps(*numtaps)?;
            validate_cutoff("cutoff_nyquist", *cutoff_nyquist)?;
            let window = hanning(*numtaps);
            let taps = fir_highpass(*numtaps, *cutoff_nyquist, &window)?;
            Ok(lfilter(&taps, &[1.0], observed)?)
        }
        FrequencySelectiveFilterSpec::FirBandpass {
            numtaps,
            low_cutoff_nyquist,
            high_cutoff_nyquist,
        } => {
            validate_fir_taps(*numtaps)?;
            validate_band_edges(*low_cutoff_nyquist, *high_cutoff_nyquist)?;
            let window = hanning(*numtaps);
            let highpass = fir_highpass(*numtaps, *low_cutoff_nyquist, &window)?;
            let lowpass = fir_lowpass(*numtaps, *high_cutoff_nyquist, &window)?;
            let after_high = lfilter(&highpass, &[1.0], observed)?;
            Ok(lfilter(&lowpass, &[1.0], &after_high)?)
        }
        FrequencySelectiveFilterSpec::ButterLowpass {
            order,
            cutoff_nyquist,
        } => {
            validate_butter_order(*order)?;
            validate_cutoff("cutoff_nyquist", *cutoff_nyquist)?;
            let sections = butter_lowpass_sos(*order, *cutoff_nyquist)?;
            Ok(sos_filter(&sections, observed))
        }
        FrequencySelectiveFilterSpec::ButterHighpass {
            order,
            cutoff_nyquist,
        } => {
            validate_butter_order(*order)?;
            validate_cutoff("cutoff_nyquist", *cutoff_nyquist)?;
            let sections = butter_highpass_sos(*order, *cutoff_nyquist)?;
            Ok(sos_filter(&sections, observed))
        }
        FrequencySelectiveFilterSpec::ButterBandpass {
            order,
            low_cutoff_nyquist,
            high_cutoff_nyquist,
        } => {
            validate_butter_order(*order)?;
            validate_band_edges(*low_cutoff_nyquist, *high_cutoff_nyquist)?;
            let high_sections = butter_highpass_sos(*order, *low_cutoff_nyquist)?;
            let low_sections = butter_lowpass_sos(*order, *high_cutoff_nyquist)?;
            let after_high = sos_filter(&high_sections, observed);
            Ok(sos_filter(&low_sections, &after_high))
        }
    }
}

/// Measure histogram MI before and after a preregistered frequency-selective filter.
///
/// Equal-width bin edges are taken once from the raw observation and reused for
/// the filtered series. Values outside the raw range clamp to the edge bins.
/// Filter parameters and the binning policy must be frozen before outcome
/// inspection. A negative `information_delta_bits` means only that this
/// estimator reports less association after `T`; it is not causal evidence and
/// does not prove that filter removed useful physical information.
pub fn filtered_histogram_mutual_information_bits(
    observed: &[f64],
    hidden_state: &[usize],
    bins: usize,
    filter: &FrequencySelectiveFilterSpec,
) -> Result<FilteredMutualInformationAudit, FilteredInformationError> {
    if observed.is_empty() || hidden_state.is_empty() {
        return Err(InformationError::EmptyInput.into());
    }
    if observed.len() != hidden_state.len() {
        return Err(InformationError::LengthMismatch {
            observed: observed.len(),
            hidden_state: hidden_state.len(),
        }
        .into());
    }
    if bins < 2 {
        return Err(InformationError::TooFewBins { bins }.into());
    }
    for (index, &value) in observed.iter().enumerate() {
        if !value.is_finite() {
            return Err(InformationError::NonFiniteObservation { index, value }.into());
        }
    }

    let (reference_min, reference_max) = observation_range(observed)?;
    let raw_quantized =
        quantize_equal_width_with_reference_range(observed, bins, reference_min, reference_max)?;
    let raw_mutual_information_bits =
        mutual_information_from_quantized(&raw_quantized, hidden_state, bins);

    let filtered = apply_frequency_selective_filter(observed, filter)?;
    if filtered.len() != observed.len() {
        return Err(FilteredInformationError::Upstream(format!(
            "filter changed length from {} to {}",
            observed.len(),
            filtered.len()
        )));
    }
    for (index, &value) in filtered.iter().enumerate() {
        if !value.is_finite() {
            return Err(FilteredInformationError::Upstream(format!(
                "filter produced non-finite sample at index {index}: {value}"
            )));
        }
    }

    let filtered_quantized =
        quantize_equal_width_with_reference_range(&filtered, bins, reference_min, reference_max)?;
    let filtered_mutual_information_bits =
        mutual_information_from_quantized(&filtered_quantized, hidden_state, bins);

    Ok(FilteredMutualInformationAudit {
        raw_mutual_information_bits,
        filtered_mutual_information_bits,
        information_delta_bits: filtered_mutual_information_bits - raw_mutual_information_bits,
        binning: FrozenEqualWidthBinning {
            bins,
            reference_min,
            reference_max,
        },
        binning_policy: FrozenEqualWidthBinning::POLICY,
        filter_family: filter.family_label(),
        filter_provenance: filter.provenance_detail(),
        filter: filter.clone(),
    })
}

fn validate_fir_taps(numtaps: usize) -> Result<(), FilteredInformationError> {
    if numtaps == 0 {
        return Err(FilteredInformationError::Upstream(
            "numtaps must be > 0".into(),
        ));
    }
    if numtaps > MAX_FILTERED_INFORMATION_FIR_TAPS {
        return Err(FilteredInformationError::TooManyFirTaps {
            requested: numtaps,
            maximum: MAX_FILTERED_INFORMATION_FIR_TAPS,
        });
    }
    Ok(())
}

fn validate_butter_order(order: usize) -> Result<(), FilteredInformationError> {
    if order == 0 {
        return Err(FilteredInformationError::Upstream(
            "Butterworth order must be > 0".into(),
        ));
    }
    if order > MAX_FILTERED_INFORMATION_BUTTER_ORDER {
        return Err(FilteredInformationError::TooManyButterOrder {
            requested: order,
            maximum: MAX_FILTERED_INFORMATION_BUTTER_ORDER,
        });
    }
    Ok(())
}

fn validate_cutoff(name: &'static str, value: f64) -> Result<(), FilteredInformationError> {
    if !(value.is_finite() && value > 0.0 && value < 1.0) {
        return Err(FilteredInformationError::InvalidCutoff { name, value });
    }
    Ok(())
}

fn validate_band_edges(low: f64, high: f64) -> Result<(), FilteredInformationError> {
    validate_cutoff("low_cutoff_nyquist", low)?;
    validate_cutoff("high_cutoff_nyquist", high)?;
    if low >= high {
        return Err(FilteredInformationError::InvalidBandEdges { low, high });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn assert_close(actual: f64, expected: f64, tol: f64) {
        assert!(
            (actual - expected).abs() <= tol,
            "actual={actual}, expected={expected}, tol={tol}"
        );
    }

    fn low_frequency_state_coded(n: usize) -> (Vec<f64>, Vec<usize>) {
        // One bit coded as block-constant DC amplitude (pure low frequency),
        // plus a small mid-band distractor independent of Z. A high-pass
        // removes the DC code while retaining the distractor.
        let block = 64usize;
        let mut hidden = Vec::with_capacity(n);
        let mut observed = Vec::with_capacity(n);
        for index in 0..n {
            let state = (index / block) % 2;
            hidden.push(state);
            let amplitude = if state == 0 { -1.0 } else { 1.0 };
            let distractor = 0.05 * (2.0 * PI * 0.35 * index as f64).sin();
            observed.push(amplitude + distractor);
        }
        (observed, hidden)
    }

    #[test]
    fn identity_filter_has_exact_zero_delta_under_frozen_bins() {
        let observed = [-2.0, -1.0, 1.0, 2.0, -2.0, -1.0, 1.0, 2.0];
        let hidden = [0, 0, 1, 1, 0, 0, 1, 1];
        let audit = filtered_histogram_mutual_information_bits(
            &observed,
            &hidden,
            4,
            &FrequencySelectiveFilterSpec::Identity,
        )
        .unwrap();
        assert_eq!(
            audit.raw_mutual_information_bits.to_bits(),
            audit.filtered_mutual_information_bits.to_bits()
        );
        assert_eq!(audit.information_delta_bits.to_bits(), 0.0f64.to_bits());
        assert_eq!(audit.binning_policy, "frozen_equal_width_from_raw");
        assert_eq!(audit.filter_family, "identity");
    }

    #[test]
    fn highpass_destroys_low_frequency_coded_bit() {
        let (observed, hidden) = low_frequency_state_coded(1024);
        let audit = filtered_histogram_mutual_information_bits(
            &observed,
            &hidden,
            8,
            &FrequencySelectiveFilterSpec::FirHighpass {
                numtaps: 65,
                cutoff_nyquist: 0.25,
            },
        )
        .unwrap();
        assert!(
            audit.raw_mutual_information_bits > 0.4,
            "raw MI too small: {}",
            audit.raw_mutual_information_bits
        );
        assert!(
            audit.information_delta_bits < -0.2,
            "expected substantial negative delta, got {}",
            audit.information_delta_bits
        );
        assert!(audit.filtered_mutual_information_bits < audit.raw_mutual_information_bits);
    }

    #[test]
    fn independence_control_stays_near_zero_before_and_after_filter() {
        // Observation independent of Z: constant-amplitude mid-band tone with
        // an unrelated alternating label sequence.
        let n = 512;
        let observed: Vec<f64> = (0..n)
            .map(|index| (2.0 * PI * 0.2 * index as f64).sin())
            .collect();
        let hidden: Vec<usize> = (0..n).map(|index| (index / 3) % 2).collect();
        let audit = filtered_histogram_mutual_information_bits(
            &observed,
            &hidden,
            8,
            &FrequencySelectiveFilterSpec::ButterLowpass {
                order: 4,
                cutoff_nyquist: 0.35,
            },
        )
        .unwrap();
        assert!(
            audit.raw_mutual_information_bits < 0.05,
            "raw MI={}",
            audit.raw_mutual_information_bits
        );
        assert!(
            audit.filtered_mutual_information_bits < 0.05,
            "filtered MI={}",
            audit.filtered_mutual_information_bits
        );
        assert_close(audit.information_delta_bits, 0.0, 0.05);
    }

    #[test]
    fn determinism_is_bit_identical_across_calls() {
        let (observed, hidden) = low_frequency_state_coded(256);
        let filter = FrequencySelectiveFilterSpec::FirBandpass {
            numtaps: 33,
            low_cutoff_nyquist: 0.05,
            high_cutoff_nyquist: 0.2,
        };
        let a = filtered_histogram_mutual_information_bits(&observed, &hidden, 6, &filter).unwrap();
        let b = filtered_histogram_mutual_information_bits(&observed, &hidden, 6, &filter).unwrap();
        assert_eq!(
            a.raw_mutual_information_bits.to_bits(),
            b.raw_mutual_information_bits.to_bits()
        );
        assert_eq!(
            a.filtered_mutual_information_bits.to_bits(),
            b.filtered_mutual_information_bits.to_bits()
        );
        assert_eq!(
            a.information_delta_bits.to_bits(),
            b.information_delta_bits.to_bits()
        );
        assert_eq!(a.filter_provenance, b.filter_provenance);
    }

    #[test]
    fn non_finite_observation_fails_closed() {
        let observed = [0.0, f64::NAN, 1.0];
        let hidden = [0, 1, 0];
        assert!(matches!(
            filtered_histogram_mutual_information_bits(
                &observed,
                &hidden,
                2,
                &FrequencySelectiveFilterSpec::Identity
            ),
            Err(FilteredInformationError::Information(
                InformationError::NonFiniteObservation { index: 1, .. }
            ))
        ));
    }

    #[test]
    fn length_mismatch_fails_closed() {
        let observed = [0.0, 1.0, 2.0];
        let hidden = [0, 1];
        assert_eq!(
            filtered_histogram_mutual_information_bits(
                &observed,
                &hidden,
                2,
                &FrequencySelectiveFilterSpec::Identity
            ),
            Err(FilteredInformationError::Information(
                InformationError::LengthMismatch {
                    observed: 3,
                    hidden_state: 2,
                }
            ))
        );
    }

    #[test]
    fn invalid_cutoff_fails_closed() {
        let observed = [0.0, 1.0, -1.0, 0.5];
        assert_eq!(
            apply_frequency_selective_filter(
                &observed,
                &FrequencySelectiveFilterSpec::FirLowpass {
                    numtaps: 5,
                    cutoff_nyquist: 1.0,
                }
            ),
            Err(FilteredInformationError::InvalidCutoff {
                name: "cutoff_nyquist",
                value: 1.0,
            })
        );
        assert_eq!(
            apply_frequency_selective_filter(
                &observed,
                &FrequencySelectiveFilterSpec::ButterBandpass {
                    order: 2,
                    low_cutoff_nyquist: 0.4,
                    high_cutoff_nyquist: 0.2,
                }
            ),
            Err(FilteredInformationError::InvalidBandEdges {
                low: 0.4,
                high: 0.2,
            })
        );
    }

    #[test]
    fn oversized_fir_and_butter_fail_closed() {
        let observed = [0.0, 1.0, -1.0, 0.5];
        assert_eq!(
            apply_frequency_selective_filter(
                &observed,
                &FrequencySelectiveFilterSpec::FirLowpass {
                    numtaps: MAX_FILTERED_INFORMATION_FIR_TAPS + 1,
                    cutoff_nyquist: 0.2,
                }
            ),
            Err(FilteredInformationError::TooManyFirTaps {
                requested: MAX_FILTERED_INFORMATION_FIR_TAPS + 1,
                maximum: MAX_FILTERED_INFORMATION_FIR_TAPS,
            })
        );
        assert_eq!(
            apply_frequency_selective_filter(
                &observed,
                &FrequencySelectiveFilterSpec::ButterLowpass {
                    order: MAX_FILTERED_INFORMATION_BUTTER_ORDER + 1,
                    cutoff_nyquist: 0.2,
                }
            ),
            Err(FilteredInformationError::TooManyButterOrder {
                requested: MAX_FILTERED_INFORMATION_BUTTER_ORDER + 1,
                maximum: MAX_FILTERED_INFORMATION_BUTTER_ORDER,
            })
        );
    }

    #[test]
    fn even_fir_highpass_is_rejected_by_upstream_contract() {
        let observed = [0.0, 1.0, -1.0, 0.5, 0.25];
        let err = apply_frequency_selective_filter(
            &observed,
            &FrequencySelectiveFilterSpec::FirHighpass {
                numtaps: 8,
                cutoff_nyquist: 0.3,
            },
        )
        .unwrap_err();
        assert!(matches!(err, FilteredInformationError::Upstream(_)));
    }

    #[test]
    fn butterworth_lowpass_preserves_low_frequency_coded_bit_better_than_highpass() {
        let (observed, hidden) = low_frequency_state_coded(1024);
        let keep = filtered_histogram_mutual_information_bits(
            &observed,
            &hidden,
            8,
            &FrequencySelectiveFilterSpec::ButterLowpass {
                order: 4,
                cutoff_nyquist: 0.15,
            },
        )
        .unwrap();
        let destroy = filtered_histogram_mutual_information_bits(
            &observed,
            &hidden,
            8,
            &FrequencySelectiveFilterSpec::ButterHighpass {
                order: 4,
                cutoff_nyquist: 0.15,
            },
        )
        .unwrap();
        assert!(
            keep.filtered_mutual_information_bits > destroy.filtered_mutual_information_bits,
            "lowpass retained={}, highpass retained={}",
            keep.filtered_mutual_information_bits,
            destroy.filtered_mutual_information_bits
        );
        assert!(destroy.information_delta_bits < keep.information_delta_bits);
    }
}
