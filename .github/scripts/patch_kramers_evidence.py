from pathlib import Path

path = Path("src/langevin.rs")
text = path.read_text()
old = '''        let calibration =
            calibrate_kramers_matching(model, run, &intensities, &seeds, 0.005).unwrap();
        let peak = calibration
            .interior_peak
            .expect("bistable switching response should be non-monotone");
        let predicted = calibration.predicted_noise_intensity.unwrap();

        assert!(
'''
new = '''        // The local detector only establishes a strict sampled interior maximum.
        // Statistical strength is checked separately below against both scan edges,
        // so the result does not depend on an arbitrary grid-local prominence.
        let calibration =
            calibrate_kramers_matching(model, run, &intensities, &seeds, 0.0).unwrap();
        let peak = calibration
            .interior_peak
            .expect("bistable switching response should be non-monotone");
        let predicted = calibration.predicted_noise_intensity.unwrap();
        let edge_evidence =
            crate::evidence::detect_edge_separated_peak(&calibration.responses, 2.0)
                .unwrap()
                .expect("interior response should remain separated from both scan edges at 2 SE");
        assert_eq!(edge_evidence.coordinate, peak.coordinate);
        assert!(
            edge_evidence.conservative_prominence > 0.0,
            "interior peak is not uncertainty-separated from both scan edges: {edge_evidence:?}"
        );

        assert!(
'''
if old not in text:
    raise SystemExit("expected Langevin test block not found; refusing blind patch")
path.write_text(text.replace(old, new, 1))
