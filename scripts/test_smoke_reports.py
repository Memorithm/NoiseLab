"""Synthetic parser/runner contract fixtures; not scientific measurements."""
# Copyright 2026 Tarek Zekriti
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0

import tempfile
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

import check_smoke_reports as checks


def fixture(kind):
    metadata = f"# report_schema\t1\n# report_kind\t{kind}\n# mode\tNonScientificSmoke\n# scientific_claim_permitted\tfalse\n"
    if kind == "U2":
        metadata += "# surrogates_per_null\t19\n"
        metadata += "".join(f"# {key}\t{value}\n" for key, value in checks.U2_SEEDS.items())
        rows = [list(pair) + ["1", "2", "-1", "0.9", "0.8", "NoObservedConvergence", ""]
                for pair in checks.PAIRS]
        header = checks.U2_HEADER
    else:
        metadata += "# summary\tNonScientificSmoke\n# burn_in_steps\t500\n"
        rows = [[steps, "InsufficientSpikes", "0.02", "11", "0", "5", ""]
                for steps in ("2000", "4000", "8000")]
        header = checks.FHN_HEADER
    return metadata + "\t".join(header) + "\n" + "\n".join("\t".join(row) for row in rows) + "\n"


class ContractTests(unittest.TestCase):
    def test_both_negative_smoke_matrices_are_accepted(self):
        for kind in ("U2", "FHN_HORIZON_STAGE0"):
            self.assertEqual(checks.validate_report(fixture(kind), kind)["mode"], "NonScientificSmoke")

    def test_full_or_missing_or_duplicated_metadata_is_rejected(self):
        original = fixture("U2")
        for text in (original.replace("NonScientificSmoke", "Scientific"),
                     original.replace("\tfalse", "\ttrue"),
                     original.replace("# report_schema\t1\n", ""),
                     "# report_schema\t1\n" + original):
            with self.subTest(text=text), self.assertRaises(checks.ContractError):
                checks.validate_report(text, "U2")

    def test_protocol_failures_remain_in_u2_matrix(self):
        text = fixture("U2").replace("1\t2\t-1\t0.9\t0.8\tNoObservedConvergence\t", "NaN\tNaN\tNaN\t\t\t\tretained failure")
        checks.validate_report(text, "U2")

    def test_missing_reordered_or_duplicate_pairs_are_rejected(self):
        lines = fixture("U2").splitlines()
        variants = [lines[:-1], lines[:-1] + [lines[-2]], lines[:-2] + lines[-2:][::-1]]
        for variant in variants:
            with self.assertRaises(checks.ContractError):
                checks.validate_report("\n".join(variant), "U2")

    def test_bad_probabilities_nonfinite_decisions_or_ragged_rows_are_rejected(self):
        for old, new in (("\t0.9\t", "\tNaN\t"), ("\t0.8\t", "\t1.5\t"),
                         ("NoObservedConvergence", "invented"), ("1\t2\t-1", "1\t2")):
            with self.subTest(old=old), self.assertRaises((checks.ContractError, ValueError)):
                checks.validate_report(fixture("U2").replace(old, new), "U2")

    def test_wrong_repetition_budget_is_rejected(self):
        with self.assertRaises(checks.ContractError):
            checks.validate_report(fixture("U2").replace("\t19", "\t199"), "U2")

    def test_full_fhn_horizons_are_not_smoke(self):
        with self.assertRaises(checks.ContractError):
            checks.validate_report(fixture("FHN_HORIZON_STAGE0").replace("2000\t", "40000\t"), "FHN_HORIZON_STAGE0")

    def test_fhn_scientific_summary_is_rejected(self):
        text = fixture("FHN_HORIZON_STAGE0").replace("# summary\tNonScientificSmoke", "# summary\tClassified")
        with self.assertRaises(checks.ContractError):
            checks.validate_report(text, "FHN_HORIZON_STAGE0")

    def test_invocation_uses_argv_and_retains_process_failure(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(checks.subprocess, "run") as run:
            run.return_value.returncode = 23
            code, stdout, stderr = checks.invoke(Path("/example"), {}, Path(directory), "case", 7)
            self.assertEqual((code, stdout, stderr), (23, b"", b""))
            self.assertEqual(run.call_args.args, (["/example"],))
            self.assertEqual(run.call_args.kwargs["timeout"], 7)
            self.assertNotIn("shell", run.call_args.kwargs)

    def test_timeout_does_not_become_success(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(checks.subprocess, "run", side_effect=subprocess.TimeoutExpired("example", 1)):
            with self.assertRaises(subprocess.TimeoutExpired):
                checks.invoke(Path("/example"), {}, Path(directory), "timeout", 1)

    def test_existing_evidence_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "same.tsv").write_text("preserve")
            with self.assertRaises(FileExistsError):
                checks.invoke(Path("/example"), {}, output, "same", 1)
            self.assertEqual((output / "same.tsv").read_text(), "preserve")

    def test_missing_or_changed_seed_metadata_is_rejected(self):
        for key, value in checks.U2_SEEDS.items():
            line = f"# {key}\t{value}\n"
            for replacement in ("", f"# {key}\t0\n", f"# {key}\t{int(value) + 1}\n"):
                with self.subTest(key=key, replacement=replacement), self.assertRaises(checks.ContractError):
                    checks.validate_report(fixture("U2").replace(line, replacement), "U2")

    def fhn_variant(self, variant, amplitude="", seed="", observed="", required="", detail=""):
        lines = fixture("FHN_HORIZON_STAGE0").splitlines()
        lines[-3:] = ["\t".join([steps, variant, amplitude, seed, observed, required, detail])
                      for steps in ("2000", "4000", "8000")]
        return "\n".join(lines) + "\n"

    def test_all_actual_fhn_stage0_variants_are_accepted(self):
        cases = [("Stage0(ProtocolMismatch)", ""), ("Stage0(NoInteriorMinimum)", ""),
                 ("Stage0(Accepted { noise_amplitude: 0.075 })", "0.075000000000"),
                 ("Stage0(OutsideAcceptance { noise_amplitude: 0.3 })", "0.300000000000")]
        for variant, amplitude in cases:
            checks.validate_report(self.fhn_variant(variant, amplitude), "FHN_HORIZON_STAGE0")

    def test_unknown_or_inconsistent_fhn_stage0_payloads_are_rejected(self):
        cases = [("Stage0(Garbage)", ""), ("Stage0(NoInteriorMinimum)trailing", ""),
                 ("Stage0(Accepted)", "0.075000000000"),
                 ("Stage0(Accepted { noise_amplitude: NaN })", "NaN"),
                 ("Stage0(Accepted { noise_amplitude: inf })", "inf"),
                 ("Stage0(Accepted { noise_amplitude: -1.0 })", "-1.000000000000"),
                 ("Stage0(Accepted { noise_amplitude: 0.075 })", "0.080000000000"),
                 ("Stage0(NoInteriorMinimum)", "0.075000000000")]
        for variant, amplitude in cases:
            with self.subTest(variant=variant), self.assertRaises(checks.ContractError):
                checks.validate_report(self.fhn_variant(variant, amplitude), "FHN_HORIZON_STAGE0")

    def test_unused_fhn_columns_are_rejected(self):
        for variant, amplitude in [("Stage0(ProtocolMismatch)", ""),
                                   ("Stage0(Accepted { noise_amplitude: 0.075 })", "0.075000000000")]:
            for extras in ({"seed": "1"}, {"observed": "1"}, {"required": "5"}, {"detail": "unexpected"}):
                with self.subTest(extras=extras), self.assertRaises(checks.ContractError):
                    checks.validate_report(self.fhn_variant(variant, amplitude, **extras), "FHN_HORIZON_STAGE0")

    def test_failure_payloads_are_retained_but_must_be_consistent(self):
        checks.validate_report(self.fhn_variant("ProtocolFailure", detail="retained reason"), "FHN_HORIZON_STAGE0")
        for text in (self.fhn_variant("ProtocolFailure"),
                     self.fhn_variant("ProtocolFailure", seed="1", detail="reason"),
                     self.fhn_variant("InsufficientSpikes", "0.02", "11", "0", "5", "unexpected")):
            with self.assertRaises(checks.ContractError):
                checks.validate_report(text, "FHN_HORIZON_STAGE0")


if __name__ == "__main__":
    unittest.main()
