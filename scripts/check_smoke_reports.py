#!/usr/bin/env python3
"""Check report plumbing, never scientific outcomes or full workloads."""
# Copyright 2026 Tarek Zekriti
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0

import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess

REPORTS = (
    ("fluctuation_u2_report", "NOISELAB_U2", "U2"),
    ("fhn_horizon_stage0_report", "NOISELAB_FHN_HORIZON", "FHN_HORIZON_STAGE0"),
)
SOURCES = ("DrivenDampedOscillator", "SemiconductorLaser", "BistableLangevin", "FitzHughNagumo")
PAIRS = [(left, right) for i, left in enumerate(SOURCES) for right in SOURCES[i + 1:]]
U2_HEADER = "left right fine_distance terminal_distance convergence p_shuffle p_phase decision protocol_error".split()
FHN_HEADER = "steps decision noise_amplitude seed observed_spikes required_spikes detail".split()
DECISIONS = {
    "NoObservedConvergence", "CompatibleWithMarginalNull",
    "SpectrumExplainedCandidate", "CrossMechanismCandidate",
}


class ContractError(ValueError):
    """Report structure, mode or repeatability is not the declared contract."""


def require(condition, message):
    if not condition:
        raise ContractError(message)


def validate_report(text, kind):
    metadata, table = {}, []
    for line in text.splitlines():
        if line.startswith("# "):
            key, separator, value = line[2:].partition("\t")
            require(separator and key not in metadata, "missing or duplicate metadata key")
            require(not table, "metadata must precede the table")
            metadata[key] = value
        else:
            require(bool(line), "unexpected empty table line")
            table.append(line)
    expected = {"report_schema": "1", "report_kind": kind,
                "mode": "NonScientificSmoke", "scientific_claim_permitted": "false"}
    require(all(metadata.get(k) == v for k, v in expected.items()), "missing or unsafe workload metadata")
    rows = list(csv.reader(io.StringIO("\n".join(table)), delimiter="\t"))
    require(bool(rows), "missing table")
    header, rows = rows[0], rows[1:]
    require(kind in ("U2", "FHN_HORIZON_STAGE0"), "unknown report kind")
    require(header == (U2_HEADER if kind == "U2" else FHN_HEADER), "unexpected table header")
    require(all(len(row) == len(header) for row in rows), "ragged report table")

    if kind == "U2":
        require(metadata.get("surrogates_per_null") == "19", "not the declared U2 smoke load")
        require([(row[0], row[1]) for row in rows] == PAIRS, "missing, reordered or duplicated U2 pair")
        for row in rows:
            if row[8]:
                # Numerical/protocol failures stay in the matrix. They are not
                # silently dropped or forced into a positive scientific result.
                continue
            numbers = [float(value) for value in row[2:7]]
            require(all(math.isfinite(value) for value in numbers), "nonfinite decided U2 row")
            require(all(0 < value <= 1 for value in numbers[-2:]), "invalid U2 p-value")
            require(row[7] in DECISIONS, "missing or unknown U2 decision")
    else:
        require(metadata.get("summary") == "NonScientificSmoke", "smoke has a scientific summary")
        require(metadata.get("burn_in_steps") == "500", "wrong smoke burn-in")
        require([row[0] for row in rows] == ["2000", "4000", "8000"], "wrong FHN smoke horizons")
        for row in rows:
            require(row[1] == "InsufficientSpikes" or row[1] == "ProtocolFailure"
                    or row[1].startswith("Stage0("), "unknown FHN decision")
            if row[1] == "ProtocolFailure":
                require(bool(row[6]), "protocol failure needs a recorded reason")
            elif row[1] == "InsufficientSpikes":
                require(math.isfinite(float(row[2])) and float(row[2]) >= 0, "invalid noise amplitude")
                require(int(row[3]) >= 0 and 0 <= int(row[4]) < int(row[5]), "invalid spike counts")
    return metadata


def invoke(binary, env, output, label, timeout):
    stdout_path, stderr_path = output / (label + ".tsv"), output / (label + ".stderr")
    with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
        result = subprocess.run([str(binary)], env=env, stdout=stdout, stderr=stderr,
                                timeout=timeout, check=False)
    return result.returncode, stdout_path.read_bytes(), stderr_path.read_bytes()


def check_reports(bin_dir, output):
    output.mkdir(parents=True, exist_ok=True)
    clean_env = {key: value for key, value in os.environ.items() if not key.startswith("NOISELAB_")}
    records = []
    for executable, prefix, kind in REPORTS:
        binary = (bin_dir / executable).resolve(strict=True)
        # Invalid/conflicting flags must fail *before* executing a panel.
        for label, full, smoke, marker in (
            ("conflict", "1", "1", b"ConflictingFlags"),
            ("invalid", "yes", "1", b"InvalidFlag"),
        ):
            env = dict(clean_env, **{prefix + "_FULL": full, prefix + "_SMOKE": smoke})
            code, stdout, stderr = invoke(binary, env, output, executable + "-" + label, 10)
            require(code != 0 and not stdout and marker in stderr, "unsafe CLI flag handling")
        env = dict(clean_env, **{prefix + "_FULL": "0", prefix + "_SMOKE": "1"})
        outputs = []
        for repeat in (1, 2):
            code, stdout, _ = invoke(binary, env, output, f"{executable}-{repeat}", 300)
            require(code == 0, f"{executable}: smoke process failed; see retained stderr")
            validate_report(stdout.decode("utf-8"), kind)
            outputs.append(stdout)
        require(outputs[0] == outputs[1], f"{executable}: repeated reports differ")
        records.append({"executable": executable, "kind": kind,
                        "sha256": hashlib.sha256(outputs[0]).hexdigest(), "identical_repeats": 2})
    manifest = {"schema_version": 1, "validation": "non-scientific-report-contract",
                "scientific_evidence": False, "reports": records}
    with (output / "report-contract.json").open("x") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin-dir", type=Path, default=Path("target/debug/examples"))
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/smoke"))
    args = parser.parse_args()
    result = check_reports(args.bin_dir, args.output_dir)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
