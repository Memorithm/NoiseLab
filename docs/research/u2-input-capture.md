# U2 input capture and reproducibility bundles

This is an execution/integrity increment, not a U2 scientific result. It extends
the runner qualified in PRs #44 and #45 without changing the SciRust revision,
source parameters, seed roots, descriptor, scales, alpha or surrogate count.
The frozen U2 preregistration remains authoritative.

## Capture before comparison

`run_stage_u2_panel_with_capture` validates readiness and the supplied panel
configuration, generates the four sources once, materializes the job manifest,
and calls an immutable input sink **before any pair analysis**. The exact same
residual arrays and job descriptors are then supplied to both null families.
The original `run_stage_u2_panel` delegates to this path with a no-op sink.

An invalid configuration never reaches the sink. A sink error becomes
`StageU2PanelError::Capture` and stops execution before pair analysis. A source
error still follows the existing generator's fail-closed behavior; capture is
not called until all four sources have been produced. Incremental recovery of
partly generated trajectories is not implemented here.

The example enables file capture only when `NOISELAB_U2_CAPTURE_DIR` is set to a
nonempty path. The path's parent must exist; the destination itself must not.
An existing directory or file is rejected, even when it appears empty. Files
are exclusively created, flushed and synchronized, never appended to an old run.
Use a stable, trusted local output directory, not an adversarial shared directory.

## Exactly what is retained

| File | Content and meaning |
| --- | --- |
| `capture.tsv` | Schema 1, workload mode, exact alpha bits, seed roots, scales, repetition count and control seed rules |
| Four `*.f64le` files | Exact binary64 little-endian values from the four arrays passed to pair analysis |
| `sources.tsv` | Family identity, retained length, burn-in, seed, SciRust revision, extraction rule and model parameter summary |
| `descriptors.tsv` | Existing multiscale descriptors encoded as exact binary64 bit patterns |
| `descriptor_errors.tsv` | Numerical descriptor rejections; no failed source is silently removed from this snapshot |
| `surrogate_jobs.tsv` | Complete stable-order pair/null/repetition/seed manifest written before pair analysis |
| `INPUTS_COMPLETE` | All input snapshot files were successfully written; not completion of the experiment |
| `surrogate_scores.tsv` | Post-analysis convergence score for every successfully realized preregistered surrogate job, encoded as exact binary64 bits and revalidated against the frozen pair/null/repetition/seed identity |
| `SURROGATE_SCORES_COMPLETE` | All successful-pair score rows passed identity/count/finiteness checks and were synchronized; realized surrogate arrays are still not retained |

The four binary files hold **post-burn-in extracted residuals, before the
multiscale analysis's normalization/coarse graining**. They are not the full
pre-burn-in states, pre-extraction trajectories or physical measurements. Each
canonical series has 8,192 binary64 values (65,536 bytes). Signed zero and
subnormal values are preserved bit-for-bit.

`surrogate_jobs.tsv` remains the outcome-blind manifest produced before analysis.
After panel analysis completes, `surrogate_scores.tsv` binds every retained score
to the exact preregistered job identity, including the deterministic seed derived
from the recorded surrogate seed root. The shuffled null initializes one RNG per
job and shuffles the left series then the right using that stream. The spectral
null uses `job.seed` for the left series and the recorded XOR tag for the right
series. Do not substitute one rule for the other.

The library now also exposes `materialize_u2_surrogate_pair` as the canonical
capture/replay primitive for one frozen `U2SurrogateJob`. It rejects unequal
lengths and non-finite inputs, uses the same shuffled-marginal RNG ordering and
spectral seed split as the existing analysis path, and is regression-tested by
recomputing exact convergence-score bits for both null families against the
current pair-analysis implementation. This is a prerequisite for later byte-level
surrogate-array retention; it does **not** mean those arrays are persisted yet.

The score export does **not** retain the realized surrogate arrays themselves and
does not turn smoke output, a scientific-load invocation or a completion marker
into a positive universality result. A malformed, reordered, duplicated,
mis-seeded or non-finite score set is rejected before the completion marker is
written. Partial files from a failed export remain diagnostic material only.

The existing TSV report retains its schema and 12-decimal rendering. Its error
rows remain error rows, not zero-valued successful observations. The snapshot,
source revision and exact score bits provide a route to replay and audit; they do
not increase the displayed report's precision retroactively.

## Non-scientific smoke example

From the repository root, after the normal build prerequisites are available:

```bash
cargo generate-lockfile
cargo build --locked --example fluctuation_u2_report
bundle="$(mktemp -d /tmp/noiselab-u2-smoke.XXXXXX)"
cp Cargo.lock "$bundle/Cargo.lock"
git rev-parse HEAD > "$bundle/source-revision.txt"
rustc -Vv > "$bundle/toolchain.txt"
NOISELAB_U2_FULL=0 NOISELAB_U2_SMOKE=1 \
  NOISELAB_U2_CAPTURE_DIR="$bundle/inputs" \
  target/debug/examples/fluctuation_u2_report \
  > "$bundle/report.tsv" 2> "$bundle/diagnostics.txt" && \
python3 scripts/evidence_bundle.py seal "$bundle" && \
python3 scripts/evidence_bundle.py verify "$bundle"
```

The full workload selector remains separately explicit as documented in the
frozen U2 protocol. It does not turn an input snapshot, score export or checksum
into a scientific verdict. No full scientific campaign is run by this PR's CI.

## Byte-integrity sealing

After all producers close their outputs, `scripts/evidence_bundle.py seal DIR`
writes a new `bundle-manifest.json` with sorted paths, byte lengths and SHA-256
hashes. `verify DIR` rechecks the complete regular-file set, rejecting modified,
missing and additional files. The script uses the Python standard library; it
is also used for the already-retained FHN smoke reports in the same CI bundle.

Symlinks, special files, unreadable directories, duplicate JSON keys, duplicate
or escaping paths, unsupported schema and attempts to claim scientific evidence
are rejected. The fixed limits are 10,000 files, 10,000 visited directories,
depth 32, 128 MiB per file, 512 MiB total and 4 MiB for the manifest. A failed
seal can leave a partial manifest for diagnosis; it must not be overwritten or
silently resealed. Create a new bundle after investigating the failure.

This is **byte integrity, not authentication, proof or scientific acceptance**.
Keep the printed manifest SHA-256 in an independent trusted record. An attacker
who can replace both data and manifest can otherwise produce a different valid
bundle. Empty directories are not authenticated. The tool is not a filesystem
sandbox and does not guarantee safety against concurrent hostile path changes.
The machine-readable `scientific_evidence=false` describes the sealer's lack of
authority, even when used to retain bytes from a separately authorized experiment.

## Validation and retention

The Rust regression suite checks capture ordering, sink errors, bit preservation,
no overwrite, invalid-input rejection and replay of the actual captured smoke
arrays/jobs through all six pair analyses. It compares failure NaNs by bits. The
score-export regressions additionally reject a tampered preregistered seed and
require a valid score set to produce the completion marker. The canonical
surrogate-realization regressions additionally require deterministic replay and
bit-exact score agreement with the existing U2 analysis for both null families.
The Python suite tests integrity, corruption, omitted/extra files and rejection
paths using explicitly synthetic fixtures.

The existing `Research report smoke contracts` workflow also runs U2 capture
twice, requires identical reports/snapshots, checks all four byte lengths and
all 228 smoke job entries, and verifies a repeated destination fails without
modifying the first snapshot. The example now also writes the exact per-surrogate
score table and its completion marker after analysis. The workflow retains
source/protocol copies, Cargo.lock, resolved Cargo metadata, toolchain/host
identity, executable hashes, the original U2/FHN smoke reports and diagnostics,
then seals the whole artifact tree. Its checks are about execution and integrity,
not a positive universality label.

The workflow uploads diagnostics even on failure, so not every uploaded artifact
has a complete snapshot, score export or valid seal. Inspect job completion and
verify the manifest. GitHub artifact retention is configured to **14 days**, not
permanent storage: archive a downloaded bundle and its independently recorded
digest before expiry when long-term retention is needed.

## Deliberate next boundaries

A complete future scientific dossier still needs its approved full-load run,
reviewed six-pair outcomes, any source-generation failures, realized surrogate
arrays when required for independent byte-level replay, and a durable archive.
The canonical realization primitive now removes one semantics ambiguity from the
array-persistence step, but persistence itself remains a separate increment.
Replay directly from an on-disk bundle remains a separate increment. Do not claim
those capabilities from this work. U0/U1 results and the exploratory U2
multiplicity policy remain unchanged.
