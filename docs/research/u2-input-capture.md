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

When `NOISELAB_U2_CAPTURE_DIR` is enabled, the report runner now persists both the
input snapshot and every preregistered realized surrogate array in that sink,
before pair analysis starts. The array export uses the same canonical
`materialize_u2_surrogate_pair` primitive that is differentially checked against
the current U2 scoring path. If either capture fails, panel analysis does not
start. A complete array marker therefore means the declared surrogate bytes were
written and synchronized; it is **not** a panel-completion or scientific-result
marker.

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
| `surrogate_arrays/` | One file per preregistered job; exact left surrogate followed by exact right surrogate as IEEE-754 binary64 little-endian values |
| `surrogate_arrays.tsv` | Pair/null/repetition/seed binding, file name, samples per side and exact byte count for every realized surrogate file |
| `SURROGATE_ARRAYS_COMPLETE` | Every preregistered surrogate file and the index were written and synchronized; not completion of pair analysis |
| `surrogate_scores.tsv` | Post-analysis convergence score for every successfully realized preregistered surrogate job, encoded as exact binary64 bits and revalidated against the frozen pair/null/repetition/seed identity |
| `SURROGATE_SCORES_COMPLETE` | All successful-pair score rows passed identity/count/finiteness checks and were synchronized; it also records whether the separate array-completion marker is present |

The four source binary files hold **post-burn-in extracted residuals, before the
multiscale analysis's normalization/coarse graining**. They are not the full
pre-burn-in states, pre-extraction trajectories or physical measurements. Each
canonical series has 8,192 binary64 values (65,536 bytes). Signed zero and
subnormal values are preserved bit-for-bit.

Each surrogate file concatenates the exact left and right arrays for one job;
`samples_per_side` in the index supplies the split point. At the frozen 8,192
samples per source, each file contains 131,072 payload bytes. The non-scientific
19-surrogate smoke panel contains 228 jobs and therefore 29,884,416 surrogate
payload bytes. The 199-surrogate full panel contains 2,388 jobs and therefore
312,999,936 surrogate payload bytes. These are deterministic payload sizes, not
runtime memory, storage-bandwidth or performance measurements.

`surrogate_jobs.tsv` remains the outcome-blind manifest produced before analysis.
The array exporter reconstructs the complete expected manifest from the frozen
pair ordering, configured surrogate count and recorded seed root, and rejects any
missing, reordered, duplicated or mis-seeded job set before creating the
`surrogate_arrays` directory. The shuffled null initializes one RNG per job and
shuffles the left series then the right using that stream. The spectral null uses
`job.seed` for the left series and the recorded XOR tag for the right series. Do
not substitute one rule for the other.

After panel analysis completes, `surrogate_scores.tsv` binds every retained score
to the exact preregistered job identity. The score exporter does not generate
surrogates; `SURROGATE_SCORES_COMPLETE` only records whether the separately
written `SURROGATE_ARRAYS_COMPLETE` marker is present in the same capture bundle.
A malformed, reordered, duplicated, mis-seeded or non-finite score set is
rejected before the score-completion marker is written. Partial files from a
failed export remain diagnostic material only.

The existing TSV report retains its schema and 12-decimal rendering. Its error
rows remain error rows, not zero-valued successful observations. The snapshot,
source revision, exact realized-surrogate bytes and exact score bits provide a
route to independent replay and audit; they do not increase the displayed
report's precision retroactively.

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
frozen U2 protocol. It does not turn an input snapshot, array export, score export
or checksum into a scientific verdict. No full scientific campaign is run by
this PR's CI.

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
surrogate-realization regressions require deterministic replay and bit-exact
score agreement with the existing U2 analysis for both null families. Array
capture regressions use explicitly synthetic short residuals to require the full
preregistered manifest, verify a stored surrogate file byte-for-byte against the
canonical materializer, check exact row/byte accounting and reject a mis-seeded
manifest before the array directory or completion marker is created. The Python
suite tests integrity, corruption, omitted/extra files and rejection paths using
explicitly synthetic fixtures.

The existing `Research report smoke contracts` workflow runs U2 capture twice,
requires identical reports/snapshots, checks all four source byte lengths and all
228 smoke job entries, and verifies a repeated destination fails without
modifying the first snapshot. With array retention enabled, the capture bundle
also contains all 228 realized surrogate files plus their index and completion
marker before score export. The workflow retains source/protocol copies,
Cargo.lock, resolved Cargo metadata, toolchain/host identity, executable hashes,
the original U2/FHN smoke reports and diagnostics, then seals the whole artifact
tree. Its checks are about execution and integrity, not a positive universality
label.

The workflow uploads diagnostics even on failure, so not every uploaded artifact
has a complete snapshot, array export, score export or valid seal. Inspect job
completion and verify the manifest. GitHub artifact retention is configured to
**14 days**, not permanent storage: archive a downloaded bundle and its
independently recorded digest before expiry when long-term retention is needed.

## Deliberate next boundaries

A complete future scientific dossier still needs its approved full-load run,
reviewed six-pair outcomes, any source-generation failures and a durable archive.
The next reproducibility increment is direct replay from retained on-disk
surrogate arrays without regenerating either source trajectories or null
transformations, plus verification that stored score bits match that replay.
Do not claim that capability from array persistence alone. U0/U1 results and the
exploratory U2 multiplicity policy remain unchanged.
