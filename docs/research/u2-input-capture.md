# U2 input capture and reproducibility bundles

This is an execution/integrity surface, not a U2 scientific result. It extends
the qualified U2 runner without changing the SciRust revision, source parameters,
seed roots, descriptors, scales, alpha, surrogate count, p-value rule or frozen
U2 preregistration.

## Capture before comparison

`run_stage_u2_panel_with_capture` validates readiness and the supplied panel
configuration, generates the four sources once, materializes the outcome-blind
job manifest, and calls an immutable input sink **before any pair analysis**. The
exact same residual arrays and job descriptors are then supplied to both null
families. The original `run_stage_u2_panel` delegates to this path with a no-op
sink.

When `NOISELAB_U2_CAPTURE_DIR` is enabled, the report runner persists the input
snapshot and every preregistered realized surrogate array before pair analysis.
The array export uses the canonical `materialize_u2_surrogate_pair` primitive
that is differentially checked against the scoring path. If capture fails, panel
analysis does not start. A complete array marker therefore means the declared
surrogate bytes were written and synchronized; it is **not** a panel-completion
or scientific-result marker.

After analysis, the runner persists the exact per-job convergence scores. It then
persists the six post-analysis pair rows only after the canonical score-only
replay reproduces the stored p-values bit-for-bit and reproduces the frozen
decision. These post-analysis files remain reproducibility evidence, not
scientific acceptance.

An invalid configuration never reaches the sink. A sink error becomes
`StageU2PanelError::Capture` and stops execution before pair analysis. A source
error follows the existing fail-closed generator behavior; capture is not called
until all four sources have been produced. Incremental recovery of partly
generated trajectories is not implemented.

The destination must not already exist when the input snapshot is created. Files
are exclusively created, flushed and synchronized rather than appended to an
old run. Use a stable trusted output directory, not an adversarial shared path.

## Exactly what is retained

| File | Content and meaning |
| --- | --- |
| `capture.tsv` | Schema, workload mode, exact alpha bits, seed roots, scales, repetition count and control seed rules |
| Four `*.f64le` files | Exact binary64 little-endian values from the four post-burn-in residual arrays supplied to pair analysis |
| `sources.tsv` | Family identity, retained length, burn-in, seed, SciRust revision, extraction rule and model parameter summary |
| `descriptors.tsv` | Existing multiscale descriptors encoded as exact binary64 bit patterns |
| `descriptor_errors.tsv` | Numerical descriptor rejections; no failed source is silently removed from the snapshot |
| `surrogate_jobs.tsv` | Complete stable-order pair/null/repetition/seed manifest written before pair analysis |
| `INPUTS_COMPLETE` | All input snapshot files were successfully written; not experiment completion |
| `surrogate_arrays/` | One file per preregistered job; exact left surrogate followed by exact right surrogate as binary64 little-endian values |
| `surrogate_arrays.tsv` | Pair/null/repetition/seed binding, file name, samples per side and exact byte count for every realized surrogate file |
| `SURROGATE_ARRAYS_COMPLETE` | Every preregistered surrogate file and the index were written and synchronized; not completion of pair analysis |
| `surrogate_scores.tsv` | Post-analysis convergence score for every successfully realized preregistered surrogate job, encoded as exact binary64 bits and bound to the frozen job identity |
| `SURROGATE_SCORES_COMPLETE` | Successful score rows passed identity/count/finiteness checks and were synchronized; also records whether array capture is present |
| `pair_results.tsv` | Six frozen pair rows with exact binary64 bits for observed distances, convergence and p-values plus the frozen decision or retained protocol error |
| `PAIR_RESULTS_COMPLETE` | Pair identity/order and successful-row score replay agreed; records exact alpha bits and explicitly remains `scientific_evidence=false` |

The four source binary files contain **post-burn-in extracted residuals, before
multiscale normalization/coarse graining**. They are not full pre-burn-in states,
pre-extraction trajectories or physical measurements. Each canonical source has
8,192 binary64 values (65,536 bytes). Signed zero and subnormal values are
preserved bit-for-bit.

Each surrogate file concatenates exact left and right arrays for one job;
`samples_per_side` is the split point. At the frozen 8,192 samples per source,
each file contains 131,072 payload bytes. The non-scientific 19-surrogate smoke
panel has 228 jobs and therefore 29,884,416 surrogate payload bytes. The
199-surrogate full panel has 2,388 jobs and therefore 312,999,936 surrogate
payload bytes. These are deterministic payload sizes, not runtime-memory,
storage-bandwidth or performance measurements.

`surrogate_jobs.tsv` remains outcome-blind. The array exporter reconstructs the
complete expected manifest from the frozen pair ordering, configured surrogate
count and seed root, and rejects missing, reordered, duplicated or mis-seeded
jobs before publishing a complete array export. The shuffled null initializes
one RNG per job and shuffles left then right with that stream. The spectral null
uses `job.seed` for the left series and the frozen XOR tag for the right series.

After pair analysis, `surrogate_scores.tsv` binds every retained score to its
preregistered job identity. The score exporter does not generate surrogates. A
malformed, reordered, duplicated, mis-seeded or non-finite score set is rejected
before the score-completion marker is written. Partial files from a failed export
remain diagnostic material only.

## Direct replay from retained surrogate arrays

`verify_archived_surrogate_scores` verifies retained score computation without
regenerating source trajectories or null transformations. It consumes the
retained array index, binary64 surrogate files and retained score table together
with the frozen panel configuration.

The replay path fails closed unless:

- the array and score completion markers have the expected non-scientific
  integrity fields;
- the complete `surrogate_arrays.tsv` index exactly matches the preregistered
  pair/null/repetition/seed sequence;
- each deterministic archive file name, split count and byte count agree with its
  index row;
- referenced files are regular files with the exact declared bounded size;
- decoded binary64 values and retained convergence scores are finite;
- score rows reference valid preregistered jobs without duplicates;
- every successful pair contributes either its complete dual-null score set or
  no retained scores at all when protocol execution failed; and
- recomputing `observed_multiscale_comparison` from each retained left/right
  surrogate pair yields a convergence score with **exactly the same IEEE-754 bit
  pattern** as the archived score row.

Its return value is only the number of retained score rows whose exact replay was
verified. Byte-level replay tests use synthetic fixtures and are not U2 outcomes.
This verifies reproducibility of the retained surrogate-score computation. It
does not authenticate the bundle against an attacker who can replace all inputs,
nor does it establish the scientific validity of the frozen protocol.

## Replay of frozen p-values and decisions from retained scores

`replay_u2_statistics_from_scores` is a separate core primitive. For one frozen
pair, it first revalidates the complete preregistered job sequence, including
pair/null/repetition/seed identity and finite score values. It then recomputes the
one-sided empirical `+1` p-values
`(1 + count(surrogate >= observed)) / (R + 1)` independently for both null
families and reapplies the existing frozen `classify_stage_u2` decision rule.

This replay does **not** regenerate surrogate arrays and does not authorize a new
statistical rule. It is useful for checking that retained score evidence still
maps to the same frozen p-values and decision. Successful replay remains an
integrity result, not a scientific claim or permission to consume a protected
holdout.

## Replay-bound retention of pair results

`capture_pair_results` binds the post-analysis six-pair rows to the retained
score evidence before publishing them. It requires `SURROGATE_SCORES_COMPLETE`,
requires the exact frozen pair order, and for each successful pair invokes
`replay_u2_statistics_from_scores`. The retained `p_shuffle` and `p_phase` must
match the replayed binary64 bits exactly and the retained decision must be
identical. Successful rows also require finite observed distances/convergence.

Protocol-failure rows are preserved explicitly. They may not contain p-values,
a decision or partial surrogate-score evidence. The exporter writes
`PAIR_RESULTS_COMPLETE` only after all six rows satisfy these rules. A failed
export may leave a partial `pair_results.tsv` as diagnostic material, but never a
completion marker.

The current capture path writes `pair_results.tsv`; direct parsing/replay of that
persisted TSV back into a typed pair-result object is still a separate boundary.
Likewise, capture does not authenticate the bundle against replacement of both
data and manifest and does not turn exploratory U2 output into scientific
acceptance.

## Non-scientific smoke example

From the repository root, after normal build prerequisites are available:

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

The full workload remains separately explicit under the frozen U2 protocol. A
snapshot, array export, score export, score/statistics replay, pair-result export
or checksum does not by itself become a scientific verdict. No full scientific
campaign is run merely by this reproducibility surface.

## Byte-integrity sealing

After producers close their outputs, `scripts/evidence_bundle.py seal DIR` writes
a new `bundle-manifest.json` with sorted paths, byte lengths and SHA-256 hashes.
`verify DIR` rechecks the complete regular-file set and rejects modified, missing
or additional files. The script uses the Python standard library and is shared
with other retained smoke evidence.

Symlinks, special files, unreadable directories, duplicate JSON keys, duplicate
or escaping paths, unsupported schema and attempts to claim scientific evidence
are rejected. Fixed limits constrain file count, visited directories, depth,
per-file bytes, total bytes and manifest size. A failed seal may leave a partial
manifest for diagnosis; create a new bundle rather than silently resealing it.

This is **byte integrity, not authentication, proof or scientific acceptance**.
Keep the printed manifest SHA-256 in an independent trusted record. An attacker
who can replace both data and manifest can otherwise produce a different valid
bundle.

## Validation and retention

The Rust regressions check capture ordering, sink errors, bit preservation,
no-overwrite behavior, invalid-input rejection and replay of captured smoke
inputs/jobs. Score-export regressions reject a tampered preregistered seed.
Surrogate-array regressions require deterministic materialization, exact
row/byte accounting and manifest binding.

The archived-array replay regressions require a complete retained smoke score set
to reproduce every convergence-score bit directly from disk, reject a modified
stored score bit and reject non-finite tampering in a retained surrogate payload.
Score-statistics replay regressions check exact `+1` p-values, the decision
boundary and malformed/reordered/non-finite score evidence. Pair-result capture
regressions require successful rows to match that replay exactly and reject a
one-bit p-value mutation. These fixtures are synthetic and explicitly
non-scientific.

The `Research report smoke contracts` workflow retains source/protocol copies,
including the pair-result exporter, resolved dependency information,
toolchain/host identity, executable hashes, reports and diagnostics, then seals
the artifact tree. The smoke capture requires `PAIR_RESULTS_COMPLETE` and exactly
six retained pair rows. Those checks concern execution and integrity, not a
positive universality label. GitHub artifact retention is finite; durable
scientific evidence requires separate archival policy and independent digest
retention.

## Deliberate next boundaries

A complete future scientific dossier still needs an approved full-load run,
reviewed six-pair outcomes, retained source-generation failures where applicable
and durable archival outside transient CI storage. Array replay, score replay,
frozen-statistics replay and replay-bound pair-result retention are integrity
capabilities; none is a new scientific decision path. Direct typed replay from
the persisted pair-result artifact and any eventual dossier acceptance remain
separately reviewable boundaries. U0/U1 results and the exploratory U2
multiplicity policy remain unchanged.
