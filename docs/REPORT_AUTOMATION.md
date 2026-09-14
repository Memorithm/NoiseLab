# Automated research-report contracts

This execution layer checks that the existing U2 and FHN observation-horizon
reporters run reproducibly in their **non-scientific smoke mode**. It does not
run a confirmatory experiment, retune source parameters or establish universality.

## Explicit workload selection

Both executables use the same `report_mode::ReportMode` selector. Their existing
variable prefixes are retained:

- U2: `NOISELAB_U2_FULL` and `NOISELAB_U2_SMOKE`;
- FHN: `NOISELAB_FHN_HORIZON_FULL` and `NOISELAB_FHN_HORIZON_SMOKE`.

Only `1`, `0`, `true` and `false` are accepted (case-insensitive words). Unset
variables mean false; without FULL the default remains non-scientific smoke.
Empty strings, malformed/non-Unicode values and simultaneous FULL=true and
SMOKE=true fail **before** constructing or running the panel. Previously the
SMOKE flag was read and then ignored whenever FULL was true.

Explicit full-load selection still does not certify an experimental result or
authorize access to a protected holdout. Existing protocol prerequisites,
source provenance and scientific review remain necessary.

## TSV report schema 1

Each report starts with metadata lines `# key<TAB>value`, followed by its original
TSV table. Common mandatory metadata is `report_schema`, `report_kind`, `mode`
and the actual boolean `scientific_claim_permitted`. FHN additionally records
`summary` and `burn_in_steps`; U2 records its surrogate count and seed roots.

This corrects ambiguous old preambles: U2 did not serialize the actual mode;
FHN printed field names without the actual scientific-workload boolean.
Consumers must now parse the comment metadata, then parse the following TSV
header and rows. This is a report-format change, not a change to a frozen
scientific decision rule. Error text is kept on one table line by replacing tabs,
CR and LF. Every frozen pair/horizon is retained, including protocol failures.

The legacy `scientific_claim_permitted` field means the scientific workload was
selected; it does not mean that a positive claim passed validation. In smoke
reports it must be false. No result is promoted by parsing this flag alone.

## Automated verification

The `Research report smoke contracts` GitHub Actions workflow runs on relevant
PRs and pushes to main and also supports manual dispatch. There is no additional
cron: recurring development remains with the consolidated Memorithm autopilot.
The workflow uses read-only repository permissions and a GitHub-hosted CPU runner.

It checks the exact source SHA, builds the two existing report executables, then:

1. verifies contradictory and invalid flags fail with no stdout report;
2. executes each explicit smoke workload twice;
3. checks the non-scientific metadata, both U2 seed roots, table structure, all six
   U2 pairs and all three abbreviated FHN horizons;
4. requires byte-identical stdout reports within the same execution environment;
5. writes a machine-readable manifest containing report SHA-256 digests and
   `scientific_evidence=false`.

FHN Stage0 rows must name an actual enum variant with consistent amplitude and
empty unused fields. A scientific negative or a recorded protocol failure is not
turned into a test failure solely for being negative. Missing rows, false mode
labels, malformed rows, process failure, timeout or divergent repeats do fail the
contract check. The validator does not recompute the scientific decision from
the data.

```bash
python3 -m unittest discover -s scripts -p 'test_smoke_reports.py' -v
cargo build --examples
python3 scripts/check_smoke_reports.py --output-dir /tmp/noiselab-smoke-new-run
```

Use a new output directory per run: stdout/stderr and the final manifest use
exclusive file creation and refuse to overwrite prior evidence files. The
validator removes inherited `NOISELAB_*` variables from its child environment
and explicitly selects smoke. It does not invoke a shell to run report binaries.

## Provenance and limitations

The artifact contains source SHA, toolchain/host information, generated
`Cargo.lock`, the actual resolved Cargo dependency graph, stdout, stderr and the
validation manifest when all contract checks succeed. The current repository
has no checked-in Cargo.lock, so dependencies are resolved once and captured for
that run; builds then use `--locked`. Recording a lockfile is not a claim that
separate runs always resolve the same registry graph.

Artifacts are retained for 14 days and uploaded even after a failed step when
files exist. This is temporary CI diagnostic retention, not a permanent archive
of research results. A future scientific evidence package needs durable raw
series, descriptors, surrogate details, protocol identity and complete provenance.

The Python tests use explicitly synthetic contract fixtures. The Rust selector
tests do not mutate process environment or execute full workloads. Neither those
tests nor successful smoke reports constitute a U2 universality or FHN robustness
result. U0/U1 results, U2 exploratory multiplicity policy, FHN acceptance rules,
SciRust pins and FieldLab-facing mathematical APIs remain unchanged.
