# U2 trusted bundle-digest verification

Status: **byte-integrity/authentication boundary only; not scientific evidence**.

`scripts/evidence_bundle.py` already seals a retained artifact directory with a canonical `bundle-manifest.json` and verifies the complete regular-file set, file lengths and SHA-256 digests. That local verification cannot detect wholesale replacement when an actor replaces both the payload files and the manifest consistently.

The verifier therefore also accepts an optional independently retained manifest digest:

```bash
python3 scripts/evidence_bundle.py verify artifacts/u2-run \
  --expected-manifest-sha256 "$TRUSTED_MANIFEST_SHA256"
```

The supplied value must be exactly 64 lowercase hexadecimal characters. Verification first performs the existing complete bundle checks, computes the SHA-256 of the manifest bytes, then requires constant-time equality with the expected digest. A mismatch fails closed.

The expected digest is a trust input, not something this repository can manufacture for itself. It is useful only when the caller obtained and retained it through a channel independent of the artifact directory being verified. Saving the expected digest beside the bundle and allowing both to be replaced does not add authentication.

This mechanism does not authenticate authorship, establish who controlled the independent channel, protect a concurrently hostile filesystem, prove correctness of the U2 protocol, open a protected holdout, or grant any scientific verdict. `scientific_evidence=false` remains part of the bundle manifest contract.

For a new retained bundle, the intended sequence is:

```bash
manifest_sha256="$(python3 scripts/evidence_bundle.py seal artifacts/u2-run \
  | sed -n 's/.*manifest_sha256=//p')"
# Store $manifest_sha256 in an independent trusted record.

python3 scripts/evidence_bundle.py verify artifacts/u2-run \
  --expected-manifest-sha256 "$manifest_sha256"
```

The second command is a local demonstration only unless the expected digest has actually crossed an independent trust boundary. Durable archival and attestation policy remain separate work.