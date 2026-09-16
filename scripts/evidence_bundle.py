#!/usr/bin/env python3
"""Seal or verify a stable local artifact directory; never grant a scientific verdict.

Examples:
    python3 scripts/evidence_bundle.py seal artifacts/smoke
    python3 scripts/evidence_bundle.py verify artifacts/smoke
    python3 scripts/evidence_bundle.py verify artifacts/smoke \
      --expected-manifest-sha256 "$TRUSTED_MANIFEST_SHA256"

SHA-256 checks bind retained bytes, not their truth or authorship. Keep the
printed manifest digest in an independent trusted record. Supplying that digest
back through ``--expected-manifest-sha256`` detects wholesale replacement of a
bundle plus its manifest, but only to the extent that the caller's expected
digest is itself obtained from an independent trusted channel. Do not use this
tool as an isolation boundary against a concurrently hostile filesystem.
"""
from __future__ import annotations

import argparse
import hashlib
import hmac
import json
import os
from pathlib import Path
import re
import stat

MANIFEST = "bundle-manifest.json"
MAX_FILES = 10_000
MAX_FILE_BYTES = 128 * 1024 * 1024
MAX_TOTAL_BYTES = 512 * 1024 * 1024
MAX_MANIFEST_BYTES = 4 * 1024 * 1024


class BundleError(ValueError):
    """The artifact tree or its manifest is invalid or inconsistent."""


def _regular_bytes(path: Path, limit: int) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    if path.is_symlink():
        raise BundleError(f"symlink rejected: {path}")
    fd = os.open(path, flags)
    with os.fdopen(fd, "rb") as stream:
        before = os.fstat(stream.fileno())
        if not stat.S_ISREG(before.st_mode) or before.st_size > limit:
            raise BundleError(f"not a bounded regular file: {path}")
        data = stream.read(limit + 1)
        after = os.fstat(stream.fileno())
        identity = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
        if len(data) > limit or len(data) != before.st_size or identity(before) != identity(after):
            raise BundleError(f"file changed or exceeded budget: {path}")
        return data


def _paths(root: Path) -> list[str]:
    if root.is_symlink() or not root.is_dir():
        raise BundleError("artifact root must be a real directory")
    paths: list[str] = []
    directory_count = 0

    def reject_walk_error(error: OSError) -> None:
        raise BundleError(f"cannot inventory directory: {error}") from error

    for directory, dirs, files in os.walk(root, followlinks=False, onerror=reject_walk_error):
        directory_count += 1
        if directory_count > MAX_FILES or len(Path(directory).relative_to(root).parts) > 32:
            raise BundleError("directory budget exceeded")
        for name in dirs + files:
            path = Path(directory) / name
            mode = path.lstat().st_mode
            if not (stat.S_ISREG(mode) or stat.S_ISDIR(mode)):
                raise BundleError(f"non-regular entry rejected: {path}")
        for name in files:
            relative = (Path(directory) / name).relative_to(root).as_posix()
            if relative != MANIFEST:
                paths.append(relative)
            if len(paths) > MAX_FILES:
                raise BundleError("file budget exceeded")
    return sorted(paths)


def _inventory(root: Path) -> list[dict[str, object]]:
    paths = _paths(root)
    if not paths:
        raise BundleError("empty artifact bundle")
    records: list[dict[str, object]] = []
    total = 0
    for relative in paths:
        data = _regular_bytes(root / relative, min(MAX_FILE_BYTES, MAX_TOTAL_BYTES - total))
        total += len(data)
        records.append({"path": relative, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()})
    if _paths(root) != paths:
        raise BundleError("artifact file set changed during collection")
    return records


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise BundleError(f"duplicate manifest key: {key}")
        result[key] = value
    return result


def _validate_manifest(document: object) -> list[dict[str, object]]:
    if not isinstance(document, dict) or set(document) != {"schema_version", "purpose", "scientific_evidence", "files"}:
        raise BundleError("invalid manifest shape")
    if type(document["schema_version"]) is not int or document["schema_version"] != 1:
        raise BundleError("unsupported manifest schema")
    if document["purpose"] != "byte-integrity-only" or document["scientific_evidence"] is not False:
        raise BundleError("manifest must not grant a scientific verdict")
    records = document["files"]
    if not isinstance(records, list) or not 0 < len(records) <= MAX_FILES:
        raise BundleError("invalid manifest file count")
    names: list[str] = []
    for record in records:
        if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
            raise BundleError("invalid file record")
        name, size, digest = record["path"], record["bytes"], record["sha256"]
        if not isinstance(name, str) or not name or "\\" in name or "\x00" in name:
            raise BundleError("invalid path")
        if name == MANIFEST or any(part in ("", ".", "..") for part in name.split("/")) or Path(name).is_absolute():
            raise BundleError("non-canonical or escaping path")
        if type(size) is not int or not 0 <= size <= MAX_FILE_BYTES:
            raise BundleError("invalid byte count")
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise BundleError("invalid SHA-256")
        names.append(name)
    if names != sorted(set(names)) or sum(record["bytes"] for record in records) > MAX_TOTAL_BYTES:
        raise BundleError("duplicate/unordered records or total budget exceeded")
    return records


def _validate_expected_manifest_sha256(digest: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise BundleError("expected manifest SHA-256 must be 64 lowercase hexadecimal characters")
    return digest


def verify(root: Path, expected_manifest_sha256: str | None = None) -> str:
    """Verify the complete regular-file set; return the manifest's SHA-256.

    This detects omissions, extra files and byte changes relative to the saved
    manifest. When ``expected_manifest_sha256`` is supplied, the manifest digest
    must also match that independently retained lowercase SHA-256 value. This
    detects wholesale manifest+payload replacement only when the expected digest
    itself comes from a separate trusted channel.
    """
    if expected_manifest_sha256 is not None:
        expected_manifest_sha256 = _validate_expected_manifest_sha256(expected_manifest_sha256)
    _paths(root)  # Reject a symlink root/manifest before reading JSON.
    raw = _regular_bytes(root / MANIFEST, MAX_MANIFEST_BYTES)
    expected = _validate_manifest(json.loads(raw, object_pairs_hook=_unique_object))
    if _inventory(root) != expected:
        raise BundleError("bundle files, lengths or SHA-256 digests do not match")
    digest = hashlib.sha256(raw).hexdigest()
    if expected_manifest_sha256 is not None and not hmac.compare_digest(
        digest, expected_manifest_sha256
    ):
        raise BundleError("manifest SHA-256 does not match trusted expected digest")
    return digest


def seal(root: Path) -> str:
    """Write a new manifest without overwriting one, then verify it.

    Example: ``seal(Path("artifacts/smoke"))`` after all producers close files.
    A failed/partial collection must not be promoted to scientific evidence.
    """
    _paths(root)
    if os.path.lexists(root / MANIFEST):
        raise BundleError("manifest already exists; refusing to reseal")
    document = {"schema_version": 1, "purpose": "byte-integrity-only",
                "scientific_evidence": False, "files": _inventory(root)}
    raw = (json.dumps(document, sort_keys=True, indent=2, ensure_ascii=True) + "\n").encode()
    if len(raw) > MAX_MANIFEST_BYTES:
        raise BundleError("manifest budget exceeded")
    with (root / MANIFEST).open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return verify(root)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("seal", "verify"))
    parser.add_argument("directory", type=Path)
    parser.add_argument(
        "--expected-manifest-sha256",
        metavar="HEX",
        help="trusted independently retained manifest SHA-256; valid only for verify",
    )
    args = parser.parse_args()
    if args.action == "seal" and args.expected_manifest_sha256 is not None:
        parser.error("--expected-manifest-sha256 is valid only with verify")
    try:
        digest = (
            seal(args.directory)
            if args.action == "seal"
            else verify(args.directory, args.expected_manifest_sha256)
        )
    except (OSError, ValueError, RecursionError) as error:
        parser.exit(1, f"bundle rejected: {error}\n")
    print(f"byte_integrity=verified scientific_evidence=false manifest_sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
