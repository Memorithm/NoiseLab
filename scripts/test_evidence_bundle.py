"""Synthetic byte-integrity tests; none of these fixtures is experiment evidence."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import evidence_bundle as bundle


class BundleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "inputs").mkdir()
        (self.root / "inputs" / "raw.f64le").write_bytes(b"\0" * 16)
        (self.root / "report.tsv").write_text("negative or protocol-failure rows retained\n")

    def rewrite_manifest(self, change):
        path = self.root / bundle.MANIFEST
        document = json.loads(path.read_text())
        change(document)
        path.write_text(json.dumps(document))

    def test_seal_verify_and_repeat_digest(self):
        digest = bundle.seal(self.root)
        self.assertEqual(len(digest), 64)
        self.assertEqual(digest, bundle.verify(self.root))
        self.assertEqual(digest, bundle.verify(self.root))

    def test_identical_trees_have_identical_manifest(self):
        digest = bundle.seal(self.root)
        with tempfile.TemporaryDirectory() as directory:
            other = Path(directory)
            (other / "inputs").mkdir()
            (other / "inputs" / "raw.f64le").write_bytes(b"\0" * 16)
            (other / "report.tsv").write_bytes((self.root / "report.tsv").read_bytes())
            self.assertEqual(digest, bundle.seal(other))

    def test_tampered_file_same_size_fails(self):
        bundle.seal(self.root)
        (self.root / "inputs" / "raw.f64le").write_bytes(b"\1" * 16)
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_missing_file_fails(self):
        bundle.seal(self.root)
        (self.root / "report.tsv").unlink()
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_extra_file_fails(self):
        bundle.seal(self.root)
        (self.root / "extra").write_text("extra")
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_no_overwrite_even_when_manifest_is_corrupt(self):
        (self.root / bundle.MANIFEST).write_bytes(b"partial")
        with self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)
        self.assertEqual((self.root / bundle.MANIFEST).read_bytes(), b"partial")

    def test_symlink_file_rejected(self):
        (self.root / "link").symlink_to(self.root / "report.tsv")
        with self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)

    def test_symlink_directory_rejected(self):
        (self.root / "link").symlink_to(self.root / "inputs", target_is_directory=True)
        with self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)

    def test_symlink_root_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            link = Path(directory) / "root"
            link.symlink_to(self.root, target_is_directory=True)
            with self.assertRaises(bundle.BundleError):
                bundle.seal(link)

    def test_symlink_manifest_rejected(self):
        (self.root / bundle.MANIFEST).symlink_to(self.root / "report.tsv")
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_special_file_rejected_without_blocking(self):
        import os
        os.mkfifo(self.root / "fifo")
        with self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)

    def test_empty_bundle_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(bundle.BundleError):
                bundle.seal(Path(directory))

    def test_duplicate_json_keys_rejected(self):
        bundle.seal(self.root)
        path = self.root / bundle.MANIFEST
        path.write_text(path.read_text().replace('"schema_version": 1', '"schema_version": 1, "schema_version": 1'))
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_forged_scientific_verdict_rejected(self):
        bundle.seal(self.root)
        self.rewrite_manifest(lambda d: d.update(scientific_evidence=True))
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_escaping_paths_rejected(self):
        for name in ("../outside", "/absolute", "./x", "a//b", "a/../b", "a\\b", bundle.MANIFEST):
            document = {"schema_version": 1, "purpose": "byte-integrity-only", "scientific_evidence": False,
                        "files": [{"path": name, "bytes": 1, "sha256": "a" * 64}]}
            with self.subTest(name=name), self.assertRaises(bundle.BundleError):
                bundle._validate_manifest(document)

    def test_duplicate_records_rejected(self):
        bundle.seal(self.root)
        self.rewrite_manifest(lambda d: d["files"].append(d["files"][0]))
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_removed_record_cannot_hide_a_file(self):
        bundle.seal(self.root)
        self.rewrite_manifest(lambda d: d["files"].pop())
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_invalid_digest_and_boolean_size_rejected(self):
        for update in ({"sha256": "x" * 64}, {"bytes": True}, {"bytes": -1}):
            document = {"schema_version": 1, "purpose": "byte-integrity-only", "scientific_evidence": False,
                        "files": [{"path": "a", "bytes": 1, "sha256": "a" * 64} | update]}
            with self.subTest(update=update), self.assertRaises(bundle.BundleError):
                bundle._validate_manifest(document)

    def test_unreadable_directory_is_not_silently_omitted(self):
        def failed_walk(*args, **kwargs):
            kwargs["onerror"](PermissionError("denied"))
            return iter(())
        with mock.patch.object(bundle.os, "walk", side_effect=failed_walk):
            with self.assertRaises(bundle.BundleError):
                bundle.seal(self.root)
        self.assertFalse((self.root / bundle.MANIFEST).exists())

    def test_bounded_file_count(self):
        with mock.patch.object(bundle, "MAX_FILES", 1), self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)

    def test_bounded_file_bytes(self):
        with mock.patch.object(bundle, "MAX_FILE_BYTES", 8), self.assertRaises(bundle.BundleError):
            bundle.seal(self.root)

    def test_bounded_manifest_bytes(self):
        bundle.seal(self.root)
        with mock.patch.object(bundle, "MAX_MANIFEST_BYTES", 8), self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)

    def test_reordered_records_rejected(self):
        bundle.seal(self.root)
        self.rewrite_manifest(lambda d: d["files"].reverse())
        with self.assertRaises(bundle.BundleError):
            bundle.verify(self.root)


if __name__ == "__main__":
    unittest.main()
