"""Tests for src/dataset_validator.py (integrity, format and size checks)."""

from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path

import helpers  # noqa: F401  (sets up sys.path)
from dataset_validator import (
    ValidationError,
    is_fasta_name,
    sha256_file,
    validate_dataset,
)


class FastaNameTests(unittest.TestCase):
    def test_recognised_names(self):
        for name in ("a.fa", "a.fna", "a.fasta", "a.fas", "A.FASTA", "a.fasta.gz", "a.fa.gz"):
            self.assertTrue(is_fasta_name(name), name)

    def test_rejected_names(self):
        for name in ("a.txt", "a.csv", "a.zip", "a.fasta.bz2", "README.md"):
            self.assertFalse(is_fasta_name(name), name)


class ValidatorTests(unittest.TestCase):
    @staticmethod
    def _expected_nt() -> int:
        return sum(
            len(line.strip())
            for line in helpers.TWO_RECORDS.splitlines()
            if line and not line.startswith(">")
        )

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_accepts_plain_fasta(self):
        path = helpers.write_fasta(self.tmp / "ok.fasta")
        report = validate_dataset(path, min_mb=0, max_mb=1024, size_policy="hard")
        self.assertTrue(report["ok"])
        self.assertEqual(report["records"], 2)
        self.assertEqual(report["total_nt"], self._expected_nt())
        self.assertEqual(report["size_policy"], "hard")

    def test_accepts_gzipped_fasta(self):
        import gzip

        plain = helpers.write_fasta(self.tmp / "ok.fasta")
        gz = self.tmp / "ok.fasta.gz"
        with plain.open("rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        report = validate_dataset(gz, min_mb=0, max_mb=1024)
        self.assertEqual(report["records"], 2)
        self.assertEqual(report["total_nt"], self._expected_nt())

    def test_accepts_zip_of_fasta(self):
        archive = helpers.make_zip(
            self.tmp / "ok.zip", {"inside.fasta": helpers.TWO_RECORDS.encode()}
        )
        report = validate_dataset(archive, min_mb=0, max_mb=1024)
        self.assertEqual(report["records"], 2)
        self.assertIn("zip", report)

    def test_rejects_non_fasta_file(self):
        path = self.tmp / "notes.txt"
        path.write_text("this is not a sequence", encoding="utf-8")
        with self.assertRaises(ValidationError):
            validate_dataset(path, min_mb=0)

    def test_rejects_missing_path(self):
        with self.assertRaises(ValidationError):
            validate_dataset(self.tmp / "missing.fasta", min_mb=0)

    def test_rejects_empty_file(self):
        path = self.tmp / "empty.fasta"
        path.write_text("", encoding="utf-8")
        with self.assertRaises(ValidationError):
            validate_dataset(path, min_mb=0)

    def test_rejects_fasta_without_records(self):
        path = self.tmp / "header_only.fasta"
        path.write_text("not a header at all\nACGT\n", encoding="utf-8")
        with self.assertRaises(ValidationError):
            validate_dataset(path, min_mb=0)

    def test_rejects_non_iupac_characters(self):
        path = helpers.write_fasta(self.tmp / "bad.fasta", ">bad\nACGTZ\n")
        with self.assertRaises(ValidationError) as ctx:
            validate_dataset(path, min_mb=0)
        self.assertIn("non-IUPAC", str(ctx.exception))

    def test_rejects_empty_sequence_record(self):
        path = helpers.write_fasta(self.tmp / "empty_rec.fasta", ">a\n>b\nACGT\n")
        with self.assertRaises(ValidationError) as ctx:
            validate_dataset(path, min_mb=0)
        self.assertIn("empty sequence", str(ctx.exception))

    def test_zip_without_fasta_is_rejected(self):
        archive = helpers.make_zip(self.tmp / "nofasta.zip", {"readme.txt": b"hello"})
        with self.assertRaises(ValidationError):
            validate_dataset(archive, min_mb=0)

    def test_hard_policy_fails_below_minimum_size(self):
        path = helpers.write_fasta(self.tmp / "small.fasta")
        with self.assertRaises(ValidationError) as ctx:
            validate_dataset(path, min_mb=100, max_mb=1024, size_policy="hard")
        self.assertIn("below the required minimum", str(ctx.exception))

    def test_warn_policy_records_size_warning(self):
        path = helpers.write_fasta(self.tmp / "small.fasta")
        report = validate_dataset(path, min_mb=100, max_mb=1024, size_policy="warn")
        self.assertTrue(report["ok"])
        self.assertEqual(len(report["warnings"]), 1)
        self.assertIn("below the required minimum", report["warnings"][0])

    def test_hard_policy_fails_above_maximum_size(self):
        path = helpers.write_fasta(self.tmp / "small.fasta")
        with self.assertRaises(ValidationError) as ctx:
            validate_dataset(path, min_mb=0, max_mb=0.000001, size_policy="hard")
        self.assertIn("above the allowed maximum", str(ctx.exception))

    def test_invalid_size_policy_rejected(self):
        path = helpers.write_fasta(self.tmp / "ok.fasta")
        with self.assertRaises(ValueError):
            validate_dataset(path, min_mb=0, size_policy="whatever")


class ChecksumTests(unittest.TestCase):
    def test_sha256_matches_hashlib(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = helpers.write_fasta(Path(tmp) / "x.fasta")
            expected = hashlib.sha256(path.read_bytes()).hexdigest()
            self.assertEqual(sha256_file(path), expected)

    def test_checksum_file_uses_relative_paths(self):
        from dataset_validator import write_checksums

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = helpers.write_fasta(root / "data" / "x.fasta")
            out = root / "checksums.sha256"
            write_checksums([target], out, base=root)
            line = out.read_text(encoding="utf-8").strip().splitlines()[-1]
            digest, shown = line.split("  ", 1)
            self.assertEqual(shown, "data/x.fasta")
            self.assertFalse(Path(shown).is_absolute())
            self.assertEqual(digest, sha256_file(target))


if __name__ == "__main__":
    unittest.main()
