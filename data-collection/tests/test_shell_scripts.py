"""Tests for the shell entry points (syntax, CLI contract, safety).

These tests never hit the network and never replace the real dataset.
"""

from __future__ import annotations

import subprocess
import unittest
from pathlib import Path

import helpers

BASE_DIR = helpers.BASE_DIR
SCRIPTS = [
    BASE_DIR.parent / "download.sh",
    BASE_DIR / "download.sh",
    BASE_DIR / "scripts" / "validate_dataset.sh",
]
PROCESSED_MANIFEST = BASE_DIR / "data" / "processed" / "manifest.json"


def run(args: list[str], timeout: int = 120) -> subprocess.CompletedProcess:
    return subprocess.run(
        args, cwd=BASE_DIR.parent, capture_output=True, text=True, timeout=timeout
    )


class SyntaxTests(unittest.TestCase):
    def test_scripts_are_executable(self):
        for script in SCRIPTS:
            self.assertTrue(script.exists(), script)
            self.assertTrue(script.stat().st_mode & 0o111, f"{script} is not executable")

    def test_bash_syntax(self):
        for script in SCRIPTS:
            result = subprocess.run(["bash", "-n", str(script)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, f"{script}: {result.stderr}")


class CliContractTests(unittest.TestCase):
    def test_download_help_exits_zero(self):
        result = run(["bash", "download.sh", "--help"])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--taxon", result.stdout)
        self.assertIn("--file", result.stdout)

    def test_unknown_argument_is_rejected(self):
        result = run(["bash", "download.sh", "--not-a-flag"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unknown argument", result.stderr)

    def test_missing_local_file_is_rejected_without_touching_dataset(self):
        before = PROCESSED_MANIFEST.stat().st_mtime_ns if PROCESSED_MANIFEST.exists() else None
        result = run(["bash", "download.sh", "--file", "/nonexistent/hbv.fasta"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("input file not found", result.stderr)
        after = PROCESSED_MANIFEST.stat().st_mtime_ns if PROCESSED_MANIFEST.exists() else None
        self.assertEqual(before, after, "a failed run must not touch the processed dataset")

    def test_validator_requires_dataset(self):
        result = run(["bash", "data-collection/scripts/validate_dataset.sh"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--dataset", result.stdout + result.stderr)

    def test_validator_accepts_fixture(self):
        result = run(
            [
                "bash",
                "data-collection/scripts/validate_dataset.sh",
                "--dataset",
                str(helpers.SAMPLE_FASTA),
                "--min-mb",
                "0",
            ]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('"records": 2', result.stdout)


class PythonModuleCliTests(unittest.TestCase):
    def test_ncbi_client_help(self):
        result = run(
            ["python3", "data-collection/src/ncbi_client.py", "--help"]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--taxon", result.stdout)

    def test_file_loader_requires_a_source(self):
        result = run(["python3", "data-collection/src/file_loader.py"])
        self.assertNotEqual(result.returncode, 0)

    def test_validator_help(self):
        result = run(
            ["python3", "data-collection/src/dataset_validator.py", "--help"]
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--size-policy", result.stdout)


if __name__ == "__main__":
    unittest.main()
