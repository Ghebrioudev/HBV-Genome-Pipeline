"""Tests for src/ncbi_client.py.

No test performs network I/O: urllib is mocked.
"""

from __future__ import annotations

import hashlib
import io
import json
import tempfile
import unittest
import urllib.error
import zipfile
from pathlib import Path
from unittest import mock

import helpers  # noqa: F401  (sets up sys.path)
import ncbi_client
from ncbi_client import (
    API_URL,
    DownloadError,
    build_payload,
    download_package,
    extract_zip,
    sha256_file,
    validate_zip,
)


def make_package_zip(path: Path, members: dict[str, bytes] | None = None) -> Path:
    if members is None:
        members = {
            "ncbi_dataset/data/genomic.fna": helpers.TWO_RECORDS.encode(),
            "ncbi_dataset/data/data_report.jsonl": b'{"accession": "NC_TEST.1"}\n',
            "ncbi_dataset/data/dataset_catalog.json": b"{}",
            "md5sum.txt": b"d41d8cd98f00b204e9800998ecf8427e  x\n",
        }
    return helpers.make_zip(path, members)


class FakeResponse:
    def __init__(self, payload: bytes, headers: dict | None = None):
        self._stream = io.BytesIO(payload)
        self.headers = headers or {}

    def read(self, size: int = -1) -> bytes:
        return self._stream.read(size)

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        return False


class PayloadTests(unittest.TestCase):
    def test_default_payload_shape(self):
        payload = build_payload(["10404"], ["genome"])
        self.assertEqual(payload["taxons"], ["10404"])
        self.assertEqual(payload["include_sequence"], ["GENOME"])
        self.assertFalse(payload["complete_only"])
        self.assertFalse(payload["refseq_only"])
        self.assertEqual(payload["released_since"], "0001-01-01T00:00:00Z")
        self.assertEqual(payload["updated_since"], "0001-01-01T00:00:00Z")

    def test_include_is_upper_cased_and_stripped(self):
        payload = build_payload(["10407"], [" genome ", "protein"])
        self.assertEqual(payload["include_sequence"], ["GENOME", "PROTEIN"])

    def test_taxons_are_stringified(self):
        payload = build_payload([10407], ["genome"])
        self.assertEqual(payload["taxons"], ["10407"])

    def test_empty_taxons_rejected(self):
        with self.assertRaises(ValueError):
            build_payload([], ["genome"])


class ZipTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_valid_package_passes(self):
        archive = make_package_zip(self.tmp / "ok.zip")
        validate_zip(archive)  # must not raise

    def test_non_zip_is_rejected(self):
        path = self.tmp / "not_a.zip"
        path.write_bytes(b"plain text")
        with self.assertRaises(DownloadError):
            validate_zip(path)

    def test_missing_required_member_is_rejected(self):
        archive = make_package_zip(
            self.tmp / "incomplete.zip", {"ncbi_dataset/data/genomic.fna": b">a\nACGT\n"}
        )
        with self.assertRaises(DownloadError) as ctx:
            validate_zip(archive)
        self.assertIn("missing required members", str(ctx.exception))

    def test_extract_creates_expected_tree(self):
        archive = make_package_zip(self.tmp / "ok.zip")
        destination = self.tmp / "extracted"
        extract_zip(archive, destination)
        self.assertTrue((destination / "ncbi_dataset" / "data" / "genomic.fna").exists())

    def test_extract_replaces_existing_directory(self):
        archive = make_package_zip(self.tmp / "ok.zip")
        destination = self.tmp / "extracted"
        destination.mkdir()
        (destination / "stale.txt").write_text("old", encoding="utf-8")
        extract_zip(archive, destination)
        self.assertFalse((destination / "stale.txt").exists())

    def test_path_traversal_is_blocked(self):
        archive = helpers.make_zip(self.tmp / "evil.zip", {"../escaped.txt": b"x"})
        with self.assertRaises(DownloadError) as ctx:
            extract_zip(archive, self.tmp / "out")
        self.assertIn("unsafe path", str(ctx.exception))


class HashTests(unittest.TestCase):
    def test_sha256_matches_hashlib(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "f.bin"
            path.write_bytes(b"hbv")
            self.assertEqual(sha256_file(path), hashlib.sha256(b"hbv").hexdigest())


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.destination = self.tmp / "pkg.zip"
        self.payload = build_payload(["10404"], ["genome"])

    def tearDown(self):
        self._tmp.cleanup()

    @mock.patch("ncbi_client.urllib.request.urlopen")
    def test_successful_download_writes_file(self, urlopen):
        urlopen.return_value = FakeResponse(b"PK\x03\x04payload")
        written = download_package(
            API_URL, self.payload, self.destination, show_progress=False
        )
        self.assertEqual(written, len(b"PK\x03\x04payload"))
        self.assertEqual(self.destination.read_bytes(), b"PK\x03\x04payload")
        self.assertFalse(Path(str(self.destination) + ".part").exists())

    @mock.patch("ncbi_client.time.sleep")
    @mock.patch("ncbi_client.urllib.request.urlopen")
    def test_rate_limit_is_retried(self, urlopen, sleep):
        error = urllib.error.HTTPError(
            API_URL, 429, "Too Many Requests", {"Retry-After": "1"}, None
        )
        urlopen.side_effect = [error, FakeResponse(b"data")]
        written = download_package(
            API_URL, self.payload, self.destination, attempts=3, show_progress=False
        )
        self.assertEqual(written, 4)
        self.assertTrue(self.destination.exists())
        sleep.assert_called()
        self.assertGreaterEqual(sleep.call_args[0][0], 1)

    @mock.patch("ncbi_client.time.sleep")
    @mock.patch("ncbi_client.urllib.request.urlopen")
    def test_non_retryable_error_fails_fast(self, urlopen, sleep):
        urlopen.side_effect = urllib.error.HTTPError(API_URL, 404, "Not Found", {}, None)
        with self.assertRaises(DownloadError):
            download_package(API_URL, self.payload, self.destination, attempts=3, show_progress=False)
        self.assertEqual(urlopen.call_count, 1)
        sleep.assert_not_called()
        self.assertFalse(self.destination.exists())

    @mock.patch("ncbi_client.time.sleep")
    @mock.patch("ncbi_client.urllib.request.urlopen")
    def test_empty_response_is_retried_then_fails(self, urlopen, sleep):
        urlopen.return_value = FakeResponse(b"")
        with self.assertRaises(DownloadError):
            download_package(API_URL, self.payload, self.destination, attempts=2, show_progress=False)
        self.assertEqual(urlopen.call_count, 2)
        self.assertFalse(self.destination.exists())

    @mock.patch("ncbi_client.urllib.request.urlopen")
    def test_request_body_and_headers(self, urlopen):
        urlopen.return_value = FakeResponse(b"data")
        download_package(API_URL, self.payload, self.destination, api_key="KEY", show_progress=False)
        request = urlopen.call_args[0][0]
        self.assertEqual(request.get_method(), "POST")
        self.assertEqual(request.get_header("Api-key"), "KEY")
        body = json.loads(request.data.decode("utf-8"))
        self.assertEqual(body["taxons"], ["10404"])


if __name__ == "__main__":
    unittest.main()
