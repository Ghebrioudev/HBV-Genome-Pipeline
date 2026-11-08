"""Tests for src/file_loader.py (local + NCBI normalisation)."""

from __future__ import annotations

import csv
import gzip
import io
import tempfile
import unittest
from pathlib import Path

import helpers  # noqa: F401  (sets up sys.path)
from file_loader import (
    METADATA_COLUMNS,
    LoaderError,
    normalise_local,
    normalise_ncbi_report,
    normalise_package,
    parse_fasta_headers,
)


class LocalInputTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.out = self.root / "processed"
        self.work = self.root / "work"

    def tearDown(self):
        self._tmp.cleanup()

    def _assert_outputs(self, result, expected_records=2):
        self.assertEqual(result["records"], expected_records)
        sequences = self.out / "sequences.fasta"
        metadata = self.out / "metadata.csv"
        manifest = self.out / "manifest.json"
        self.assertTrue(sequences.exists())
        self.assertTrue(metadata.exists())
        self.assertTrue(manifest.exists())
        self.assertEqual(sequences.read_text(encoding="utf-8").count(">"), expected_records)
        rows = list(csv.DictReader(io.StringIO(metadata.read_text(encoding="utf-8"))))
        self.assertEqual(len(rows), expected_records)
        self.assertEqual(list(rows[0].keys()), METADATA_COLUMNS)
        return rows

    def test_plain_fasta(self):
        source = helpers.write_fasta(self.root / "in.fasta")
        rows = self._assert_outputs(normalise_local(source, self.out, self.work))
        self.assertEqual(rows[0]["accession"], "NC_TEST.1")
        self.assertEqual(rows[1]["accession"], "NC_TEST.2")

    def test_gzipped_fasta(self):
        plain = helpers.write_fasta(self.root / "in.fasta")
        source = self.root / "in.fasta.gz"
        with plain.open("rb") as src, gzip.open(source, "wb") as dst:
            dst.write(src.read())
        rows = self._assert_outputs(normalise_local(source, self.out, self.work))
        self.assertEqual(rows[0]["accession"], "NC_TEST.1")

    def test_zip_of_fasta(self):
        source = helpers.make_zip(
            self.root / "in.zip", {"inner.fasta": helpers.TWO_RECORDS.encode()}
        )
        self._assert_outputs(normalise_local(source, self.out, self.work))

    def test_unsupported_format_is_rejected(self):
        source = self.root / "in.txt"
        source.write_text("nope", encoding="utf-8")
        with self.assertRaises(LoaderError) as ctx:
            normalise_local(source, self.out, self.work)
        self.assertIn("unsupported format", str(ctx.exception))

    def test_missing_input_is_rejected(self):
        with self.assertRaises(LoaderError):
            normalise_local(self.root / "absent.fasta", self.out, self.work)

    def test_zip_without_fasta_is_rejected(self):
        source = helpers.make_zip(self.root / "empty.zip", {"readme.txt": b"hi"})
        with self.assertRaises(LoaderError):
            normalise_local(source, self.out, self.work)

    def test_header_parsing(self):
        source = helpers.write_fasta(self.root / "in.fasta")
        rows = parse_fasta_headers(source)
        self.assertEqual([row["accession"] for row in rows], ["NC_TEST.1", "NC_TEST.2"])


class NcbiReportTests(unittest.TestCase):
    def test_camel_case_record_is_normalised(self):
        record = {
            "accession": "AB123.1",
            "length": 3215,
            "completeness": "COMPLETE",
            "isAnnotated": True,
            "sourceDatabase": "GenBank",
            "releaseDate": "2020-02-02T00:00:00Z",
            "updateDate": "2020-03-03T00:00:00Z",
            "isolate": {"name": "isolate-1"},
            "submitter": {"names": ["Doe J.", "Roe R."]},
            "virus": {
                "organism_name": "Hepatitis B virus",
                "tax_id": 10407,
                "lineage": [
                    {"name": "Hepadnaviridae", "tax_id": 10404},
                    {"name": "Hepatitis B virus", "tax_id": 10407},
                ],
            },
        }
        row = normalise_ncbi_report(record)
        self.assertEqual(row["accession"], "AB123.1")
        self.assertEqual(row["source_database"], "GenBank")
        self.assertEqual(row["release_date"], "2020-02-02T00:00:00Z")
        self.assertEqual(row["organism_name"], "Hepatitis B virus")
        self.assertEqual(row["tax_id"], 10407)
        self.assertEqual(row["isolate"], "isolate-1")
        self.assertEqual(row["lineage"], "Hepadnaviridae; Hepatitis B virus")
        self.assertEqual(row["submitter"], "Doe J.; Roe R.")

    def test_snake_case_record_is_normalised(self):
        record = {
            "accession": "CD456.1",
            "length": 100,
            "completeness": "PARTIAL",
            "source_database": "RefSeq",
            "release_date": "2019-05-05T00:00:00Z",
            "virus": {"organism_name": "Hepatitis B virus", "tax_id": 10407},
        }
        row = normalise_ncbi_report(record)
        self.assertEqual(row["source_database"], "RefSeq")
        self.assertEqual(row["release_date"], "2019-05-05T00:00:00Z")
        self.assertEqual(row["lineage"], "")

    def test_missing_nested_fields_default_to_empty(self):
        row = normalise_ncbi_report({"accession": "EF789.1"})
        self.assertEqual(row["isolate"], "")
        self.assertEqual(row["organism_name"], "")
        self.assertEqual(row["submitter"], "")


class NcbiPackageTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name)
        self.package = helpers.make_ncbi_package(self.root / "package")
        self.out = self.root / "processed"

    def tearDown(self):
        self._tmp.cleanup()

    def test_package_is_normalised(self):
        result = normalise_package(self.package, self.out)
        self.assertEqual(result["records"], 2)
        self.assertEqual(result["source"], "ncbi_datasets_virus")

        rows = list(csv.DictReader(io.StringIO((self.out / "metadata.csv").read_text(encoding="utf-8"))))
        self.assertEqual([row["accession"] for row in rows], ["NC_TEST.1", "NC_TEST.2"])
        self.assertEqual(rows[0]["completeness"], "COMPLETE")
        self.assertEqual(rows[0]["isolate"], "isolate-A")
        self.assertEqual(rows[1]["completeness"], "PARTIAL")
        self.assertEqual(list(rows[0].keys()), METADATA_COLUMNS)

        sequences = (self.out / "sequences.fasta").read_text(encoding="utf-8")
        self.assertEqual(sequences.count(">"), 2)

    def test_missing_genomic_fna_is_rejected(self):
        (self.package / "ncbi_dataset" / "data" / "genomic.fna").unlink()
        with self.assertRaises(LoaderError):
            normalise_package(self.package, self.out)

    def test_missing_report_is_rejected(self):
        (self.package / "ncbi_dataset" / "data" / "data_report.jsonl").unlink()
        with self.assertRaises(LoaderError):
            normalise_package(self.package, self.out)


if __name__ == "__main__":
    unittest.main()
