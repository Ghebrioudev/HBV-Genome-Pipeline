"""Offline tests for Task 2 preprocessing (Python standard library only)."""

import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import preprocess  # noqa: E402


METADATA_COLUMNS = [
    "accession", "length", "completeness", "is_annotated", "gene_count",
    "protein_count", "source_database", "release_date", "update_date",
    "organism_name", "tax_id", "isolate", "lineage", "submitter",
]
HBV_LINEAGE = (
    "Viruses; Riboviria; Hepadnaviridae; Orthohepadnavirus; Hepatitis B virus"
)
NON_HBV_LINEAGE = "Viruses; Riboviria; Hepadnaviridae; Avihepadnavirus; Duck hepatitis B virus"


def make_dataset(root, records, metadata_rows=None):
    """Create a small Task 1-like dataset and a matching integrity manifest."""
    input_dir = Path(root) / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = input_dir / "sequences.fasta"
    metadata_path = input_dir / "metadata.csv"
    manifest_path = input_dir / "manifest.json"

    with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence in records:
            handle.write(">{}\n".format(header))
            if sequence:
                handle.write(sequence + "\n")

    if metadata_rows is None:
        metadata_rows = []
        for header, sequence in records:
            accession = header.split()[0]
            metadata_rows.append({
                "accession": accession,
                "length": str(len(sequence)),
                "completeness": "PARTIAL" if len(sequence) < 3000 else "COMPLETE",
                "is_annotated": "True",
                "gene_count": "1",
                "protein_count": "1",
                "source_database": "GenBank",
                "release_date": "2025-01-01T00:00:00Z",
                "update_date": "2025-01-01T00:00:00Z",
                "organism_name": "Hepatitis B virus",
                "tax_id": "10407",
                "isolate": "",
                "lineage": HBV_LINEAGE,
                "submitter": "",
            })

    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_COLUMNS)
        writer.writeheader()
        writer.writerows(metadata_rows)

    manifest = {
        "source": "test_fixture",
        "created_at": "2026-01-01T00:00:00+00:00",
        "records": len(records),
        "files": [
            {
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in (fasta_path, metadata_path)
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return input_dir


def run_fixture(tmp_path, records, metadata_rows=None):
    input_dir = make_dataset(tmp_path, records, metadata_rows)
    return preprocess.process_dataset(
        input_dir / "sequences.fasta",
        input_dir / "metadata.csv",
        input_dir / "manifest.json",
    )[0]


class PreprocessingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_valid_sequence_is_retained(self):
        result = run_fixture(self.root, [("AB123456.1 sample", "ACGTACGT")])
        self.assertEqual(result.counts["retained_records"], 1)
        self.assertEqual(result.retained_records[0][0].sequence, "ACGTACGT")

    def test_valid_iupac_ambiguity_is_preserved_and_counted(self):
        result = run_fixture(self.root, [("AB123456.1 sample", "ACGTRYSWKMBDHVN")])
        self.assertEqual(result.counts["retained_records"], 1)
        self.assertEqual(result.retained_records[0][0].sequence, "ACGTRYSWKMBDHVN")
        self.assertEqual(result.retained_records[0][2]["ambiguous_count"], 11)
        self.assertGreater(result.ambiguity_stats["ambiguous_percentage"], 0)

    def test_empty_sequence_is_rejected(self):
        result = run_fixture(self.root, [("AB123456.1 empty", "")])
        self.assertEqual(result.counts["empty_sequence_records"], 1)
        self.assertEqual(result.counts["rejected_fasta_records"], 1)
        self.assertEqual(result.rejected[0]["reason"], "empty_sequence")

    def test_unsupported_characters_are_rejected(self):
        result = run_fixture(self.root, [("AB123456.1 invalid", "ACGTX?")])
        self.assertEqual(result.counts["unsupported_symbol_records"], 1)
        self.assertEqual(result.rejected[0]["symbols"], "?X")
        self.assertIn("unsupported_nucleotide_symbols", result.rejected[0]["reason"])

    def test_uracil_and_gap_are_quarantined_for_review(self):
        result = run_fixture(self.root, [
            ("AB123456.1 uracil", "ACGTU"),
            ("AB123457.1 gap", "ACG-T"),
        ])
        self.assertEqual(result.counts["review_required_records"], 2)
        self.assertEqual(result.counts["records_containing_U"], 1)
        self.assertEqual(result.counts["records_containing_gap"], 1)
        self.assertEqual(result.counts["retained_records"], 0)

    def test_identical_duplicate_accession_retains_one_record(self):
        result = run_fixture(self.root, [
            ("AB123456.1 first", "ACGTACGT"),
            ("AB123456.1 repeated", "acgtacgt"),
        ])
        self.assertEqual(result.counts["retained_records"], 1)
        self.assertEqual(result.counts["duplicate_accession_identical_records_removed"], 1)
        self.assertEqual(result.duplicate_accessions[0]["action"], "deduplicated_identical_copy")

    def test_conflicting_duplicate_accession_quarantines_all_records(self):
        result = run_fixture(self.root, [
            ("AB123456.1 first", "ACGTACGT"),
            ("AB123456.1 conflicting", "ACGTACGA"),
        ])
        self.assertEqual(result.counts["retained_records"], 0)
        self.assertEqual(result.counts["duplicate_accession_conflict_groups"], 1)
        self.assertEqual(result.counts["duplicate_accession_conflict_records"], 2)
        self.assertEqual(len(result.rejected), 2)

    def test_different_accession_versions_are_not_merged(self):
        result = run_fixture(self.root, [
            ("AB123456.1 old", "ACGTACGT"),
            ("AB123456.2 new", "ACGTACGA"),
        ])
        self.assertEqual(result.counts["retained_records"], 2)
        self.assertEqual(result.counts["duplicate_accession_groups"], 0)

    def test_identical_sequences_with_different_accessions_are_preserved(self):
        result = run_fixture(self.root, [
            ("AB123456.1 one", "ACGTACGT"),
            ("CD789012.1 two", "ACGTACGT"),
        ])
        self.assertEqual(result.counts["retained_records"], 2)
        self.assertEqual(result.counts["identical_sequence_groups"], 1)
        self.assertEqual(len(result.identical_sequences), 2)
        self.assertEqual(
            {record[0].accession for record in result.retained_records},
            {"AB123456.1", "CD789012.1"},
        )

    def test_partial_sequences_are_retained_without_minimum_length(self):
        result = run_fixture(self.root, [("AB123456.1 partial", "ACGTAC")])
        self.assertEqual(result.counts["retained_records"], 1)
        self.assertEqual(result.completeness_counts["PARTIAL"], 1)
        self.assertEqual(result.retained_records[0][0].sequence, "ACGTAC")

    def test_taxonomy_filters_non_hbv_and_reports_unresolved(self):
        records = [
            ("AB123456.1 hbv", "ACGT"),
            ("CD123456.1 other", "TGCA"),
            ("EF123456.1 unknown", "AAAA"),
        ]
        metadata = []
        for header, sequence in records:
            accession = header.split()[0]
            if accession.startswith("AB"):
                tax_id, lineage, name = "10407", HBV_LINEAGE, "Hepatitis B virus"
            elif accession.startswith("CD"):
                tax_id, lineage, name = "10406", NON_HBV_LINEAGE, "Duck hepatitis B virus"
            else:
                tax_id, lineage, name = "", "", ""
            metadata.append({
                "accession": accession, "length": str(len(sequence)), "completeness": "PARTIAL",
                "is_annotated": "False", "gene_count": "0", "protein_count": "0",
                "source_database": "GenBank", "release_date": "", "update_date": "",
                "organism_name": name, "tax_id": tax_id, "isolate": "", "lineage": lineage,
                "submitter": "",
            })
        result = run_fixture(self.root, records, metadata)
        self.assertEqual(result.counts["retained_records"], 1)
        self.assertEqual(result.counts["excluded_non_hbv_records"], 1)
        self.assertEqual(result.counts["unresolved_taxonomy_records"], 1)
        self.assertEqual(result.excluded_non_hbv[0]["accession"], "CD123456.1")
        self.assertEqual(result.unresolved_taxonomy[0]["accession"], "EF123456.1")

    def test_hbv_descendant_is_recognized_from_lineage(self):
        metadata = {
            "tax_id": "999999",
            "lineage": HBV_LINEAGE + "; HBV subtype example",
        }
        self.assertEqual(preprocess.taxonomy_classification(metadata)[0], "hbv")

    def test_metadata_is_joined_by_accession_not_row_position(self):
        records = [("AB123456.1 first", "AAAA"), ("CD123456.1 second", "CCCC")]
        base = []
        for accession, seq, isolate in [
            ("AB123456.1", "AAAA", "isolate-A"),
            ("CD123456.1", "CCCC", "isolate-C"),
        ]:
            base.append({
                "accession": accession, "length": str(len(seq)), "completeness": "PARTIAL",
                "is_annotated": "True", "gene_count": "1", "protein_count": "1",
                "source_database": "GenBank", "release_date": "", "update_date": "",
                "organism_name": "Hepatitis B virus", "tax_id": "10407", "isolate": isolate,
                "lineage": HBV_LINEAGE, "submitter": "",
            })
        # Reverse metadata order to prove row position is not used as a join key.
        result = run_fixture(self.root, records, list(reversed(base)))
        isolates = {row[0].accession: row[1]["isolate"] for row in result.retained_records}
        self.assertEqual(isolates["AB123456.1"], "isolate-A")
        self.assertEqual(isolates["CD123456.1"], "isolate-C")

    def test_metadata_length_mismatch_is_sent_to_review(self):
        records = [("AB123456.1 mismatch", "ACGTAC")]
        metadata = [{
            "accession": "AB123456.1", "length": "99", "completeness": "PARTIAL",
            "is_annotated": "True", "gene_count": "1", "protein_count": "1",
            "source_database": "GenBank", "release_date": "", "update_date": "",
            "organism_name": "Hepatitis B virus", "tax_id": "10407", "isolate": "",
            "lineage": HBV_LINEAGE, "submitter": "",
        }]
        result = run_fixture(self.root, records, metadata)
        self.assertEqual(result.counts["review_required_records"], 1)
        self.assertIn("metadata_length_does_not_match", result.review_records[0]["reason"])

    def test_reverse_complements_are_not_grouped_as_identical(self):
        result = run_fixture(self.root, [
            ("AB123456.1 forward", "AACCGTTA"),
            ("CD789012.1 reverse_complement", "AACCGTTA"[::-1].translate(str.maketrans("ACGT", "TGCA"))),
        ])
        self.assertEqual(result.counts["retained_records"], 2)
        self.assertEqual(result.counts["identical_sequence_groups"], 0)

    def test_malformed_fasta_sequence_before_header_is_rejected(self):
        input_dir = make_dataset(self.root, [("AB123456.1 sample", "ACGT")])
        fasta = input_dir / "sequences.fasta"
        fasta.write_text("ACGT\n>AB123456.1 sample\nACGT\n", encoding="utf-8")
        manifest_path = input_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        manifest["files"][0]["bytes"] = fasta.stat().st_size
        manifest["files"][0]["sha256"] = hashlib.sha256(fasta.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaises(preprocess.PreprocessingError):
            preprocess.process_dataset(
                fasta, input_dir / "metadata.csv", manifest_path
            )

    def test_missing_metadata_for_fasta_record_is_unresolved(self):
        records = [("AB123456.1 sample", "ACGT")]
        result = run_fixture(self.root, records, metadata_rows=[])
        self.assertEqual(result.counts["unresolved_taxonomy_records"], 1)
        self.assertEqual(result.unresolved_taxonomy[0]["reason"], "no_matching_metadata_row")

    def test_empty_fasta_dataset_fails_clearly(self):
        input_dir = make_dataset(self.root, [], metadata_rows=[])
        with self.assertRaisesRegex(preprocess.PreprocessingError, "no sequence records"):
            preprocess.process_dataset(
                input_dir / "sequences.fasta", input_dir / "metadata.csv", input_dir / "manifest.json"
            )

    def test_manifest_checksum_mismatch_fails_before_processing(self):
        input_dir = make_dataset(self.root, [("AB123456.1 sample", "ACGT")])
        with (input_dir / "sequences.fasta").open("a", encoding="utf-8") as handle:
            handle.write("\n")
        with self.assertRaises(preprocess.PreprocessingError):
            preprocess.process_dataset(
                input_dir / "sequences.fasta", input_dir / "metadata.csv", input_dir / "manifest.json"
            )

    def test_output_files_and_manifest_are_generated(self):
        input_dir = make_dataset(self.root, [("AB123456.1 sample", "ACGTRYS")])
        result, input_manifest, input_hashes = preprocess.process_dataset(
            input_dir / "sequences.fasta", input_dir / "metadata.csv", input_dir / "manifest.json"
        )
        output_dir = self.root / "data" / "processed"
        reports_dir = self.root / "results" / "preprocessing"
        _metadata, columns, _count = preprocess.load_metadata(input_dir / "metadata.csv")
        manifest = preprocess.write_outputs(
            result, input_manifest, input_hashes, self.root, output_dir, reports_dir, columns
        )
        for name in ("hbv_sequences.fasta", "hbv_metadata.csv", "manifest.json"):
            self.assertTrue((output_dir / name).is_file())
        for name in preprocess.REPORT_COLUMNS:
            self.assertTrue((reports_dir / name).is_file())
        self.assertTrue((reports_dir / "cleaning_report.json").is_file())
        self.assertEqual(manifest["record_counts"]["retained_records"], 1)
        self.assertIn("manifest.json", manifest["input_files"])
        output_fasta = (output_dir / "hbv_sequences.fasta").read_text(encoding="utf-8")
        self.assertIn("ACGTRYS", output_fasta)
        with (output_dir / "hbv_metadata.csv").open(newline="", encoding="utf-8") as handle:
            row = next(csv.DictReader(handle))
        self.assertEqual(row["sequence_length"], "7")
        self.assertEqual(row["contains_ambiguity"], "true")
        manifest_on_disk = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
        self.assertTrue(any(item["path"].endswith("cleaning_report.json") for item in manifest_on_disk["outputs"]))


if __name__ == "__main__":
    unittest.main()
