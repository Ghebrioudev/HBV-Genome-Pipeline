"""Offline unit tests for Task 3 mutation calling and orchestration."""

import csv
import json
import random
import shutil
import tempfile
import unittest
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import analyze_mutations as task3  # noqa: E402


class AlignmentCallingTests(unittest.TestCase):
    def call(self, reference, sample, **kwargs):
        return task3.call_variants_from_alignment(
            "SAMPLE.1", "NC_003977.2", reference, sample,
            min_aligned_bases=kwargs.pop("min_aligned_bases", 1),
            min_identity=kwargs.pop("min_identity", 0.0),
            **kwargs,
        )

    def test_identical_reference_and_sample_has_no_variants(self):
        rows, metrics, reason = self.call("ACGTACGT", "ACGTACGT")
        self.assertEqual(rows, [])
        self.assertIsNone(reason)
        self.assertEqual(metrics.aligned_bases, 8)
        self.assertEqual(metrics.percent_identity, 100.0)

    def test_known_substitution_has_one_based_reference_position(self):
        rows, _, reason = self.call("ACGTACGT", "ACGTTCGT")
        self.assertIsNone(reason)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample_id"], "SAMPLE.1")
        self.assertEqual(rows[0]["reference_accession"], "NC_003977.2")
        self.assertEqual(rows[0]["position"], "5")
        self.assertEqual(rows[0]["reference_base"], "A")
        self.assertEqual(rows[0]["observed_base"], "T")
        self.assertEqual(rows[0]["variation_type"], "substitution")
        self.assertEqual(rows[0]["call_status"], "definitive")

    def test_insertion_uses_between_base_coordinate(self):
        rows, _, reason = self.call("AC-GT", "ACTGT")
        self.assertIsNone(reason)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position"], "2^3")
        self.assertEqual(rows[0]["reference_base"], "-")
        self.assertEqual(rows[0]["observed_base"], "T")
        self.assertEqual(rows[0]["variation_type"], "insertion")

    def test_deletion_is_called_only_inside_observed_coverage(self):
        rows, _, reason = self.call("ACGTAC", "AC-TAC")
        self.assertIsNone(reason)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position"], "3")
        self.assertEqual(rows[0]["reference_base"], "G")
        self.assertEqual(rows[0]["observed_base"], "-")
        self.assertEqual(rows[0]["variation_type"], "deletion")

    def test_partial_terminal_gaps_are_not_reported_as_deletions(self):
        rows, metrics, reason = self.call("AACCGGTT", "--CCGG--", min_aligned_bases=4)
        self.assertIsNone(reason)
        self.assertEqual(rows, [])
        self.assertEqual(metrics.first_reference_position, 3)
        self.assertEqual(metrics.last_reference_position, 6)

    def test_ambiguous_symbol_is_uncertain_not_substitution(self):
        rows, _, reason = self.call("ACGT", "ACNT", min_aligned_bases=3)
        self.assertIsNone(reason)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position"], "3")
        self.assertEqual(rows[0]["variation_type"], "ambiguous")
        self.assertEqual(rows[0]["call_status"], "uncertain")
        self.assertEqual(rows[0]["observed_base"], "N")

    def test_ambiguous_inserted_bases_are_uncertain(self):
        rows, _, reason = self.call("AC-GT", "ACRGT")
        self.assertIsNone(reason)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position"], "2^3")
        self.assertEqual(rows[0]["variation_type"], "ambiguous")
        self.assertEqual(rows[0]["call_status"], "uncertain")

    def test_short_or_low_identity_alignment_is_rejected(self):
        rows, _, reason = self.call(
            "ACGTACGT", "ACGTACGT", min_aligned_bases=9
        )
        self.assertEqual(rows, [])
        self.assertEqual(reason, "insufficient_aligned_bases")
        rows, _, reason = self.call(
            "AAAAAAAA", "CCCCCCCC", min_aligned_bases=8, min_identity=70.0
        )
        self.assertEqual(rows, [])
        self.assertEqual(reason, "identity_below_threshold")

    def test_unequal_alignment_lengths_are_invalid(self):
        with self.assertRaises(ValueError):
            self.call("ACGT", "ACG")


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_fake_mafft(self):
        executable = self.root / "mafft-test"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "if '--version' in sys.argv:\n"
            "    print('v7.526-test-fixture')\n"
            "    raise SystemExit(0)\n"
            "if '--6merpair' not in sys.argv:\n"
            "    print('Expected the low-cost --6merpair strategy', file=sys.stderr)\n"
            "    raise SystemExit(2)\n"
            "fragments = sys.argv[sys.argv.index('--addfragments') + 1]\n"
            "reference = sys.argv[-1]\n"
            "for path in (reference, fragments):\n"
            "    header = None\n"
            "    parts = []\n"
            "    with open(path, encoding='utf-8') as handle:\n"
            "        for line in handle:\n"
            "            line = line.strip()\n"
            "            if not line: continue\n"
            "            if line.startswith('>'):\n"
            "                if header is not None: print(header); print(''.join(parts))\n"
            "                header = line; parts = []\n"
            "            else: parts.append(line)\n"
            "    if header is not None: print(header); print(''.join(parts))\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return str(executable)

    def write_inputs(self, fasta_text, metadata_rows):
        fasta = self.root / "input.fasta"
        metadata = self.root / "metadata.csv"
        fasta.write_text(fasta_text, encoding="utf-8")
        with metadata.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["accession", "completeness"])
            writer.writeheader()
            writer.writerows(metadata_rows)
        return fasta, metadata

    def test_end_to_end_with_mafft_test_double_writes_expected_outputs(self):
        reference = "ACGT" * 30
        sample = reference[:9] + ("T" if reference[9] != "T" else "A") + reference[10:]
        fasta, metadata = self.write_inputs(
            f">NC_003977.2 reference\n{reference}\n>SAMPLE_001.1 test sample\n{sample}\n",
            [
                {"accession": "NC_003977.2", "completeness": "COMPLETE"},
                {"accession": "SAMPLE_001.1", "completeness": "PARTIAL"},
            ],
        )
        output = self.root / "results" / "mutations.csv"
        skipped = self.root / "results" / "skipped.csv"
        report_path = self.root / "results" / "analysis.json"

        report = task3.analyze_dataset(
            input_fasta=fasta,
            metadata_csv=metadata,
            output_csv=output,
            skipped_csv=skipped,
            report_json=report_path,
            mafft=self.write_fake_mafft(),
            min_sequence_length=50,
            min_aligned_bases=50,
            min_identity=70.0,
            batch_size=1,
            threads=1,
        )

        with output.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sample_id"], "SAMPLE_001.1")
        self.assertEqual(rows[0]["position"], "10")
        self.assertEqual(rows[0]["completeness"], "PARTIAL")
        self.assertEqual(report["alignment"]["version"], "v7.526-test-fixture")
        self.assertEqual(report["counts"]["input_sequences"], 2)
        self.assertEqual(report["counts"]["analyzed_sequences"], 1)
        self.assertEqual(report["counts"]["variants_by_type"], {"substitution": 1})
        self.assertTrue(skipped.exists())
        self.assertTrue(report_path.exists())
        self.assertEqual(json.loads(report_path.read_text())["counts"]["mutation_rows"], 1)

    def write_flaky_mafft(self):
        """Fake MAFFT that simulates pairlocalalign crashing on multi-record batches."""
        executable = self.root / "mafft-flaky-test"
        executable.write_text(
            "#!/usr/bin/env python3\n"
            "import sys\n"
            "if '--version' in sys.argv:\n"
            "    print('v7.525-test-fixture')\n"
            "    raise SystemExit(0)\n"
            "if '--6merpair' not in sys.argv:\n"
            "    print('Expected --6merpair', file=sys.stderr)\n"
            "    raise SystemExit(2)\n"
            "fragments = sys.argv[sys.argv.index('--addfragments') + 1]\n"
            "reference = sys.argv[-1]\n"
            "def read(path):\n"
            "    out=[]; header=None; parts=[]\n"
            "    for raw in open(path, encoding='utf-8'):\n"
            "        line=raw.strip()\n"
            "        if not line: continue\n"
            "        if line.startswith('>'):\n"
            "            if header is not None: out.append((header, ''.join(parts)))\n"
            "            header=line; parts=[]\n"
            "        else: parts.append(line)\n"
            "    if header is not None: out.append((header, ''.join(parts)))\n"
            "    return out\n"
            "fragment_records=read(fragments)\n"
            "if len(fragment_records) > 1:\n"
            "    print('simulated pairlocalalign segmentation fault', file=sys.stderr)\n"
            "    raise SystemExit(139)\n"
            "for header, sequence in read(reference) + fragment_records:\n"
            "    print(header); print(sequence)\n",
            encoding="utf-8",
        )
        executable.chmod(0o755)
        return str(executable)

    def test_failed_batch_is_split_and_retried_without_losing_samples(self):
        reference = "ACGT" * 30
        sample_a = reference[:9] + ("T" if reference[9] != "T" else "A") + reference[10:]
        sample_b = reference[:19] + ("C" if reference[19] != "C" else "G") + reference[20:]
        fasta, metadata = self.write_inputs(
            f">NC_003977.2 reference\n{reference}\n"
            f">SAMPLE_A.1 first\n{sample_a}\n"
            f">SAMPLE_B.1 second\n{sample_b}\n",
            [
                {"accession": "NC_003977.2", "completeness": "COMPLETE"},
                {"accession": "SAMPLE_A.1", "completeness": "PARTIAL"},
                {"accession": "SAMPLE_B.1", "completeness": "PARTIAL"},
            ],
        )
        output = self.root / "results" / "mutations.csv"
        skipped = self.root / "results" / "skipped.csv"
        report_path = self.root / "results" / "analysis.json"
        report = task3.analyze_dataset(
            input_fasta=fasta,
            metadata_csv=metadata,
            output_csv=output,
            skipped_csv=skipped,
            report_json=report_path,
            mafft=self.write_flaky_mafft(),
            min_sequence_length=50,
            min_aligned_bases=50,
            min_identity=70.0,
            batch_size=2,
            threads=1,
        )
        with output.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({row["sample_id"] for row in rows}, {"SAMPLE_A.1", "SAMPLE_B.1"})
        self.assertEqual({row["position"] for row in rows}, {"10", "20"})
        self.assertEqual(report["counts"]["analyzed_sequences"], 2)
        self.assertEqual(report["counts"]["skipped_sequences"], 0)
        self.assertGreater(report["counts"]["alignment_batches"], 1)
        self.assertEqual(report["alignment"]["method"], "--6merpair --addfragments")

    def test_unverified_header_only_csv_templates_are_removed(self):
        mutations = self.root / "results" / "mutations.csv"
        skipped = self.root / "results" / "skipped.csv"
        report = self.root / "results" / "analysis.json"
        mutations.parent.mkdir(parents=True)
        mutations.write_text("sample_id,position\n", encoding="utf-8")
        skipped.write_text("sample_id,reason\n", encoding="utf-8")
        task3.remove_unverified_header_only_outputs(mutations, skipped, report)
        self.assertFalse(mutations.exists())
        self.assertFalse(skipped.exists())

    def test_valid_zero_row_csv_is_preserved_when_success_report_exists(self):
        mutations = self.root / "results" / "mutations.csv"
        skipped = self.root / "results" / "skipped.csv"
        report = self.root / "results" / "analysis.json"
        mutations.parent.mkdir(parents=True)
        mutations.write_text("sample_id,position\n", encoding="utf-8")
        skipped.write_text("sample_id,reason\n", encoding="utf-8")
        report.write_text('{"status":"completed"}\n', encoding="utf-8")
        task3.remove_unverified_header_only_outputs(mutations, skipped, report)
        self.assertTrue(mutations.exists())
        self.assertTrue(skipped.exists())

    def test_empty_fasta_fails_without_creating_outputs(self):
        fasta, metadata = self.write_inputs("", [])
        output = self.root / "mutations.csv"
        with self.assertRaisesRegex(ValueError, "not present"):
            task3.analyze_dataset(
                input_fasta=fasta,
                metadata_csv=metadata,
                output_csv=output,
                skipped_csv=self.root / "skipped.csv",
                report_json=self.root / "report.json",
                mafft="missing-mafft-is-not-reached",
            )
        self.assertFalse(output.exists())

    def test_metadata_missing_sample_is_reported_in_skipped_csv(self):
        reference = "ACGT" * 30
        sample = "A" * 120
        fasta, metadata = self.write_inputs(
            f">NC_003977.2 reference\n{reference}\n>NO_META.1 sample\n{sample}\n",
            [{"accession": "NC_003977.2", "completeness": "COMPLETE"}],
        )
        fake_mafft = self.write_fake_mafft()
        output = self.root / "mutations.csv"
        skipped = self.root / "skipped.csv"
        task3.analyze_dataset(
            input_fasta=fasta,
            metadata_csv=metadata,
            output_csv=output,
            skipped_csv=skipped,
            report_json=self.root / "report.json",
            mafft=fake_mafft,
            min_sequence_length=50,
            min_aligned_bases=50,
        )
        with skipped.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["sample_id"], "NO_META.1")
        self.assertEqual(rows[0]["reason"], "missing_metadata")

    def test_invalid_sequence_is_skipped_and_audited(self):
        reference = "ACGT" * 30
        fasta, metadata = self.write_inputs(
            f">NC_003977.2 reference\n{reference}\n>INVALID.1 bad record\n"
            + ("A" * 59) + "?" + ("A" * 60) + "\n",
            [
                {"accession": "NC_003977.2", "completeness": "COMPLETE"},
                {"accession": "INVALID.1", "completeness": "PARTIAL"},
            ],
        )
        skipped = self.root / "skipped.csv"
        task3.analyze_dataset(
            input_fasta=fasta,
            metadata_csv=metadata,
            output_csv=self.root / "mutations.csv",
            skipped_csv=skipped,
            report_json=self.root / "report.json",
            mafft=self.write_fake_mafft(),
            min_sequence_length=50,
            min_aligned_bases=50,
        )
        with skipped.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["sample_id"], "INVALID.1")
        self.assertEqual(rows[0]["reason"], "unsupported_symbols:?")

    def test_missing_mafft_has_actionable_error(self):
        with self.assertRaisesRegex(RuntimeError, "Install MAFFT"):
            task3.get_mafft_version("definitely-not-a-real-mafft-command")


@unittest.skipUnless(shutil.which("mafft"), "real MAFFT executable is not installed")
class RealMafftIntegrationTests(unittest.TestCase):
    """Small optional integration test against the installed MAFFT executable."""

    def test_mafft_addfragments_retains_reference_and_aligns_sample(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generator = random.Random(31415)
            reference_sequence = "".join(generator.choice("ACGT") for _ in range(120))
            replacement = "A" if reference_sequence[39] != "A" else "C"
            sample_sequence = reference_sequence[:39] + replacement + reference_sequence[40:]
            reference_path = root / "reference.fasta"
            fragments_path = root / "fragments.fasta"
            alignment_path = root / "alignment.fasta"
            reference_path.write_text(
                f">NC_003977.2\n{reference_sequence}\n", encoding="utf-8"
            )
            fragments_path.write_text(
                f">SAMPLE.1\n{sample_sequence}\n", encoding="utf-8"
            )
            task3.run_mafft_batch(
                "mafft", reference_path, fragments_path, alignment_path, 1
            )
            aligned = task3.validate_alignment(
                alignment_path,
                "NC_003977.2",
                {"SAMPLE.1"},
                reference_sequence,
            )
            rows, metrics, reason = task3.call_variants_from_alignment(
                "SAMPLE.1",
                "NC_003977.2",
                aligned["NC_003977.2"],
                aligned["SAMPLE.1"],
                min_aligned_bases=100,
                min_identity=70.0,
            )
            self.assertIsNone(reason)
            self.assertEqual(metrics.aligned_bases, 120)
            substitutions = [row for row in rows if row["variation_type"] == "substitution"]
            self.assertEqual(len(substitutions), 1)
            self.assertEqual(substitutions[0]["position"], "40")


if __name__ == "__main__":
    unittest.main()
