"""Tests for Task 4 (statistics.py) using tiny artificial data."""
import importlib.util
import tempfile
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# Loaded under another name so it never clashes with the stdlib 'statistics'.
spec = importlib.util.spec_from_file_location("task4_statistics", ROOT / "statistics.py")
stats = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stats)

MUTATIONS = pd.DataFrame(
    {
        "sequence_accession": ["S1", "S2", "S3", "S1", "S2"],
        "reference_accession": ["NC_003977.2"] * 5,
        "position": [100, 100, 100, 200, 300],
        "reference_nucleotide": ["A", "A", "A", "C", "G"],
        "sample_nucleotide": ["G", "G", "G", "T", "R"],  # R = ambiguous
        "mutation_type": ["substitution"] * 5,
    }
)


class TestStatistics(unittest.TestCase):
    def test_ambiguous_calls_are_separated(self):
        confirmed, uncertain = stats.split_uncertain(stats.standardize(MUTATIONS))
        self.assertEqual(len(confirmed), 4)
        self.assertEqual(len(uncertain), 1)

    def test_frequency_uses_given_denominator(self):
        confirmed, _ = stats.split_uncertain(stats.standardize(MUTATIONS))
        freq = stats.compute_frequencies(confirmed, denominator=10)
        top = freq.iloc[0]
        self.assertEqual(top["mutation"], "100:A>G")
        self.assertEqual(top["n_sequences"], 3)
        self.assertAlmostEqual(top["frequency"], 0.3)

    def test_missing_column_raises_clear_error(self):
        with self.assertRaises(ValueError):
            stats.standardize(MUTATIONS.drop(columns=["position"]))

    def test_end_to_end_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            MUTATIONS.to_csv(tmp / "mutations.csv", index=False)
            code = stats.main(
                ["--mutations", str(tmp / "mutations.csv"),
                 "--metadata", str(tmp / "none.csv"),
                 "--outdir", str(tmp / "out"), "--total-sequences", "10"]
            )
            self.assertEqual(code, 0)
            self.assertTrue((tmp / "out" / "summary.csv").is_file())
            self.assertTrue((tmp / "out" / "statistics_summary.json").is_file())
            self.assertTrue((tmp / "out" / "figures" / "top_mutations.png").is_file())


if __name__ == "__main__":
    unittest.main()
