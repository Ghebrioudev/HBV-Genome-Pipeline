"""Clean-start integration test: run all four pipeline stages on a tiny fixture.

Unlike the orchestration-only check it replaces, this test executes the real
`preprocess.py`, `analyze_mutations.py` and `statistics.py`. MAFFT is replaced by
a deterministic stub that emits a valid right-padded alignment, so the test is
offline and does not need the real binary.

The fixture is a 240 nt synthetic reference plus three samples: one identical,
one with a single substitution (position 61, A>T) and one partial (180 nt).
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ("run.sh", "download.sh", "preprocess.py", "analyze_mutations.py", "statistics.py")
REFERENCE = "ACGT" * 60
SAMPLES = {
    "S1.1": REFERENCE,
    "S2.1": REFERENCE[:60] + "T" + REFERENCE[61:],   # 61 A>T
    "S3.1": REFERENCE[:180],                          # partial
}
COMPLETENESS = {"S1.1": "COMPLETE", "S2.1": "COMPLETE", "S3.1": "PARTIAL"}
METADATA_HEADER = [
    "accession", "length", "completeness", "is_annotated", "gene_count",
    "protein_count", "source_database", "release_date", "update_date",
    "organism_name", "tax_id", "isolate", "lineage", "submitter",
]

MAFFT_STUB = r'''#!/usr/bin/env python3
"""Deterministic MAFFT stand-in: right-pad every record to one width."""
import sys


def read_fasta(path):
    records, name, buf = [], None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        records.append((name, "".join(buf)))
    return records


def main():
    args = sys.argv[1:]
    if "--version" in args:
        print("v7.999 (integration stub)")
        return 0
    fragments = args[args.index("--addfragments") + 1]
    reference = args[-1]
    records = read_fasta(reference) + read_fasta(fragments)
    width = max(len(sequence) for _, sequence in records)
    for name, sequence in records:
        print(">" + name)
        print(sequence + "-" * (width - len(sequence)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_fixture(project):
    """Create the Task 1 output contract and a MAFFT stub for the project."""
    for name in SCRIPTS:
        shutil.copy2(ROOT / name, project / name)
    (project / "run.sh").chmod(0o755)

    processed = project / "data-collection" / "data" / "processed"
    processed.mkdir(parents=True)

    fasta = processed / "sequences.fasta"
    fasta.write_text(
        ">NC_003977.2 Hepatitis B virus reference (synthetic)\n"
        + REFERENCE + "\n"
        + "".join(f">{name} synthetic sample\n{seq}\n" for name, seq in SAMPLES.items()),
        encoding="utf-8",
    )

    metadata = processed / "metadata.csv"
    rows = [",".join(METADATA_HEADER)]
    lineage = (
        "Viruses; Riboviria; Pararnavirae; Artverviricota; Revtraviricetes; "
        "Blubervirales; Hepadnaviridae; Orthohepadnavirus; Hepatitis B virus"
    )
    records = [("NC_003977.2", REFERENCE, "COMPLETE")] + [
        (name, sequence, COMPLETENESS[name]) for name, sequence in SAMPLES.items()
    ]
    for name, sequence, completeness in records:
        rows.append(",".join([
            name, str(len(sequence)), completeness, "True", "4", "7",
            "RefSeq", "2026-01-01T00:00:00Z", "2026-01-01T00:00:00Z",
            "Hepatitis B virus", "10407", "synthetic", lineage, "integration-test",
        ]))
    metadata.write_text("\n".join(rows) + "\n", encoding="utf-8")

    manifest = processed / "manifest.json"
    manifest.write_text(json.dumps({
        "source": "integration-fixture",
        "records": 1 + len(SAMPLES),
        "files": [
            {"name": fasta.name, "bytes": fasta.stat().st_size, "sha256": sha256(fasta)},
            {"name": metadata.name, "bytes": metadata.stat().st_size, "sha256": sha256(metadata)},
        ],
    }), encoding="utf-8")

    bindir = project / "bin"
    bindir.mkdir()
    stub = bindir / "mafft"
    stub.write_text(MAFFT_STUB, encoding="utf-8")
    stub.chmod(0o755)
    return bindir


class PipelineIntegrationTest(unittest.TestCase):
    def setUp(self):
        try:
            import matplotlib  # noqa: F401
            import pandas  # noqa: F401
        except ImportError as exc:
            raise unittest.SkipTest(f"Task 4 dependencies unavailable: {exc}")
        self.project = Path(tempfile.mkdtemp(prefix="hbv_integration_"))
        self.addCleanup(shutil.rmtree, self.project, ignore_errors=True)
        self.bindir = build_fixture(self.project)

    def run_pipeline(self, *arguments, project=None, bindir=None):
        project = project or self.project
        env = dict(os.environ)
        env["PYTHON"] = sys.executable
        env["PATH"] = f"{bindir or self.bindir}:{env['PATH']}"
        return subprocess.run(
            [str(project / "run.sh"), *arguments],
            cwd=project, env=env, text=True, capture_output=True,
        )

    def test_clean_start_runs_all_four_stages(self):
        result = self.run_pipeline("--skip-download")
        self.assertEqual(
            result.returncode, 0,
            f"pipeline failed\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )

        processed = self.project / "data" / "processed"
        self.assertTrue((processed / "hbv_sequences.fasta").is_file())
        self.assertTrue((processed / "hbv_metadata.csv").is_file())

        mutations = self.project / "results" / "mutations.csv"
        lines = mutations.read_text(encoding="utf-8").splitlines()
        self.assertEqual(
            lines[0],
            "sample_id,reference_accession,position,reference_base,"
            "observed_base,variation_type,call_status,aligned_bases,"
            "percent_identity,completeness",
        )
        rows = lines[1:]
        self.assertEqual(len(rows), 1, f"expected one substitution, got: {rows}")
        self.assertEqual(
            rows[0].split(",")[:7],
            ["S2.1", "NC_003977.2", "61", "A", "T", "substitution", "definitive"],
        )

        for name in ("summary.csv", "mutation_types.csv", "statistics_summary.json"):
            self.assertTrue((self.project / "results" / name).is_file(), name)
        for name in ("top_mutations.png", "mutations_along_genome.png", "mutation_types.png"):
            self.assertTrue((self.project / "results" / "figures" / name).is_file(), name)

        summary = json.loads((self.project / "results" / "statistics_summary.json").read_text())
        self.assertEqual(summary["denominator"], 3)
        self.assertEqual(summary["n_distinct_mutations"], 1)

    def test_unknown_option_is_rejected(self):
        result = self.run_pipeline("--not-a-real-option")
        self.assertEqual(result.returncode, 2)
        self.assertIn("unknown option", result.stderr.lower())

    def test_missing_input_fails_loudly(self):
        empty = Path(tempfile.mkdtemp(prefix="hbv_empty_"))
        self.addCleanup(shutil.rmtree, empty, ignore_errors=True)
        for name in SCRIPTS:
            shutil.copy2(ROOT / name, empty / name)
        (empty / "run.sh").chmod(0o755)
        result = self.run_pipeline("--skip-download", project=empty, bindir=self.bindir)
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(
            "does not exist" in result.stderr.lower()
            or "no such file" in result.stderr.lower(),
            result.stderr,
        )


if __name__ == "__main__":
    unittest.main()
