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
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ("run.sh", "download.sh", "preprocess.py", "analyze_mutations.py", "statistics.py",
           "report/make_dashboard.py", "scripts/lib.sh", "scripts/shell_analysis.sh")
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
    import os, time
    time.sleep(float(os.environ.get("HBV_STUB_SLEEP", "0")))
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


def remove_tree(path):
    """rmtree that first restores write permission (raw data are read-only)."""
    for root, dirs, files in os.walk(path):
        for name in dirs + files:
            try:
                os.chmod(os.path.join(root, name), 0o755)
            except OSError:
                pass
    shutil.rmtree(path, ignore_errors=True)


def copy_scripts(project):
    for name in SCRIPTS:
        (project / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / name, project / name)
    (project / "run.sh").chmod(0o755)


def build_fixture(project):
    """Create the Task 1 output contract and a MAFFT stub for the project."""
    copy_scripts(project)

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

    # Raw NCBI-like package with its md5sum.txt, verified by run.sh then locked.
    raw = project / "data-collection" / "data" / "raw" / "extracted"
    (raw / "ncbi_dataset" / "data").mkdir(parents=True)
    genomic = raw / "ncbi_dataset" / "data" / "genomic.fna"
    shutil.copy2(fasta, genomic)
    md5 = hashlib.md5(genomic.read_bytes()).hexdigest()
    (raw / "md5sum.txt").write_text(f"{md5}  ncbi_dataset/data/genomic.fna\n", encoding="utf-8")

    # Task 1 checksum list (paths relative to data-collection/), verified by run.sh.
    checksums = project / "data-collection" / "data" / "checksums"
    checksums.mkdir(parents=True)
    (checksums / "dataset.sha256").write_text(
        "# integration fixture\n" + "".join(
            f"{sha256(path)}  data/processed/{path.name}\n" for path in (manifest, metadata, fasta)
        ),
        encoding="utf-8",
    )

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
        self.addCleanup(remove_tree, self.project)
        self.bindir = build_fixture(self.project)

    def run_pipeline(self, *arguments, project=None, bindir=None, extra_env=None):
        project = project or self.project
        env = dict(os.environ)
        env.update(extra_env or {})
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
        copy_scripts(empty)
        result = self.run_pipeline("--skip-download", project=empty, bindir=self.bindir)
        self.assertNotEqual(result.returncode, 0)
        self.assertTrue(
            "does not exist" in result.stderr.lower()
            or "no such file" in result.stderr.lower(),
            result.stderr,
        )

    # ---------------------------------------------------------- run.sh contract
    def test_help_lists_getopts_options(self):
        result = self.run_pipeline("-h")
        self.assertEqual(result.returncode, 0, result.stderr)
        for option in ("-h", "-s", "-c", "-f", "-n", "-t TAXON", "-j THREADS"):
            self.assertIn(option, result.stdout)

    def test_option_value_is_validated(self):
        missing_value = self.run_pipeline("-j")
        self.assertEqual(missing_value.returncode, 2)
        self.assertIn("requires a value", missing_value.stderr)
        bad_value = self.run_pipeline("-s", "-j", "zero")
        self.assertEqual(bad_value.returncode, 2)
        self.assertIn("positive integer", bad_value.stderr)

    def test_full_run_writes_log_checksums_archive_and_locks_raw_data(self):
        result = self.run_pipeline("-s", "-k", "5")
        self.assertEqual(result.returncode, 0, result.stderr)

        logs = sorted((self.project / "logs").glob("run_*.log"))
        self.assertEqual(len(logs), 1)
        first_line = logs[0].read_text(encoding="utf-8").splitlines()[0]
        self.assertRegex(first_line, r"^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\]")

        archives = list((self.project / "archives").glob("hbv_results_*.tar.gz"))
        self.assertEqual(len(archives), 1)
        self.assertRegex(archives[0].name, r"^hbv_results_\d{8}_\d{6}\.tar\.gz$")
        listing = subprocess.run(["tar", "-tzf", str(archives[0])], text=True,
                                 capture_output=True, check=True).stdout
        self.assertIn("results/summary.csv", listing)
        self.assertIn("report/dashboard.html", listing)

        checksums = (self.project / "results" / "checksums.sha256").read_text()
        self.assertIn("results/summary.csv", checksums)
        self.assertTrue((self.project / "report" / "dashboard.html").is_file())

        for name in ("sequence_profile.tsv", "length_distribution.tsv",
                     "crosscheck_top_mutations.tsv", "events_by_type.tsv"):
            self.assertTrue((self.project / "results" / "shell" / name).is_file(), name)
        crosscheck = (self.project / "results" / "shell" / "crosscheck_top_mutations.tsv").read_text()
        self.assertIn("61:A>T\t1\t1\tyes", crosscheck)

        genomic = (self.project / "data-collection" / "data" / "raw" / "extracted"
                   / "ncbi_dataset" / "data" / "genomic.fna")
        self.assertEqual(genomic.stat().st_mode & 0o222, 0, "raw data must be read-only")

    def test_rerun_without_manual_cleanup(self):
        first = self.run_pipeline("-s", "-n")
        self.assertEqual(first.returncode, 0, first.stderr)
        second = self.run_pipeline("-s", "-n")
        self.assertEqual(second.returncode, 0, second.stderr)
        rows = (self.project / "results" / "mutations.csv").read_text().splitlines()
        self.assertEqual(len(rows), 2)

    def test_checksum_mismatch_stops_the_pipeline(self):
        fasta = self.project / "data-collection" / "data" / "processed" / "sequences.fasta"
        with fasta.open("a", encoding="utf-8") as handle:
            handle.write(">TAMPERED.1\nACGT\n")
        result = self.run_pipeline("-s")
        self.assertEqual(result.returncode, 5, result.stderr)
        self.assertIn("checksum verification failed", result.stderr)
        self.assertFalse((self.project / "results" / "mutations.csv").exists())

    def test_shell_top_n_on_large_input_does_not_sigpipe(self):
        """Regression: `sort | head` under pipefail died with SIGPIPE (exit 141)."""
        work = self.project / "big"
        work.mkdir()
        mutations = work / "mutations.csv"
        summary = work / "summary.csv"
        header = ("sample_id,reference_accession,position,reference_base,observed_base,"
                  "variation_type,call_status,aligned_bases,percent_identity,completeness\n")
        with mutations.open("w", encoding="utf-8") as handle:
            handle.write(header)
            # 20,000 distinct mutations: the sorted stream (~300 KB) is larger
            # than a pipe buffer (64 KB), as on the real 10-million-row table.
            for i in range(100_000):
                handle.write(f"S{i}.1,NC_003977.2,{1 + i % 20_000},A,T,substitution,definitive,3182,99.0,COMPLETE\n")
        summary.write_text("mutation,position,ref,alt,type,n_sequences,frequency\n", encoding="utf-8")
        result = subprocess.run(
            ["bash", str(self.project / "scripts" / "shell_analysis.sh"), "-m", "mutations",
             "-i", str(mutations), "-s", str(summary), "-o", str(work / "out"), "-n", "3"],
            text=True, capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        top = (work / "out" / "top_mutations_shell.tsv").read_text().splitlines()
        self.assertEqual(len(top), 4)  # header + 3 rows

    def test_sigterm_removes_temporary_files_then_rerun_succeeds(self):
        if not shutil.which("pkill"):
            self.skipTest("pkill not available")
        tmp_root = Path(tempfile.mkdtemp(prefix="hbv_tmproot_"))
        self.addCleanup(shutil.rmtree, tmp_root, ignore_errors=True)
        env = dict(os.environ, PYTHON=sys.executable, TMPDIR=str(tmp_root),
                   HBV_STUB_SLEEP="60", PATH=f"{self.bindir}:{os.environ['PATH']}")
        process = subprocess.Popen([str(self.project / "run.sh"), "-s", "-n"], cwd=self.project,
                                   env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        deadline = time.time() + 60
        while time.time() < deadline:
            logs = list((self.project / "logs").glob("run_*.log"))
            if logs and "START  5/10" in logs[0].read_text(encoding="utf-8"):
                break
            time.sleep(0.2)
        else:
            process.kill()
            self.fail("pipeline never reached the alignment stage")
        time.sleep(1)
        self.assertTrue(any(tmp_root.iterdir()), "run.sh should have created a temp dir")
        process.send_signal(signal.SIGTERM)
        _, stderr = process.communicate(timeout=30)
        self.assertEqual(process.returncode, 143, stderr)
        self.assertIn("SIGTERM", stderr)
        self.assertEqual(list(tmp_root.iterdir()), [], "temporary files were not removed")
        self.assertFalse((self.project / "logs" / ".run.pid").exists())

        rerun = self.run_pipeline("-s", "-n")
        self.assertEqual(rerun.returncode, 0, rerun.stderr)


if __name__ == "__main__":
    unittest.main()
