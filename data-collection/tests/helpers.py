"""Shared helpers for the Phase 1 test suite.

Kept dependency-free: only the Python standard library is used, matching the
project rule that the pipeline itself must not require third-party packages.
"""

from __future__ import annotations

import json
import sys
import zipfile
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parents[1]
SRC_DIR = BASE_DIR / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

SAMPLE_FASTA = BASE_DIR / "tests" / "data" / "sample.fasta"
SAMPLE_GZ = BASE_DIR / "tests" / "data" / "sample.fasta.gz"
SAMPLE_ZIP = BASE_DIR / "tests" / "data" / "sample.zip"

TWO_RECORDS = (
    ">NC_TEST.1 Hepatitis B virus isolate test\n"
    "ATGCGTACGTACGTAGCTAGCTAGCTAGCNNNN\n"
    ">NC_TEST.2 Hepatitis B virus isolate test 2\n"
    "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n"
)


def write_fasta(path: Path, content: str = TWO_RECORDS) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def make_zip(path: Path, members: dict[str, bytes]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in members.items():
            archive.writestr(name, data)
    return path


def make_ncbi_package(root: Path, records: list[dict] | None = None) -> Path:
    """Build a minimal NCBI-shaped package: ncbi_dataset/data/{genomic.fna,data_report.jsonl}."""
    data_dir = root / "ncbi_dataset" / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "genomic.fna").write_text(TWO_RECORDS, encoding="utf-8")

    if records is None:
        records = [
            {
                "accession": "NC_TEST.1",
                "length": 34,
                "completeness": "COMPLETE",
                "isAnnotated": True,
                "sourceDatabase": "GenBank",
                "releaseDate": "2020-01-01T00:00:00Z",
                "updateDate": "2021-01-01T00:00:00Z",
                "isolate": {"name": "isolate-A"},
                "virus": {
                    "organism_name": "Hepatitis B virus",
                    "tax_id": 10407,
                    "lineage": [
                        {"name": "Viruses", "tax_id": 10239},
                        {"name": "Hepadnaviridae", "tax_id": 10404},
                        {"name": "Hepatitis B virus", "tax_id": 10407},
                    ],
                },
                "submitter": {"names": ["Doe J."]},
            },
            {
                "accession": "NC_TEST.2",
                "length": 32,
                "completeness": "PARTIAL",
                "source_database": "RefSeq",
                "release_date": "2019-05-05T00:00:00Z",
                "virus": {"organism_name": "Hepatitis B virus", "tax_id": 10407},
            },
        ]

    with (data_dir / "data_report.jsonl").open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")
    return root
