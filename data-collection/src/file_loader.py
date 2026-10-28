#!/usr/bin/env python3
"""Normalise both input pathways to one downstream contract.

Regardless of where the data came from, this module writes::

    <out-dir>/sequences.fasta   nucleotide sequences (FASTA)
    <out-dir>/metadata.csv      one row per record
    <out-dir>/manifest.json     provenance + counts + checksums

Supported inputs
----------------
* NCBI data package (``--package <raw>/extracted``) — uses ``genomic.fna``
  and the rich ``data_report.jsonl`` metadata.
* Local FASTA: ``.fasta``, ``.fa``, ``.fna``, ``.fas``, optionally ``.gz``.
* Local ZIP archive containing one or more FASTA files.

For local FASTA the metadata CSV is derived from the headers only; for NCBI
packages every available metadata field is preserved.

Usage
-----
    python3 src/file_loader.py --package data/raw/extracted
    python3 src/file_loader.py --input ~/upload/hbv.fasta --size-policy warn
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import shutil
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

CHUNK_SIZE = 1 << 20
FASTA_SUFFIXES = {".fa", ".fna", ".fasta", ".fas"}

METADATA_COLUMNS = [
    "accession",
    "length",
    "completeness",
    "is_annotated",
    "gene_count",
    "protein_count",
    "source_database",
    "release_date",
    "update_date",
    "organism_name",
    "tax_id",
    "isolate",
    "lineage",
    "submitter",
]


class LoaderError(RuntimeError):
    """Raised when the input cannot be read or normalised."""


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in FASTA_SUFFIXES


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _copy_text(source: Path, destination: Path) -> None:
    with _open_text(source) as src, destination.open("w", encoding="utf-8") as dst:
        for line in src:
            dst.write(line)


def _extract_fasta_members(archive_path: Path, workdir: Path) -> list[Path]:
    if not zipfile.is_zipfile(archive_path):
        raise LoaderError(f"{archive_path}: not a valid ZIP archive")
    extracted: list[Path] = []
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            if member.is_dir():
                continue
            member_path = Path(member.filename)
            if not is_fasta(member_path):
                continue
            target = workdir / member_path.name
            with archive.open(member) as src, target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            extracted.append(target)
    return extracted


def prepare_local_fasta(input_path: Path, workdir: Path) -> list[Path]:
    """Return one or more plain FASTA files for a user supplied dataset."""
    if not input_path.exists():
        raise LoaderError(f"{input_path}: does not exist")
    workdir.mkdir(parents=True, exist_ok=True)

    suffix = input_path.suffix.lower()
    if suffix == ".zip":
        files = _extract_fasta_members(input_path, workdir)
        if not files:
            raise LoaderError(f"{input_path}: archive contains no FASTA files")
        return files
    if is_fasta(input_path):
        if input_path.name.lower().endswith(".gz"):
            target = workdir / Path(input_path.name[:-3]).name
            _copy_text(input_path, target)
            return [target]
        return [input_path]
    raise LoaderError(
        f"{input_path}: unsupported format (expected .fasta/.fa/.fna/.fas, optionally .gz or .zip)"
    )


def parse_fasta_headers(path: Path) -> list[dict]:
    """Minimal metadata derived from FASTA headers only."""
    rows: list[dict] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            accession = header.split()[0] if header else ""
            rows.append({"accession": accession, "description": header})
    return rows


def normalise_ncbi_report(row: dict) -> dict:
    """Map an NCBI data report record (camelCase or snake_case) to our columns."""
    virus = row.get("virus") or {}
    isolate = row.get("isolate") or {}
    submitter = row.get("submitter") or {}
    lineage = virus.get("lineage") or []
    if isinstance(lineage, list):
        lineage_text = "; ".join(item.get("name", "") for item in lineage if isinstance(item, dict))
    else:
        lineage_text = str(lineage)
    if isinstance(isolate, dict):
        isolate_text = isolate.get("name") or isolate.get("isolate_name") or ""
    else:
        isolate_text = str(isolate)
    if isinstance(submitter, dict):
        names = submitter.get("names") or submitter.get("name") or ""
        submitter_text = "; ".join(names) if isinstance(names, list) else str(names)
    else:
        submitter_text = str(submitter)

    return {
        "accession": row.get("accession", ""),
        "length": row.get("length", ""),
        "completeness": row.get("completeness", ""),
        "is_annotated": row.get("is_annotated", row.get("isAnnotated", "")),
        "gene_count": row.get("gene_count", row.get("geneCount", "")),
        "protein_count": row.get("protein_count", row.get("proteinCount", "")),
        "source_database": row.get("source_database", row.get("sourceDatabase", "")),
        "release_date": row.get("release_date", row.get("releaseDate", "")),
        "update_date": row.get("update_date", row.get("updateDate", "")),
        "organism_name": virus.get("organism_name", virus.get("organismName", "")),
        "tax_id": virus.get("tax_id", virus.get("taxId", "")),
        "isolate": isolate_text,
        "lineage": lineage_text,
        "submitter": submitter_text,
    }


def read_ncbi_report(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            rows.append(normalise_ncbi_report(json.loads(line)))
    return rows


def write_metadata(rows: list[dict], destination: Path) -> None:
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=METADATA_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in METADATA_COLUMNS})


def concatenate_fasta(files: list[Path], destination: Path) -> int:
    records = 0
    with destination.open("w", encoding="utf-8") as out:
        for file in files:
            with _open_text(file) as handle:
                for line in handle:
                    if line.startswith(">"):
                        records += 1
                    out.write(line if line.endswith("\n") else line + "\n")
    return records


def _write_manifest(out_dir: Path, source: str, files: list[Path], records: int) -> Path:
    manifest = {
        "source": source,
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "records": records,
        "files": [
            {
                "name": file.name,
                "bytes": file.stat().st_size,
                "sha256": sha256_file(file),
            }
            for file in files
            if file.exists()
        ],
    }
    path = out_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return path


def normalise_package(package_dir: Path, out_dir: Path) -> dict:
    fasta = package_dir / "ncbi_dataset" / "data" / "genomic.fna"
    report = package_dir / "ncbi_dataset" / "data" / "data_report.jsonl"
    if not fasta.exists():
        raise LoaderError(f"{fasta}: missing from NCBI package")
    if not report.exists():
        raise LoaderError(f"{report}: missing from NCBI package")

    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "sequences.fasta"
    shutil.copyfile(fasta, target)
    rows = read_ncbi_report(report)
    write_metadata(rows, out_dir / "metadata.csv")
    records = len(rows) or concatenate_fasta([fasta], out_dir / "sequences.fasta")
    path = _write_manifest(out_dir, "ncbi_datasets_virus", [target, out_dir / "metadata.csv"], records)
    return {"source": "ncbi_datasets_virus", "records": records, "manifest": str(path)}


def normalise_local(input_path: Path, out_dir: Path, workdir: Path) -> dict:
    files = prepare_local_fasta(input_path, workdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "sequences.fasta"
    records = concatenate_fasta(files, target)

    rows: list[dict] = []
    for file in files:
        rows.extend(parse_fasta_headers(file))
    write_metadata(rows, out_dir / "metadata.csv")
    path = _write_manifest(out_dir, f"local_file:{input_path.name}", [target, out_dir / "metadata.csv"], records)
    return {"source": f"local_file:{input_path.name}", "records": records, "manifest": str(path)}


def main(argv: list[str] | None = None) -> int:
    base_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--package", help="extracted NCBI data package directory")
    group.add_argument("--input", help="local dataset file (.fasta/.fa/.fna, .gz or .zip)")
    parser.add_argument("--out-dir", default=str(base_dir / "data" / "processed"), help="normalised output directory")
    parser.add_argument("--workdir", default=str(base_dir / "data" / "raw" / "local"), help="scratch directory for archives")
    args = parser.parse_args(argv)

    out_dir = Path(args.out_dir).expanduser().resolve()
    try:
        if args.package:
            result = normalise_package(Path(args.package).expanduser().resolve(), out_dir)
        else:
            result = normalise_local(
                Path(args.input).expanduser().resolve(),
                out_dir,
                Path(args.workdir).expanduser().resolve(),
            )
    except (LoaderError, json.JSONDecodeError, OSError) as exc:
        log(f"[loader] ERROR: {exc}")
        return 1

    log(f"[loader] {result['records']} records -> {out_dir}")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
