#!/usr/bin/env python3
"""Validate a downloaded or user-supplied nucleotide dataset.

Checks performed
----------------
1. The path exists and is non-empty (file or directory).
2. ZIP archives are readable, complete and contain the expected members.
3. Every FASTA file parses: at least one record, a header per record, and only
   IUPAC nucleotide characters in the sequence bodies.
4. The dataset size is compared against the project window (default 100 MB –
   1 GB). The policy decides whether an out-of-window size is fatal (``hard``)
   or only a warning (``warn``).
5. SHA-256 checksums are calculated for every artefact and written out.

Usage
-----
    python3 src/dataset_validator.py --dataset data/raw/extracted --size-policy hard
    python3 src/dataset_validator.py --dataset ~/upload/hbv.fasta --size-policy warn
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path

CHUNK_SIZE = 1 << 20
FASTA_SUFFIXES = {".fa", ".fna", ".fasta", ".fas"}
IUPAC = set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv-")
MIN_MB = 100.0
MAX_MB = 1024.0


class ValidationError(RuntimeError):
    """Raised when a dataset fails a hard validation rule."""


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def path_size(path: Path) -> int:
    """Total bytes of a file, or of all regular files below a directory."""
    if path.is_file():
        return path.stat().st_size
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def is_fasta_name(name: str) -> bool:
    """True for .fa/.fna/.fasta/.fas, optionally followed by .gz."""
    lowered = name.lower()
    if lowered.endswith(".gz"):
        lowered = lowered[:-3]
    return Path(lowered).suffix in FASTA_SUFFIXES


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def find_fasta_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path] if is_fasta_name(path.name) else []
    return sorted(item for item in path.rglob("*") if item.is_file() and is_fasta_name(item.name))


def count_fasta_stream(lines, label: str) -> dict:
    """Count records and nucleotides in an iterable of FASTA text lines."""
    records = 0
    total_nt = 0
    invalid_chars = 0
    empty_records = 0
    current_len = 0
    longest = 0

    for line in lines:
        if line.startswith(">"):
            if records and current_len == 0:
                empty_records += 1
            records += 1
            total_nt += current_len
            longest = max(longest, current_len)
            current_len = 0
            continue
        sequence = line.strip()
        if not sequence:
            continue
        current_len += len(sequence)
        invalid_chars += sum(1 for char in sequence if char not in IUPAC)

    if records and current_len == 0:
        empty_records += 1
    total_nt += current_len
    longest = max(longest, current_len)

    if records == 0:
        raise ValidationError(f"{label}: no FASTA records found")
    if empty_records:
        raise ValidationError(f"{label}: {empty_records} record(s) have an empty sequence")
    if invalid_chars:
        raise ValidationError(f"{label}: {invalid_chars} non-IUPAC character(s) in sequence data")

    return {"file": label, "records": records, "total_nt": total_nt, "longest_record": longest}


def count_fasta(path: Path) -> dict:
    """Stream a (optionally gzipped) FASTA file and return record statistics."""
    with _open_text(path) as handle:
        result = count_fasta_stream(handle, str(path))
    result["sha256"] = sha256_file(path)
    return result


def count_fasta_in_zip(path: Path) -> list[dict]:
    """Validate and count every FASTA member inside a ZIP archive."""
    results: list[dict] = []
    with zipfile.ZipFile(path) as archive:
        members = [m for m in archive.infolist() if not m.is_dir() and is_fasta_name(m.filename)]
        if not members:
            raise ValidationError(f"{path}: archive contains no FASTA files")
        for member in members:
            with archive.open(member) as raw:
                if member.filename.lower().endswith(".gz"):
                    raw = gzip.open(raw, "rt", encoding="utf-8", errors="replace")
                else:
                    raw = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
                results.append(count_fasta_stream(raw, f"{path}!{member.filename}"))
    return results


def validate_zip(path: Path) -> dict:
    if not zipfile.is_zipfile(path):
        raise ValidationError(f"{path}: not a valid ZIP archive")
    with zipfile.ZipFile(path) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise ValidationError(f"{path}: corrupt member {bad!r}")
        members = archive.namelist()
    if not members:
        raise ValidationError(f"{path}: archive is empty")
    return {"file": str(path), "members": len(members), "sha256": sha256_file(path)}


def validate_dataset(
    path: Path,
    min_mb: float = MIN_MB,
    max_mb: float = MAX_MB,
    size_policy: str = "hard",
) -> dict:
    if size_policy not in {"hard", "warn"}:
        raise ValueError("size_policy must be 'hard' or 'warn'")
    if not path.exists():
        raise ValidationError(f"{path}: does not exist")

    report: dict = {"path": str(path), "size_policy": size_policy, "warnings": [], "fasta": []}
    size_bytes = path_size(path)
    size_mb = size_bytes / (1024 * 1024)
    report["size_bytes"] = size_bytes
    report["size_mb"] = round(size_mb, 2)

    if size_bytes == 0:
        raise ValidationError(f"{path}: dataset is empty")

    if path.is_file() and path.suffix.lower() == ".zip":
        report["zip"] = validate_zip(path)
        report["fasta"].extend(count_fasta_in_zip(path))
    else:
        fasta_files = find_fasta_files(path)
        if not fasta_files:
            raise ValidationError(f"{path}: no FASTA files found (.fa/.fna/.fasta, optionally .gz)")

        for fasta in fasta_files:
            report["fasta"].append(count_fasta(fasta))

    report["records"] = sum(item["records"] for item in report["fasta"])
    report["total_nt"] = sum(item["total_nt"] for item in report["fasta"])

    if size_mb < min_mb:
        message = (f"dataset is {size_mb:.1f} MB, below the required minimum of {min_mb:.0f} MB")
        if size_policy == "hard":
            raise ValidationError(message)
        report["warnings"].append(message)
    if size_mb > max_mb:
        message = (f"dataset is {size_mb:.1f} MB, above the allowed maximum of {max_mb:.0f} MB")
        if size_policy == "hard":
            raise ValidationError(message)
        report["warnings"].append(message)

    report["ok"] = True
    return report


def write_checksums(files: list[Path], destination: Path, base: Path | None = None) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as handle:
        for file in files:
            if file.is_file():
                shown = os.path.relpath(file, base) if base else str(file)
                handle.write(f"{sha256_file(file)}  {shown}\n")
    return destination


def main(argv: list[str] | None = None) -> int:
    base_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", required=True, help="dataset file or directory to validate")
    parser.add_argument("--min-mb", type=float, default=MIN_MB, help=f"minimum size in MB (default {MIN_MB:.0f})")
    parser.add_argument("--max-mb", type=float, default=MAX_MB, help=f"maximum size in MB (default {MAX_MB:.0f})")
    parser.add_argument("--size-policy", choices=("hard", "warn"), default="hard",
                        help="fail (hard) or only warn (warn) when the size is out of range")
    parser.add_argument("--report-json", default="", help="write the validation report to this path")
    parser.add_argument("--checksums", default="", help="write a sha256sums file to this path")
    parser.add_argument("--checksums-base", default="",
                        help="write paths in the checksum file relative to this directory")
    args = parser.parse_args(argv)

    dataset = Path(args.dataset).expanduser().resolve()
    try:
        report = validate_dataset(dataset, args.min_mb, args.max_mb, args.size_policy)
    except (ValidationError, ValueError, OSError) as exc:
        log(f"[validate] FAILED: {exc}")
        return 1

    report["validated_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    if args.checksums:
        files = [dataset] if dataset.is_file() else sorted(
            item for item in dataset.rglob("*") if item.is_file()
        )
        path = write_checksums(
            files,
            Path(args.checksums),
            Path(args.checksums_base).expanduser().resolve() if args.checksums_base else None,
        )
        report["checksums_file"] = str(path)

    if args.report_json:
        Path(args.report_json).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(report, indent=2))
    for warning in report.get("warnings", []):
        log(f"[validate] WARNING: {warning}")
    log(f"[validate] OK: {report['records']} records, {report['total_nt']} nt, {report['size_mb']} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
