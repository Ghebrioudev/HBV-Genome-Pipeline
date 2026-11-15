#!/usr/bin/env python3
"""Download virus genome data packages from the official NCBI Datasets API.

Endpoint
--------
POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download

This is the REST equivalent of::

    datasets download virus genome taxon <taxon> --include genome

The API answers with a ZIP data package containing::

    ncbi_dataset/data/genomic.fna      nucleotide sequences (FASTA)
    ncbi_dataset/data/data_report.jsonl  per-record metadata (JSON Lines)
    ncbi_dataset/data/dataset_catalog.json
    md5sum.txt

Notes
-----
* The endpoint is a single streaming POST: there is no pagination, the whole
  package arrives in one response.
* The service is rate limited to 5 requests / minute per IP. Retries therefore
  back off for at least ``RATE_LIMIT_SLEEP`` seconds.
* Downloads are written to ``<name>.part`` and only promoted to the final name
  after the ZIP validates, so a failed run can never clobber a good dataset.

Usage
-----
    python3 src/ncbi_client.py --taxon 10404
    python3 src/ncbi_client.py --taxon 10404 --out-dir data/raw --force
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import socket
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download"
DEFAULT_TAXON = "10404"  # Hepadnaviridae (HBV + related hepadnaviruses)
CHUNK_SIZE = 1 << 20  # 1 MiB
RATE_LIMIT_SLEEP = 15  # seconds; service allows 5 requests / minute
PACKAGE_ROOT = "ncbi_dataset/data"
REQUIRED_MEMBERS = (f"{PACKAGE_ROOT}/genomic.fna", f"{PACKAGE_ROOT}/data_report.jsonl")


class DownloadError(RuntimeError):
    """Raised when the data package could not be retrieved or verified."""


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def build_payload(
    taxons: list[str],
    include: list[str],
    complete_only: bool = False,
    refseq_only: bool = False,
    annotated_only: bool = False,
    host: str = "",
    geo_location: str = "",
    usa_state: str = "",
    released_since: str = "0001-01-01T00:00:00Z",
    updated_since: str = "0001-01-01T00:00:00Z",
) -> dict:
    """Build the JSON body accepted by the virus genome download endpoint."""
    if not taxons:
        raise ValueError("at least one taxon is required")
    return {
        "annotated_only": annotated_only,
        "complete_only": complete_only,
        "format": "tsv",
        "geo_location": geo_location,
        "host": host,
        "include_sequence": [item.strip().upper() for item in include if item.strip()],
        "pangolin_classification": "",
        "refseq_only": refseq_only,
        "released_since": released_since,
        "taxons": [str(t) for t in taxons],
        "updated_since": updated_since,
        "usa_state": usa_state,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def _retry_delay(exc: Exception, attempt: int) -> float:
    """Delay before retrying, honouring Retry-After when NCBI supplies it."""
    if isinstance(exc, urllib.error.HTTPError):
        retry_after = exc.headers.get("Retry-After") if exc.headers else None
        if retry_after:
            try:
                return max(float(retry_after), 1.0)
            except ValueError:
                pass
    return max(RATE_LIMIT_SLEEP, RATE_LIMIT_SLEEP * attempt)


def download_package(
    url: str,
    payload: dict,
    destination: Path,
    api_key: str | None = None,
    attempts: int = 4,
    timeout: int = 600,
    show_progress: bool = True,
) -> int:
    """Stream the ZIP package to ``destination`` atomically. Returns byte count."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    body = json.dumps(payload).encode("utf-8")
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/zip",
        "User-Agent": "HBV-Genome-Pipeline/1.0 (data-collection)",
    }
    if api_key:
        headers["api-key"] = api_key

    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        request = urllib.request.Request(url, data=body, headers=headers, method="POST")
        written = 0
        try:
            log(f"[ncbi] request attempt {attempt}/{attempts} -> {url}")
            with urllib.request.urlopen(request, timeout=timeout) as response:
                total = int(response.headers.get("Content-Length") or 0)
                with partial.open("wb") as handle:
                    while True:
                        chunk = response.read(CHUNK_SIZE)
                        if not chunk:
                            break
                        handle.write(chunk)
                        written += len(chunk)
                        if show_progress and written % (20 * CHUNK_SIZE) < CHUNK_SIZE:
                            if total:
                                pct = 100.0 * written / total
                                log(f"[ncbi] {written / 1e6:8.1f} MB / {total / 1e6:.1f} MB ({pct:.0f}%)")
                            else:
                                log(f"[ncbi] {written / 1e6:8.1f} MB")
            if written == 0:
                raise DownloadError("server returned an empty response")
            os.replace(partial, destination)
            return written
        except urllib.error.HTTPError as exc:
            partial.unlink(missing_ok=True)
            last_error = exc
            retryable = exc.code in (429, 500, 502, 503, 504)
            log(f"[ncbi] HTTP {exc.code} on attempt {attempt}: {exc.reason}")
            if not retryable or attempt == attempts:
                break
        except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError, DownloadError) as exc:
            partial.unlink(missing_ok=True)
            last_error = exc
            log(f"[ncbi] network error on attempt {attempt}: {exc}")
            if attempt == attempts:
                break
        delay = _retry_delay(last_error, attempt) if isinstance(last_error, Exception) else RATE_LIMIT_SLEEP
        log(f"[ncbi] retrying in {delay:.0f}s (rate limit is 5 requests/minute)")
        time.sleep(delay)

    partial.unlink(missing_ok=True)
    raise DownloadError(f"download failed after {attempts} attempts: {last_error}")


def validate_zip(path: Path) -> None:
    """Verify the archive is readable, complete, and contains the expected files."""
    if not zipfile.is_zipfile(path):
        raise DownloadError(f"{path} is not a valid ZIP archive")
    with zipfile.ZipFile(path) as archive:
        bad = archive.testzip()
        if bad is not None:
            raise DownloadError(f"corrupt member in archive: {bad}")
        names = set(archive.namelist())
        missing = [member for member in REQUIRED_MEMBERS if member not in names]
        if missing:
            raise DownloadError(f"archive is missing required members: {', '.join(missing)}")


def extract_zip(path: Path, destination: Path) -> Path:
    """Replace ``destination`` with a clean extraction of ``path``."""
    staging = destination.with_name(destination.name + ".tmp")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    with zipfile.ZipFile(path) as archive:
        root = staging.resolve()
        for member in archive.infolist():
            target = (staging / member.filename).resolve()
            if not str(target).startswith(str(root)):
                raise DownloadError(f"unsafe path in archive: {member.filename}")
        archive.extractall(staging)

    if destination.exists():
        shutil.rmtree(destination)
    staging.rename(destination)
    return destination


def main(argv: list[str] | None = None) -> int:
    base_dir = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxon", action="append", default=None,
                        help=f"NCBI taxonomy id or name (repeatable, default {DEFAULT_TAXON})")
    parser.add_argument("--include", default="genome",
                        help="comma separated sequence types: genome,cds,protein,annotation,biosample (default: genome)")
    parser.add_argument("--complete-only", action="store_true", help="restrict to complete sequences")
    parser.add_argument("--refseq-only", action="store_true", help="restrict to RefSeq records")
    parser.add_argument("--annotated-only", action="store_true", help="restrict to annotated records")
    parser.add_argument("--host", default="", help="filter by host species")
    parser.add_argument("--geo-location", default="", help="filter by continent or country")
    parser.add_argument("--out-dir", default=str(base_dir / "data" / "raw"), help="output directory")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="NCBI API key (optional)")
    parser.add_argument("--force", action="store_true", help="re-download even if a valid package exists")
    parser.add_argument("--attempts", type=int, default=4, help="download attempts (default 4)")
    parser.add_argument("--timeout", type=int, default=900, help="per-request timeout in seconds")
    parser.add_argument("--quiet", action="store_true", help="suppress progress output")
    args = parser.parse_args(argv)

    taxons = args.taxon or [DEFAULT_TAXON]
    include = [item for item in args.include.split(",") if item.strip()]
    out_dir = Path(args.out_dir).expanduser().resolve()
    archive = out_dir / "ncbi_dataset.zip"
    extract_dir = out_dir / "extracted"

    payload = build_payload(
        taxons=taxons,
        include=include,
        complete_only=args.complete_only,
        refseq_only=args.refseq_only,
        annotated_only=args.annotated_only,
        host=args.host,
        geo_location=args.geo_location,
    )

    try:
        if archive.exists() and not args.force:
            previous = out_dir / "retrieval.json"
            previous_request = None
            if previous.exists():
                try:
                    previous_request = json.loads(previous.read_text(encoding="utf-8")).get("request")
                except (json.JSONDecodeError, OSError):
                    previous_request = None
            if previous_request != payload:
                log("[ncbi] cached package was built with different parameters; re-downloading")
                archive.unlink(missing_ok=True)
            else:
                try:
                    validate_zip(archive)
                    log(f"[ncbi] reusing validated package {archive} (use --force to re-download)")
                except DownloadError as exc:
                    log(f"[ncbi] existing package is invalid ({exc}); re-downloading")
                    archive.unlink(missing_ok=True)

        if not archive.exists():
            download_package(
                API_URL, payload, archive,
                api_key=args.api_key, attempts=args.attempts,
                timeout=args.timeout, show_progress=not args.quiet,
            )
            log(f"[ncbi] downloaded {archive.stat().st_size / 1e6:.1f} MB")

        validate_zip(archive)
        extract_zip(archive, extract_dir)
        log(f"[ncbi] extracted to {extract_dir}")

        manifest = {
            "source": "ncbi_datasets_virus",
            "endpoint": API_URL,
            "retrieved_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "request": payload,
            "archive": str(archive),
            "archive_bytes": archive.stat().st_size,
            "archive_sha256": sha256_file(archive),
            "extracted_dir": str(extract_dir),
        }
        manifest_path = out_dir / "retrieval.json"
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(manifest, indent=2))
        return 0
    except (DownloadError, ValueError, OSError) as exc:
        log(f"[ncbi] ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
