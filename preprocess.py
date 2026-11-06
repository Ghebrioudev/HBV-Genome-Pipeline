#!/usr/bin/env python3
"""Preprocess the Task 1 HBV FASTA/metadata dataset.

Inputs (read-only):
    data-collection/data/processed/sequences.fasta
    data-collection/data/processed/metadata.csv
    data-collection/data/processed/manifest.json

Outputs:
    data/processed/hbv_sequences.fasta
    data/processed/hbv_metadata.csv
    data/processed/manifest.json
    results/preprocessing/*.csv and cleaning_report.json

The implementation uses only the Python standard library. It does not download,
modify, align, or silently discard biological sequence data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

SCRIPT_VERSION = "1.0.0"
HBV_TAXON_ID = "10407"
HEPADNAVIRIDAE_NAME = "hepadnaviridae"
HBV_NAME = "hepatitis b virus"
DNA_IUPAC = frozenset("ACGTRYSWKMBDHVN")
AMBIGUOUS_IUPAC = frozenset("RYSWKMBDHVN")
REVIEW_SYMBOLS = frozenset("U-")
UNSUPPORTED_SYMBOL_RE = re.compile(r"[^ACGTRYSWKMBDHVNU-]")
AMBIGUOUS_SYMBOL_RE = re.compile(r"[RYSWKMBDHVN]")
FASTA_WRAP = 80

REPORT_COLUMNS = {
    "rejected_records.csv": [
        "record_number", "record_type", "accession", "header", "reason",
        "sequence_length", "sequence_sha256", "symbols",
    ],
    "excluded_non_hbv.csv": [
        "record_number", "accession", "organism_name", "tax_id", "lineage", "reason",
    ],
    "unresolved_taxonomy.csv": [
        "record_number", "accession", "organism_name", "tax_id", "lineage", "reason",
    ],
    "review_records.csv": [
        "record_number", "accession", "reason", "sequence_length",
        "metadata_length", "contains_u", "contains_gap", "sequence_sha256",
    ],
    "duplicate_accessions.csv": [
        "accession", "record_number", "sequence_sha256", "action", "reason",
    ],
    "identical_sequences.csv": [
        "sequence_group", "sequence_sha256", "sequence_length", "group_size", "accession",
    ],
}

METADATA_EXTRA_COLUMNS = [
    "sequence_length",
    "sequence_sha256",
    "ambiguous_count",
    "ambiguous_percentage",
    "contains_ambiguity",
    "duplicate_sequence_group",
    "metadata_length_status",
    "validation_status",
]


class PreprocessingError(RuntimeError):
    """Raised when the Task 1 input contract or preprocessing fails."""


@dataclass(frozen=True)
class FastaRecord:
    """One FASTA record, with its original header and normalized sequence."""

    record_number: int
    header: str
    accession: str
    sequence: str


@dataclass
class ProcessingResult:
    """Collected outputs and statistics for one preprocessing run."""

    retained_records: List[Tuple[FastaRecord, Dict[str, str], Dict[str, object]]]
    rejected: List[Dict[str, object]]
    excluded_non_hbv: List[Dict[str, object]]
    unresolved_taxonomy: List[Dict[str, object]]
    review_records: List[Dict[str, object]]
    duplicate_accessions: List[Dict[str, object]]
    identical_sequences: List[Dict[str, object]]
    counts: Dict[str, object]
    length_stats: Dict[str, object]
    ambiguity_stats: Dict[str, object]
    completeness_counts: Dict[str, int]
    input_metadata_only: int
    metadata_columns: List[str]


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of a file, streaming it in bounded memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(sequence: str) -> str:
    """Hash a normalized nucleotide sequence."""
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def iter_fasta(path: Path) -> Iterable[FastaRecord]:
    """Yield FASTA records, removing formatting whitespace and uppercasing DNA.

    The header is preserved (apart from surrounding whitespace). A sequence
    before the first header is treated as a malformed FASTA input.
    """
    header: Optional[str] = None
    sequence_parts: List[str] = []
    record_number = 0

    def build_record(current_header: str, parts: List[str], number: int) -> FastaRecord:
        token = current_header.split()[0] if current_header.split() else ""
        sequence = "".join(parts).upper()
        return FastaRecord(number, current_header, token, sequence)

    with path.open("r", encoding="utf-8", errors="replace", newline=None) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    record_number += 1
                    yield build_record(header, sequence_parts, record_number)
                header = line[1:].strip()
                sequence_parts = []
            elif not line.strip():
                continue
            elif header is None:
                raise PreprocessingError(
                    "Malformed FASTA: sequence data occurs before the first header "
                    "(line {}).".format(line_number)
                )
            else:
                # FASTA whitespace is formatting, not sequence content.
                sequence_parts.append("".join(line.split()))

    if header is not None:
        record_number += 1
        yield build_record(header, sequence_parts, record_number)


def load_metadata(path: Path) -> Tuple[Dict[str, List[Dict[str, str]]], List[str], int]:
    """Load CSV metadata indexed by exact, versioned accession."""
    by_accession: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "accession" not in reader.fieldnames:
            raise PreprocessingError("metadata.csv must contain an 'accession' column.")
        fieldnames = list(reader.fieldnames)
        row_count = 0
        for row in reader:
            row_count += 1
            accession = (row.get("accession") or "").strip()
            row["accession"] = accession
            by_accession[accession].append(row)
    return by_accession, fieldnames, row_count


def verify_input_manifest(
    manifest_path: Path, sequences_path: Path, metadata_path: Path
) -> Tuple[Dict[str, object], Dict[str, Dict[str, object]]]:
    """Verify Task 1 file hashes against its manifest before processing."""
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PreprocessingError("Cannot read Task 1 manifest: {}".format(exc)) from exc

    entries = {}
    for item in manifest.get("files", []):
        if isinstance(item, dict) and item.get("name"):
            entries[str(item["name"])] = item

    actual_files = {
        "sequences.fasta": sequences_path,
        "metadata.csv": metadata_path,
    }
    verified: Dict[str, Dict[str, object]] = {}
    for name, path in actual_files.items():
        if not path.is_file():
            raise PreprocessingError("Required Task 1 input is missing: {}".format(path))
        if name not in entries:
            raise PreprocessingError(
                "Task 1 manifest does not contain a checksum entry for {}.".format(name)
            )
        digest = sha256_file(path)
        expected = str(entries[name].get("sha256", "")).lower()
        expected_bytes = entries[name].get("bytes")
        actual_bytes = path.stat().st_size
        if expected_bytes is not None and int(expected_bytes) != actual_bytes:
            raise PreprocessingError(
                "Byte-size mismatch for {}: expected {}, calculated {}.".format(
                    path, expected_bytes, actual_bytes
                )
            )
        if not expected or digest != expected:
            raise PreprocessingError(
                "Checksum mismatch for {}: expected {}, calculated {}.".format(
                    path, expected or "<missing>", digest
                )
            )
        verified[name] = {
            "path": str(path),
            "bytes": actual_bytes,
            "sha256": digest,
        }
    verified["manifest.json"] = {
        "path": str(manifest_path),
        "bytes": manifest_path.stat().st_size,
        "sha256": sha256_file(manifest_path),
    }
    return manifest, verified


def taxonomy_classification(metadata: Dict[str, str]) -> Tuple[str, str]:
    """Classify HBV using the NCBI taxid and the NCBI lineage path.

    A record is HBV when its own taxid is 10407 or its lineage contains the exact
    HBV taxon name (which also recognizes descendants). A known non-HBV member
    of Hepadnaviridae is excluded. Incomplete or contradictory taxonomy is
    reported as unresolved instead of being guessed.
    """
    tax_id = (metadata.get("tax_id") or "").strip()
    lineage = (metadata.get("lineage") or "").strip()
    parts = [part.strip().casefold() for part in lineage.split(";") if part.strip()]
    has_hbv_ancestor = HBV_NAME in parts
    has_family = HEPADNAVIRIDAE_NAME in parts

    if tax_id == HBV_TAXON_ID:
        if lineage and not has_hbv_ancestor:
            return "unresolved", "tax_id_is_HBV_but_lineage_does_not_confirm_HBV"
        return "hbv", "tax_id_10407"
    if has_hbv_ancestor:
        return "hbv", "lineage_contains_HBV_taxon_or_descendant"
    if has_family:
        return "non_hbv", "lineage_is_Hepadnaviridae_without_HBV_ancestor"
    return "unresolved", "insufficient_NCBI_taxonomy_information"


def _record_report_row(record: FastaRecord, reason: str, digest: str = "", symbols: str = "") -> Dict[str, object]:
    return {
        "record_number": record.record_number,
        "record_type": "fasta",
        "accession": record.accession,
        "header": record.header,
        "reason": reason,
        "sequence_length": len(record.sequence),
        "sequence_sha256": digest,
        "symbols": symbols,
    }


def _metadata_length_status(metadata: Dict[str, str], actual_length: int) -> Tuple[str, Optional[str]]:
    raw_length = (metadata.get("length") or "").strip()
    if not raw_length:
        return "missing", None
    try:
        expected_length = int(raw_length)
    except ValueError:
        return "invalid", "metadata_length_is_not_an_integer"
    if expected_length != actual_length:
        return "mismatch", "metadata_length_does_not_match_FASTA_sequence"
    return "match", None


def _write_csv(path: Path, columns: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _length_statistics(lengths: Sequence[int]) -> Dict[str, object]:
    if not lengths:
        return {"count": 0, "minimum": None, "maximum": None, "mean": None, "median": None}
    ordered = sorted(lengths)
    count = len(ordered)
    middle = count // 2
    median = ordered[middle] if count % 2 else (ordered[middle - 1] + ordered[middle]) / 2
    return {
        "count": count,
        "minimum": ordered[0],
        "maximum": ordered[-1],
        "mean": round(sum(ordered) / count, 4),
        "median": median,
    }


def process_dataset(
    sequences_path: Path,
    metadata_path: Path,
    manifest_path: Path,
) -> Tuple[ProcessingResult, Dict[str, object], Dict[str, Dict[str, object]]]:
    """Run taxonomy filtering, validation, deduplication and reporting in memory."""
    input_manifest, input_hashes = verify_input_manifest(manifest_path, sequences_path, metadata_path)
    metadata_by_accession, metadata_columns, metadata_row_count = load_metadata(metadata_path)
    print("[preprocess] Loaded {} metadata rows.".format(metadata_row_count), file=sys.stderr, flush=True)
    if not metadata_columns or "accession" not in metadata_columns:
        raise PreprocessingError("Invalid metadata schema: missing accession column.")

    records = list(iter_fasta(sequences_path))
    if not records:
        raise PreprocessingError("The Task 1 FASTA file contains no sequence records.")
    print("[preprocess] Parsed {} FASTA records.".format(len(records)), file=sys.stderr, flush=True)
    manifest_record_count = input_manifest.get("records")
    if manifest_record_count is not None and int(manifest_record_count) != len(records):
        raise PreprocessingError(
            "Task 1 manifest record count ({}) does not match FASTA record count ({}).".format(
                manifest_record_count, len(records)
            )
        )

    rejected: List[Dict[str, object]] = []
    excluded: List[Dict[str, object]] = []
    unresolved: List[Dict[str, object]] = []
    review: List[Dict[str, object]] = []
    duplicate_accessions: List[Dict[str, object]] = []
    retained: List[Tuple[FastaRecord, Dict[str, str], Dict[str, object]]] = []

    # Detect accession duplication before any filtering. Versioned accessions
    # are compared literally, so .1 and .2 remain distinct records.
    fasta_groups: Dict[str, List[FastaRecord]] = defaultdict(list)
    for record in records:
        if not record.accession:
            rejected.append(_record_report_row(record, "missing_accession_in_FASTA_header"))
        else:
            fasta_groups[record.accession].append(record)

    canonical_records: List[FastaRecord] = []
    duplicate_accession_groups = 0
    duplicate_accession_identical_removed = 0
    duplicate_accession_conflict_groups = 0
    duplicate_accession_conflict_records = 0
    for accession, group in fasta_groups.items():
        if len(group) == 1:
            canonical_records.append(group[0])
            continue
        duplicate_accession_groups += 1
        sequences = {item.sequence for item in group}
        if len(sequences) == 1:
            canonical_records.append(group[0])
            for duplicate in group[1:]:
                duplicate_accession_identical_removed += 1
                duplicate_accessions.append({
                    "accession": accession,
                    "record_number": duplicate.record_number,
                    "sequence_sha256": sha256_text(duplicate.sequence),
                    "action": "deduplicated_identical_copy",
                    "reason": "same_versioned_accession_and_identical_normalized_sequence",
                })
        else:
            duplicate_accession_conflict_groups += 1
            duplicate_accession_conflict_records += len(group)
            for conflicting in group:
                digest = sha256_text(conflicting.sequence)
                duplicate_accessions.append({
                    "accession": accession,
                    "record_number": conflicting.record_number,
                    "sequence_sha256": digest,
                    "action": "quarantined_conflicting_accession",
                    "reason": "same_versioned_accession_has_different_sequences",
                })
                rejected.append(_record_report_row(
                    conflicting,
                    "duplicate_accession_conflicting_sequences",
                    digest,
                ))

    # Report metadata rows that cannot be associated with a FASTA record.
    fasta_accessions = set(fasta_groups)
    metadata_only_rows = 0
    for accession, rows in metadata_by_accession.items():
        if accession not in fasta_accessions:
            for row in rows:
                metadata_only_rows += 1
                rejected.append({
                    "record_number": "",
                    "record_type": "metadata",
                    "accession": accession,
                    "header": "",
                    "reason": "metadata_row_without_matching_FASTA_record",
                    "sequence_length": "",
                    "sequence_sha256": "",
                    "symbols": "",
                })

    taxonomy_counts = Counter()
    invalid_symbol_records = 0
    empty_sequence_records = 0
    uracil_records = 0
    gap_records = 0
    ambiguous_positions = 0
    nucleotide_counts: Counter = Counter()
    retained_lengths: List[int] = []
    ambiguity_symbol_counts: Counter = Counter()
    completeness_counts: Counter = Counter()
    missing_metadata_length_count = 0

    # The original metadata CSV is accession-sorted while FASTA order can be
    # different; always join by the exact accession.version, never row index.
    for processed_index, record in enumerate(canonical_records, start=1):
        if processed_index % 20000 == 0:
            print("[preprocess] Validated {}/{} candidate records.".format(processed_index, len(canonical_records)), file=sys.stderr, flush=True)
        metadata_rows = metadata_by_accession.get(record.accession, [])
        if not metadata_rows:
            unresolved.append({
                "record_number": record.record_number,
                "accession": record.accession,
                "organism_name": "",
                "tax_id": "",
                "lineage": "",
                "reason": "no_matching_metadata_row",
            })
            taxonomy_counts["unresolved"] += 1
            continue

        if any(row != metadata_rows[0] for row in metadata_rows[1:]):
            row = metadata_rows[0]
            unresolved.append({
                "record_number": record.record_number,
                "accession": record.accession,
                "organism_name": row.get("organism_name", ""),
                "tax_id": row.get("tax_id", ""),
                "lineage": row.get("lineage", ""),
                "reason": "multiple_conflicting_metadata_rows_for_accession",
            })
            taxonomy_counts["unresolved"] += 1
            continue
        metadata = metadata_rows[0]

        tax_status, tax_reason = taxonomy_classification(metadata)
        taxonomy_counts[tax_status] += 1
        if tax_status == "non_hbv":
            excluded.append({
                "record_number": record.record_number,
                "accession": record.accession,
                "organism_name": metadata.get("organism_name", ""),
                "tax_id": metadata.get("tax_id", ""),
                "lineage": metadata.get("lineage", ""),
                "reason": tax_reason,
            })
            continue
        if tax_status == "unresolved":
            unresolved.append({
                "record_number": record.record_number,
                "accession": record.accession,
                "organism_name": metadata.get("organism_name", ""),
                "tax_id": metadata.get("tax_id", ""),
                "lineage": metadata.get("lineage", ""),
                "reason": tax_reason,
            })
            continue

        sequence = record.sequence
        digest = sha256_text(sequence)
        if not sequence:
            empty_sequence_records += 1
            rejected.append(_record_report_row(record, "empty_sequence", digest))
            continue

        unsupported = sorted(set(UNSUPPORTED_SYMBOL_RE.findall(sequence)))
        contains_u = "U" in sequence
        contains_gap = "-" in sequence
        if contains_u:
            uracil_records += 1
        if contains_gap:
            gap_records += 1

        if unsupported:
            invalid_symbol_records += 1
            reason = "unsupported_nucleotide_symbols:{}".format("".join(unsupported))
            if contains_u or contains_gap:
                reason += ";also_contains_review_symbols:{}".format(
                    "".join(sorted(({"U"} if contains_u else set()) | ({"-"} if contains_gap else set())))
                )
            rejected.append(_record_report_row(record, reason, digest, "".join(unsupported)))
            continue

        length_status, length_problem = _metadata_length_status(metadata, len(sequence))
        if length_status == "missing":
            missing_metadata_length_count += 1
        if contains_u or contains_gap or length_problem:
            reasons = []
            if contains_u:
                reasons.append("contains_U_uracil_review_required")
            if contains_gap:
                reasons.append("contains_alignment_gap_review_required")
            if length_problem:
                reasons.append(length_problem)
            review.append({
                "record_number": record.record_number,
                "accession": record.accession,
                "reason": ";".join(reasons),
                "sequence_length": len(sequence),
                "metadata_length": metadata.get("length", ""),
                "contains_u": str(contains_u).lower(),
                "contains_gap": str(contains_gap).lower(),
                "sequence_sha256": digest,
            })
            continue

        ambiguous_count = AMBIGUOUS_SYMBOL_RE.subn("", sequence)[1]
        nucleotide_counts.update(sequence)
        ambiguous_positions += ambiguous_count
        retained_lengths.append(len(sequence))
        completeness = (metadata.get("completeness") or "UNKNOWN").strip().upper() or "UNKNOWN"
        completeness_counts[completeness] += 1
        derived = {
            "sequence_length": len(sequence),
            "sequence_sha256": digest,
            "ambiguous_count": ambiguous_count,
            "ambiguous_percentage": round(100.0 * ambiguous_count / len(sequence), 6),
            "contains_ambiguity": str(ambiguous_count > 0).lower(),
            "duplicate_sequence_group": "",
            "metadata_length_status": length_status,
            "validation_status": "accepted",
        }
        retained.append((record, metadata, derived))

    # Exact sequence duplicates are grouped independently from accession
    # duplicates. Distinct accessions are all retained.
    hash_groups: Dict[str, List[int]] = defaultdict(list)
    for index, (record, _metadata, derived) in enumerate(retained):
        hash_groups[str(derived["sequence_sha256"])].append(index)
    identical_groups = []
    identical_report_rows: List[Dict[str, object]] = []
    for digest, indices in sorted(hash_groups.items()):
        if len(indices) < 2:
            continue
        # Verify actual sequence equality rather than trusting the hash alone.
        by_sequence: Dict[str, List[int]] = defaultdict(list)
        for index in indices:
            by_sequence[retained[index][0].sequence].append(index)
        for sequence, exact_indices in sorted(by_sequence.items(), key=lambda item: item[0]):
            if len(exact_indices) < 2:
                continue
            group_id = "sha256:" + digest
            accessions = sorted(retained[index][0].accession for index in exact_indices)
            identical_groups.append({
                "sequence_group": group_id,
                "sequence_sha256": digest,
                "sequence_length": len(sequence),
                "group_size": len(exact_indices),
                "accessions": accessions,
            })
            for index in exact_indices:
                record, _metadata, derived = retained[index]
                derived["duplicate_sequence_group"] = group_id
                identical_report_rows.append({
                    "sequence_group": group_id,
                    "sequence_sha256": digest,
                    "sequence_length": len(sequence),
                    "group_size": len(exact_indices),
                    "accession": record.accession,
                })

    # Preserve FASTA ordering in the cleaned FASTA and matching metadata CSV.
    lengths_total = sum(retained_lengths)
    ambiguity_symbol_counts = {symbol: nucleotide_counts.get(symbol, 0) for symbol in sorted(AMBIGUOUS_IUPAC)}
    print("[preprocess] Retained {} validated HBV records.".format(len(retained)), file=sys.stderr, flush=True)
    ambiguity_stats = {
        "ambiguous_positions": ambiguous_positions,
        "retained_nucleotide_positions": lengths_total,
        "ambiguous_percentage": round(100.0 * ambiguous_positions / lengths_total, 6) if lengths_total else 0.0,
        "records_with_ambiguity": sum(1 for _record, _metadata, derived in retained if derived["ambiguous_count"] > 0),
        "counts_by_symbol": {symbol: ambiguity_symbol_counts.get(symbol, 0) for symbol in sorted(AMBIGUOUS_IUPAC)},
    }

    # Metadata-only rows are not FASTA sequence records; they are reported
    # separately in rejected_records.csv and counted independently.
    sequence_rejected_count = sum(1 for row in rejected if row.get("record_type") == "fasta")
    counts = {
        "input_fasta_records": len(records),
        "input_metadata_rows": metadata_row_count,
        "unique_fasta_accessions": len(fasta_groups),
        "retained_records": len(retained),
        "rejected_fasta_records": sequence_rejected_count,
        "rejected_report_rows": len(rejected),
        "excluded_non_hbv_records": len(excluded),
        "unresolved_taxonomy_records": len(unresolved),
        "review_required_records": len(review),
        "duplicate_accession_groups": duplicate_accession_groups,
        "duplicate_accession_identical_records_removed": duplicate_accession_identical_removed,
        "duplicate_accession_conflict_groups": duplicate_accession_conflict_groups,
        "duplicate_accession_conflict_records": duplicate_accession_conflict_records,
        "identical_sequence_groups": len(identical_groups),
        "records_in_identical_sequence_groups": sum(int(group["group_size"]) for group in identical_groups),
        "extra_records_in_identical_sequence_groups": sum(int(group["group_size"]) - 1 for group in identical_groups),
        "metadata_rows_without_FASTA_record": metadata_only_rows,
        "empty_sequence_records": empty_sequence_records,
        "unsupported_symbol_records": invalid_symbol_records,
        "records_containing_U": uracil_records,
        "records_containing_gap": gap_records,
        "records_with_missing_metadata_length": missing_metadata_length_count,
        "taxonomy_classification": {
            "hbv_candidates": taxonomy_counts.get("hbv", 0),
            "non_hbv_excluded": taxonomy_counts.get("non_hbv", 0),
            "unresolved": taxonomy_counts.get("unresolved", 0),
        },
    }

    result = ProcessingResult(
        retained_records=retained,
        rejected=rejected,
        excluded_non_hbv=excluded,
        unresolved_taxonomy=unresolved,
        review_records=review,
        duplicate_accessions=duplicate_accessions,
        identical_sequences=identical_report_rows,
        counts=counts,
        length_stats=_length_statistics(retained_lengths),
        ambiguity_stats=ambiguity_stats,
        completeness_counts=dict(sorted(completeness_counts.items())),
        input_metadata_only=metadata_only_rows,
        metadata_columns=metadata_columns,
    )
    return result, input_manifest, input_hashes


def write_outputs(
    result: ProcessingResult,
    input_manifest: Dict[str, object],
    input_hashes: Dict[str, Dict[str, object]],
    repo_root: Path,
    output_dir: Path,
    reports_dir: Path,
    metadata_columns: Sequence[str],
) -> Dict[str, object]:
    """Write all outputs to staging and promote each completed file atomically."""
    output_dir = output_dir.resolve()
    reports_dir = reports_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".hbv-preprocess-", dir=str(repo_root)) as temp_name:
        stage_root = Path(temp_name)
        stage_data = stage_root / "data" / "processed"
        stage_reports = stage_root / "results" / "preprocessing"
        stage_data.mkdir(parents=True)
        stage_reports.mkdir(parents=True)

        fasta_out = stage_data / "hbv_sequences.fasta"
        metadata_out = stage_data / "hbv_metadata.csv"
        with fasta_out.open("w", encoding="utf-8", newline="\n") as fasta_handle:
            for record, _metadata, _derived in result.retained_records:
                fasta_handle.write(">{}\n".format(record.header))
                sequence = record.sequence
                for offset in range(0, len(sequence), FASTA_WRAP):
                    fasta_handle.write(sequence[offset:offset + FASTA_WRAP] + "\n")

        metadata_output_columns = list(metadata_columns) + [
            column for column in METADATA_EXTRA_COLUMNS if column not in metadata_columns
        ]
        with metadata_out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=metadata_output_columns, extrasaction="ignore")
            writer.writeheader()
            for _record, metadata, derived in result.retained_records:
                row = dict(metadata)
                row.update(derived)
                writer.writerow(row)

        reports = {
            "rejected_records.csv": result.rejected,
            "excluded_non_hbv.csv": result.excluded_non_hbv,
            "unresolved_taxonomy.csv": result.unresolved_taxonomy,
            "review_records.csv": result.review_records,
            "duplicate_accessions.csv": result.duplicate_accessions,
            "identical_sequences.csv": result.identical_sequences,
        }
        for name, rows in reports.items():
            _write_csv(stage_reports / name, REPORT_COLUMNS[name], rows)

        relative_input = {}
        for name, info in input_hashes.items():
            input_path = Path(str(info["path"]))
            try:
                display_path = input_path.resolve().relative_to(repo_root.resolve()).as_posix()
            except ValueError:
                display_path = input_path.name
            relative_input[name] = {
                "path": display_path,
                "bytes": info["bytes"],
                "sha256": info["sha256"],
            }

        report = {
            "schema_version": 1,
            "project": "HBV Genome Pipeline",
            "stage": "Task 2 - Data Preprocessing",
            "input": relative_input,
            "source_manifest_created_at": input_manifest.get("created_at"),
            "counts": result.counts,
            "retained_sequence_length_statistics": result.length_stats,
            "retained_ambiguity_statistics": result.ambiguity_stats,
            "retained_completeness_counts": result.completeness_counts,
            "policies": {
                "taxonomic_scope": "tax_id 10407 or an NCBI lineage containing the exact Hepatitis B virus taxon; known non-HBV Hepadnaviridae are excluded; incomplete/contradictory taxonomy is reviewed",
                "partial_sequences": "Retain valid COMPLETE, PARTIAL and UNKNOWN records; no minimum-length threshold is applied",
                "nucleotide_alphabet": "ACGT plus IUPAC ambiguity symbols RYSWKMBDHVN",
                "ambiguous_symbols": "Preserved and counted; no consensus base is imputed",
                "uracil": "U-containing records are quarantined for manual review; no U-to-T conversion is performed",
                "gaps": "Gap-containing records are quarantined for review; no alignment gaps are removed",
                "unsupported_symbols": "Reject records containing symbols outside the DNA IUPAC alphabet, U and gap",
                "duplicate_accessions": "Compare exact versioned accession strings; identical copies retain the first record; conflicting sequences for one accession.version are quarantined; different versions remain distinct",
                "identical_sequences": "SHA-256 candidate grouping followed by exact sequence equality; different accessions remain retained; reverse complements are not merged",
                "metadata_join": "Join by exact accession.version, not by CSV/FASTA row position",
                "sequence_alignment": "Not performed; alignment belongs to Task 3",
                "raw_inputs": "Task 1 files are read-only and never modified",
            },
            "warnings": [
                "Records with missing metadata length are retained with the original field left blank; actual sequence length is recorded separately.",
                "A valid IUPAC ambiguity symbol is retained but may be unsuitable for a definitive mutation call at that position.",
            ],
        }
        report_path = stage_reports / "cleaning_report.json"
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        staged_outputs = [fasta_out, metadata_out] + [stage_reports / name for name in reports] + [report_path]
        output_entries = []
        for source in staged_outputs:
            if source.parent == stage_data:
                destination = output_dir / source.name
            else:
                destination = reports_dir / source.name
            try:
                display_path = destination.resolve().relative_to(repo_root.resolve()).as_posix()
            except ValueError:
                display_path = destination.name
            output_entries.append({
                "path": display_path,
                "bytes": source.stat().st_size,
                "sha256": sha256_file(source),
            })

        manifest = {
            "schema_version": 1,
            "source": "HBV Genome Pipeline Task 1 processed dataset",
            "source_manifest_created_at": input_manifest.get("created_at"),
            "deterministic_for_identical_inputs_and_configuration": True,
            "input_files": relative_input,
            "outputs": output_entries,
            "record_counts": result.counts,
            "processing": {
                "software": "HBV Genome Pipeline preprocess.py",
                "version": SCRIPT_VERSION,
                "taxon_id": HBV_TAXON_ID,
                "sequence_alignment_performed": False,
                "configuration": {
                    "min_sequence_length": None,
                    "accepted_dna_symbols": "ACGTRYSWKMBDHVN",
                    "preserve_partial_sequences": True,
                    "preserve_iupac_ambiguity": True,
                    "deduplicate_identical_accessions": True,
                    "preserve_identical_sequences_with_different_accessions": True,
                },
            },
        }
        staged_manifest = stage_data / "manifest.json"
        staged_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        # Only after every staged file has been successfully created and hashed
        # do we promote the completed files to their final locations.
        for source in staged_outputs:
            destination = output_dir / source.name if source.parent == stage_data else reports_dir / source.name
            os.replace(str(source), str(destination))
        os.replace(str(staged_manifest), str(output_dir / "manifest.json"))

    return manifest


def build_parser(repo_root: Path) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate and preprocess the Task 1 HBV sequence dataset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=repo_root / "data-collection" / "data" / "processed",
        help="Task 1 processed directory containing sequences.fasta, metadata.csv and manifest.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo_root / "data" / "processed",
        help="Directory for cleaned FASTA, metadata and manifest",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=repo_root / "results" / "preprocessing",
        help="Directory for preprocessing audit reports",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    repo_root = Path(__file__).resolve().parent
    args = build_parser(repo_root).parse_args(argv)
    input_dir = args.input_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    reports_dir = args.reports_dir.expanduser().resolve()
    sequences_path = input_dir / "sequences.fasta"
    metadata_path = input_dir / "metadata.csv"
    manifest_path = input_dir / "manifest.json"

    try:
        result, input_manifest, input_hashes = process_dataset(
            sequences_path, metadata_path, manifest_path
        )
        metadata_columns = result.metadata_columns
        manifest = write_outputs(
            result,
            input_manifest,
            input_hashes,
            repo_root,
            output_dir,
            reports_dir,
            metadata_columns,
        )
    except (PreprocessingError, OSError, csv.Error, ValueError, UnicodeError) as exc:
        print("[preprocess] ERROR: {}".format(exc), file=sys.stderr)
        return 1

    summary = {
        "retained_records": result.counts["retained_records"],
        "rejected_fasta_records": result.counts["rejected_fasta_records"],
        "excluded_non_hbv_records": result.counts["excluded_non_hbv_records"],
        "unresolved_taxonomy_records": result.counts["unresolved_taxonomy_records"],
        "review_required_records": result.counts["review_required_records"],
        "output_dir": str(output_dir),
        "reports_dir": str(reports_dir),
        "manifest_outputs": len(manifest["outputs"]),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
