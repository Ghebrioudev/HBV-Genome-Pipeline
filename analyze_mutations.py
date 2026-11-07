#!/usr/bin/env python3
"""Task 3: align HBV sequences to NC_003977.2 and call reference-based variants.

Inputs (Task 2 outputs):
  data/processed/hbv_sequences.fasta
  data/processed/hbv_metadata.csv

Outputs (created only after a successful run):
  results/mutations.csv
  results/mutation_analysis_skipped.csv
  results/mutation_analysis.json

MAFFT is an external system dependency. Sequences are added as fragments to a
single-reference alignment in bounded batches. The low-cost --6merpair strategy
avoids the expensive pairwise-local-alignment path that crashed in the supplied
MAFFT 7.525 run. Failed batches are retried as progressively smaller batches.
Raw sequences are never compared by unaligned character position.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Sequence


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "data" / "processed" / "hbv_sequences.fasta"
DEFAULT_METADATA = ROOT / "data" / "processed" / "hbv_metadata.csv"
DEFAULT_OUTPUT = ROOT / "results" / "mutations.csv"
DEFAULT_SKIPPED = ROOT / "results" / "mutation_analysis_skipped.csv"
DEFAULT_REPORT = ROOT / "results" / "mutation_analysis.json"
DEFAULT_REFERENCE = "NC_003977.2"
DNA_IUPAC = set("ACGTRYSWKMBDHVN")
CANONICAL = set("ACGT")

MUTATION_FIELDS = [
    "sample_id",
    "reference_accession",
    "position",
    "reference_base",
    "observed_base",
    "variation_type",
    "call_status",
    "aligned_bases",
    "percent_identity",
    "completeness",
]
SKIPPED_FIELDS = ["sample_id", "sequence_length", "completeness", "reason"]


@dataclass(frozen=True)
class FastaRecord:
    """A FASTA record; identifier is the first whitespace-delimited header token."""

    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class AlignmentMetrics:
    """Coverage and identity diagnostics for one aligned sample."""

    aligned_bases: int
    matches: int
    percent_identity: float
    first_reference_position: int
    last_reference_position: int


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
    """Read FASTA records without requiring a third-party bioinformatics package."""
    header: str | None = None
    sequence_parts: list[str] = []

    try:
        handle = path.open("r", encoding="utf-8")
    except OSError as exc:
        raise ValueError(f"Cannot open FASTA input '{path}': {exc}") from exc

    with handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield _make_record(header, sequence_parts, path)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at {path}:{line_number}")
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence data appears before the first FASTA header at "
                        f"{path}:{line_number}"
                    )
                sequence_parts.append("".join(line.split()).upper())

    if header is not None:
        yield _make_record(header, sequence_parts, path)


def _make_record(header: str, sequence_parts: list[str], path: Path) -> FastaRecord:
    identifier = header.split()[0]
    return FastaRecord(identifier, header, "".join(sequence_parts).upper())


def load_metadata(path: Path) -> dict[str, dict[str, str]]:
    """Load Task 2 metadata indexed by accession, preserving original values."""
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "accession" not in reader.fieldnames:
                raise ValueError("Metadata CSV must contain an 'accession' column")
            metadata: dict[str, dict[str, str]] = {}
            for row in reader:
                accession = (row.get("accession") or "").strip()
                if not accession:
                    continue
                if accession in metadata:
                    raise ValueError(f"Duplicate metadata accession: {accession}")
                metadata[accession] = row
            return metadata
    except OSError as exc:
        raise ValueError(f"Cannot open metadata CSV '{path}': {exc}") from exc


def find_reference(input_fasta: Path, accession: str) -> FastaRecord:
    """Find and validate the exact reference record bundled in the Task 2 FASTA."""
    found: FastaRecord | None = None
    for record in iter_fasta(input_fasta):
        if record.identifier == accession:
            if found is not None:
                raise ValueError(f"Reference accession occurs more than once: {accession}")
            found = record
    if found is None:
        raise ValueError(
            f"Reference {accession} is not present in {input_fasta}. "
            "Provide a Task 2 FASTA containing the reference record."
        )
    if not found.sequence:
        raise ValueError(f"Reference sequence {accession} is empty")
    invalid = sorted(set(found.sequence) - CANONICAL)
    if invalid:
        raise ValueError(
            f"Reference {accession} contains non-ACGT symbols: {''.join(invalid)}"
        )
    return found


def assess_alignment(
    reference_aligned: str,
    sample_aligned: str,
) -> AlignmentMetrics | None:
    """Calculate concrete-base identity and observed reference-coordinate span.

    Ambiguous sample bases count toward the observed coverage span but not toward
    the identity denominator. Alignment gaps do not count as observed bases.
    """
    if len(reference_aligned) != len(sample_aligned):
        raise ValueError("Aligned reference and sample have different lengths")

    reference_position = 0
    covered_positions: list[int] = []
    aligned_bases = 0
    matches = 0

    for ref_base, sample_base in zip(reference_aligned.upper(), sample_aligned.upper()):
        if ref_base != "-":
            reference_position += 1
        if ref_base != "-" and sample_base != "-":
            covered_positions.append(reference_position)
            if ref_base in CANONICAL and sample_base in CANONICAL:
                aligned_bases += 1
                if ref_base == sample_base:
                    matches += 1

    if not covered_positions:
        return None

    identity = (100.0 * matches / aligned_bases) if aligned_bases else 0.0
    return AlignmentMetrics(
        aligned_bases=aligned_bases,
        matches=matches,
        percent_identity=identity,
        first_reference_position=min(covered_positions),
        last_reference_position=max(covered_positions),
    )


def _variant_row(
    sample_id: str,
    reference_accession: str,
    position: str,
    reference_base: str,
    observed_base: str,
    variation_type: str,
    call_status: str,
    metrics: AlignmentMetrics,
    completeness: str,
) -> dict[str, str | int | float]:
    return {
        "sample_id": sample_id,
        "reference_accession": reference_accession,
        "position": position,
        "reference_base": reference_base,
        "observed_base": observed_base,
        "variation_type": variation_type,
        "call_status": call_status,
        "aligned_bases": metrics.aligned_bases,
        "percent_identity": f"{metrics.percent_identity:.4f}",
        "completeness": completeness,
    }


def call_variants_from_alignment(
    sample_id: str,
    reference_accession: str,
    reference_aligned: str,
    sample_aligned: str,
    *,
    min_aligned_bases: int = 80,
    min_identity: float = 70.0,
    completeness: str = "",
) -> tuple[list[dict[str, str | int | float]], AlignmentMetrics | None, str | None]:
    """Call substitutions, internal deletions, insertions and uncertain bases.

    Coordinates are 1-based on the ungapped reference. An insertion between
    reference positions N and N+1 is written as ``N^N+1``. A sample's terminal
    alignment gaps are unobserved coverage, not deletions. Only gaps between
    reference positions actually spanned by sample bases can be deletions.

    Returns (variant rows, metrics, skip_reason). Ambiguous IUPAC symbols are
    emitted as uncertain records and never as definitive substitutions.
    """
    reference_aligned = reference_aligned.upper()
    sample_aligned = sample_aligned.upper()
    metrics = assess_alignment(reference_aligned, sample_aligned)
    if metrics is None:
        return [], None, "no_reference_overlap"
    if metrics.aligned_bases < min_aligned_bases:
        return [], metrics, "insufficient_aligned_bases"
    if metrics.percent_identity < min_identity:
        return [], metrics, "identity_below_threshold"

    # Map each alignment column to its 1-based reference coordinate, or None
    # for columns that represent an insertion relative to the reference.
    column_positions: list[int | None] = []
    ref_position = 0
    for base in reference_aligned:
        if base == "-":
            column_positions.append(None)
        else:
            ref_position += 1
            column_positions.append(ref_position)

    rows: list[dict[str, str | int | float]] = []
    coverage_start = metrics.first_reference_position
    coverage_end = metrics.last_reference_position

    # Nucleotide substitutions, ambiguous calls, and internal deletions.
    for column, (ref_base, sample_base) in enumerate(
        zip(reference_aligned, sample_aligned)
    ):
        position = column_positions[column]
        if position is None or not (coverage_start <= position <= coverage_end):
            continue

        if sample_base == "-":
            rows.append(
                _variant_row(
                    sample_id, reference_accession, str(position), ref_base, "-",
                    "deletion", "definitive", metrics, completeness,
                )
            )
        elif sample_base in CANONICAL:
            if ref_base in CANONICAL and sample_base != ref_base:
                rows.append(
                    _variant_row(
                        sample_id, reference_accession, str(position), ref_base,
                        sample_base, "substitution", "definitive", metrics,
                        completeness,
                    )
                )
        elif sample_base in DNA_IUPAC:
            rows.append(
                _variant_row(
                    sample_id, reference_accession, str(position), ref_base,
                    sample_base, "ambiguous", "uncertain", metrics, completeness,
                )
            )

    # Consecutive insertions are represented as one event anchored between the
    # nearest reference bases. End-of-fragment insertions are not called because
    # one of their flanking reference coordinates is not observed.
    column = 0
    while column < len(reference_aligned):
        if reference_aligned[column] != "-" or sample_aligned[column] == "-":
            column += 1
            continue

        start = column
        inserted: list[str] = []
        while (
            column < len(reference_aligned)
            and reference_aligned[column] == "-"
            and sample_aligned[column] != "-"
        ):
            inserted.append(sample_aligned[column])
            column += 1

        left_position = next(
            (p for p in reversed(column_positions[:start]) if p is not None), None
        )
        right_position = next(
            (p for p in column_positions[column:] if p is not None), None
        )
        if (
            left_position is None
            or right_position is None
            or left_position < coverage_start
            or right_position > coverage_end
        ):
            continue

        inserted_sequence = "".join(inserted)
        if all(base in CANONICAL for base in inserted_sequence):
            variation_type, call_status = "insertion", "definitive"
        else:
            variation_type, call_status = "ambiguous", "uncertain"
        rows.append(
            _variant_row(
                sample_id,
                reference_accession,
                f"{left_position}^{right_position}",
                "-",
                inserted_sequence,
                variation_type,
                call_status,
                metrics,
                completeness,
            )
        )

    return rows, metrics, None


def get_mafft_version(executable: str) -> str:
    """Return the MAFFT version or raise an actionable dependency error."""
    resolved = shutil.which(executable)
    if resolved is None:
        raise RuntimeError(
            f"MAFFT executable '{executable}' was not found. Install MAFFT first "
            "(Ubuntu/Debian: sudo apt install mafft; conda: conda install -c "
            "bioconda mafft), then verify with: mafft --version"
        )
    try:
        result = subprocess.run(
            [resolved, "--version"], capture_output=True, text=True, check=True
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(f"Could not execute MAFFT '{resolved}': {exc}") from exc
    version_text = (result.stdout + "\n" + result.stderr).strip()
    version = version_text.splitlines()[0] if version_text else "version not reported"
    match = re.search(r"\bv?(\d+)\.(\d+)", version)
    if match:
        major, minor = int(match.group(1)), int(match.group(2))
        if major < 7 or (major == 7 and minor < 487):
            raise RuntimeError(
                f"MAFFT {version} is too old for this workflow. Use MAFFT 7.487 or newer."
            )
    return version


def run_mafft_batch(
    executable: str,
    reference_path: Path,
    fragments_path: Path,
    alignment_path: Path,
    threads: int,
) -> None:
    """Run MAFFT fragment addition against the fixed reference alignment."""
    resolved = shutil.which(executable) or executable
    command = [
        resolved,
        "--6merpair",
        "--addfragments",
        str(fragments_path),
        "--thread",
        str(threads),
        str(reference_path),
    ]
    try:
        with alignment_path.open("w", encoding="utf-8") as output_handle:
            result = subprocess.run(
                command,
                stdout=output_handle,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
    except OSError as exc:
        raise RuntimeError(f"Failed to start MAFFT: {exc}") from exc
    if result.returncode != 0:
        details = (result.stderr or "").strip()[-4000:]
        raise RuntimeError(
            f"MAFFT exited with status {result.returncode}. "
            f"Command: {' '.join(command)}\n{details}\n"
            "The batch may be too large or may contain a sequence that triggers "
            "an MAFFT failure. Task 3 can retry smaller batches automatically."
        )


def validate_alignment(
    alignment_path: Path,
    reference_accession: str,
    expected_ids: set[str],
    reference_sequence: str,
) -> dict[str, str]:
    """Read a MAFFT result and confirm it retains the reference and all samples."""
    aligned: dict[str, str] = {}
    for record in iter_fasta(alignment_path):
        if record.identifier in aligned:
            raise RuntimeError(f"MAFFT returned duplicate alignment ID: {record.identifier}")
        aligned[record.identifier] = record.sequence

    if reference_accession not in aligned:
        raise RuntimeError("MAFFT output does not contain the reference sequence")
    if aligned[reference_accession].replace("-", "") != reference_sequence:
        raise RuntimeError("MAFFT changed the reference sequence; refusing to call variants")
    missing = expected_ids - set(aligned)
    extra = set(aligned) - expected_ids - {reference_accession}
    if missing or extra:
        raise RuntimeError(
            "MAFFT output IDs do not match the input batch "
            f"(missing={len(missing)}, unexpected={len(extra)})"
        )

    lengths = {len(sequence) for sequence in aligned.values()}
    if len(lengths) != 1:
        raise RuntimeError("MAFFT output contains sequences of unequal aligned length")
    return aligned


def remove_unverified_header_only_outputs(
    output_csv: Path, skipped_csv: Path, report_json: Path
) -> None:
    """Remove header-only CSV templates when no completed run report exists.

    A header-only CSV can be a legitimate zero-row result only when accompanied
    by a successful JSON report. Without that report it is merely an unverified
    template and could be mistaken for completed analysis output.
    """
    if report_json.is_file():
        return
    for csv_path in (output_csv, skipped_csv):
        if not csv_path.is_file():
            continue
        try:
            with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader, None)
                has_data = any(row for row in reader)
            if header and not has_data:
                csv_path.unlink()
        except OSError as exc:
            raise RuntimeError(
                f"Cannot inspect existing Task 3 output '{csv_path}': {exc}"
            ) from exc


def _write_fasta(path: Path, records: Sequence[FastaRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(f">{record.description}\n")
            for start in range(0, len(record.sequence), 80):
                handle.write(record.sequence[start : start + 80] + "\n")


def _write_single_fasta(path: Path, accession: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(f">{accession}\n")
        for start in range(0, len(sequence), 80):
            handle.write(sequence[start : start + 80] + "\n")


def analyze_dataset(
    *,
    input_fasta: Path,
    metadata_csv: Path,
    output_csv: Path,
    skipped_csv: Path,
    report_json: Path,
    reference_accession: str = DEFAULT_REFERENCE,
    mafft: str = "mafft",
    min_sequence_length: int = 100,
    min_aligned_bases: int = 80,
    min_identity: float = 70.0,
    batch_size: int = 500,
    threads: int = 1,
) -> dict[str, object]:
    """Align eligible Task 2 records and atomically write Task 3 outputs."""
    for path in (input_fasta, metadata_csv):
        if not path.is_file():
            raise ValueError(f"Required input file does not exist: {path}")
    if min_sequence_length < 1 or min_aligned_bases < 1 or batch_size < 1 or threads < 1:
        raise ValueError("Length, batch-size and thread thresholds must be positive")
    if not 0.0 <= min_identity <= 100.0:
        raise ValueError("--min-identity must be between 0 and 100")

    metadata = load_metadata(metadata_csv)
    reference = find_reference(input_fasta, reference_accession)
    mafft_version = get_mafft_version(mafft)

    for path in (output_csv, skipped_csv, report_json):
        path.parent.mkdir(parents=True, exist_ok=True)

    remove_unverified_header_only_outputs(output_csv, skipped_csv, report_json)

    output_temp = output_csv.with_name(output_csv.name + ".tmp")
    skipped_temp = skipped_csv.with_name(skipped_csv.name + ".tmp")
    report_temp = report_json.with_name(report_json.name + ".tmp")
    for path in (output_temp, skipped_temp, report_temp):
        if path.exists():
            path.unlink()

    counts: Counter[str] = Counter()
    variants_by_type: Counter[str] = Counter()
    samples_with_variants: set[str] = set()
    seen_ids: set[str] = set()
    batch: list[FastaRecord] = []
    batch_count = 0
    consecutive_single_sequence_alignment_failures = 0
    run_started = datetime.now(timezone.utc).isoformat()

    try:
        with tempfile.TemporaryDirectory(prefix="hbv_task3_") as temporary_dir:
            temp_root = Path(temporary_dir)
            reference_path = temp_root / "reference.fasta"
            _write_single_fasta(reference_path, reference_accession, reference.sequence)

            with output_temp.open("w", encoding="utf-8", newline="") as mutation_handle, \
                    skipped_temp.open("w", encoding="utf-8", newline="") as skipped_handle:
                mutation_writer = csv.DictWriter(mutation_handle, fieldnames=MUTATION_FIELDS)
                skipped_writer = csv.DictWriter(skipped_handle, fieldnames=SKIPPED_FIELDS)
                mutation_writer.writeheader()
                skipped_writer.writeheader()

                def skip(record: FastaRecord, reason: str, completeness: str) -> None:
                    skipped_writer.writerow({
                        "sample_id": record.identifier,
                        "sequence_length": len(record.sequence),
                        "completeness": completeness,
                        "reason": reason,
                    })
                    counts["skipped_sequences"] += 1
                    counts[f"skipped:{reason}"] += 1

                def process_batch(records: list[FastaRecord]) -> None:
                    nonlocal batch_count, consecutive_single_sequence_alignment_failures
                    if not records:
                        return
                    batch_count += 1
                    fragments_path = temp_root / f"fragments_{batch_count:06d}.fasta"
                    alignment_path = temp_root / f"alignment_{batch_count:06d}.fasta"
                    print(
                        f"[Task 3] Aligning batch {batch_count} ({len(records)} sequences)...",
                        file=sys.stderr,
                    )
                    _write_fasta(fragments_path, records)
                    try:
                        run_mafft_batch(
                            mafft, reference_path, fragments_path, alignment_path, threads
                        )
                        expected_ids = {record.identifier for record in records}
                        aligned = validate_alignment(
                            alignment_path, reference_accession, expected_ids,
                            reference.sequence,
                        )
                    except RuntimeError as exc:
                        fragments_path.unlink(missing_ok=True)
                        alignment_path.unlink(missing_ok=True)
                        if len(records) > 1:
                            midpoint = len(records) // 2
                            print(
                                f"[Task 3] WARNING: MAFFT failed for a batch of "
                                f"{len(records)}. Retrying as {midpoint} and "
                                f"{len(records) - midpoint} sequences. Details: {exc}",
                                file=sys.stderr,
                            )
                            process_batch(records[:midpoint])
                            process_batch(records[midpoint:])
                            return

                        # A single problematic record should not prevent the other
                        # 130k+ records from being processed. Repeated singleton
                        # failures, however, indicate a systemic MAFFT problem.
                        consecutive_single_sequence_alignment_failures += 1
                        if consecutive_single_sequence_alignment_failures >= 5:
                            raise RuntimeError(
                                "MAFFT failed on five consecutive single-sequence "
                                "retries; stopping because this appears to be a "
                                "systemic alignment failure. Last error: " + str(exc)
                            ) from exc
                        record = records[0]
                        completeness = metadata[record.identifier].get("completeness", "")
                        skip(record, "mafft_alignment_failed", completeness)
                        print(
                            f"[Task 3] WARNING: skipped {record.identifier}; "
                            f"MAFFT could not align this record: {exc}",
                            file=sys.stderr,
                        )
                        return

                    consecutive_single_sequence_alignment_failures = 0
                    for record in records:
                        completeness = metadata[record.identifier].get("completeness", "")
                        rows, metrics, skip_reason = call_variants_from_alignment(
                            record.identifier,
                            reference_accession,
                            aligned[reference_accession],
                            aligned[record.identifier],
                            min_aligned_bases=min_aligned_bases,
                            min_identity=min_identity,
                            completeness=completeness,
                        )
                        if skip_reason:
                            skip(record, skip_reason, completeness)
                            continue
                        counts["analyzed_sequences"] += 1
                        if rows:
                            samples_with_variants.add(record.identifier)
                        for row in rows:
                            mutation_writer.writerow(row)
                            variants_by_type[str(row["variation_type"])] += 1
                            counts["mutation_rows"] += 1
                            if row["call_status"] == "uncertain":
                                counts["uncertain_rows"] += 1
                            else:
                                counts["definitive_rows"] += 1
                    print(
                        f"[Task 3] Finished batch {batch_count}; analyzed "
                        f"{counts['analyzed_sequences']} sequences so far.",
                        file=sys.stderr,
                    )
                    # Release parsed alignments and temporary files between batches.
                    del aligned
                    fragments_path.unlink(missing_ok=True)
                    alignment_path.unlink(missing_ok=True)

                for record in iter_fasta(input_fasta):
                    counts["input_sequences"] += 1
                    if record.identifier in seen_ids:
                        raise ValueError(f"Duplicate FASTA accession: {record.identifier}")
                    seen_ids.add(record.identifier)
                    if record.identifier == reference_accession:
                        counts["reference_records"] += 1
                        continue

                    completeness = metadata.get(record.identifier, {}).get("completeness", "")
                    if record.identifier not in metadata:
                        skip(record, "missing_metadata", completeness)
                        continue
                    if not record.sequence:
                        skip(record, "empty_sequence", completeness)
                        continue
                    if "-" in record.sequence:
                        skip(record, "input_contains_alignment_gap", completeness)
                        continue
                    invalid_symbols = sorted(set(record.sequence) - DNA_IUPAC)
                    if invalid_symbols:
                        skip(record, "unsupported_symbols:" + "".join(invalid_symbols), completeness)
                        continue
                    if len(record.sequence) < min_sequence_length:
                        skip(record, "below_min_sequence_length", completeness)
                        continue

                    batch.append(record)
                    counts["eligible_sequences"] += 1
                    if len(batch) >= batch_size:
                        process_batch(batch)
                        batch = []
                process_batch(batch)

        if counts["input_sequences"] == 0:
            raise ValueError(f"Input FASTA contains no records: {input_fasta}")
        if counts["reference_records"] != 1:
            raise ValueError(
                f"Expected exactly one reference record {reference_accession}; "
                f"found {counts['reference_records']}"
            )

        # Keep the report compact and machine-readable; skipped reasons are also
        # available as one row per accession in mutation_analysis_skipped.csv.
        report: dict[str, object] = {
            "schema_version": 1,
            "stage": "Task 3 — Sequence Comparison and Mutation Detection",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "started_at": run_started,
            "input_fasta": str(input_fasta.resolve()),
            "metadata_csv": str(metadata_csv.resolve()),
            "reference": {
                "accession": reference_accession,
                "description": reference.description,
                "length_nt": len(reference.sequence),
                "coordinate_system": "1-based, ungapped reference coordinates",
            },
            "alignment": {
                "tool": "MAFFT",
                "version": mafft_version,
                "method": "--6merpair --addfragments",
                "failure_recovery": (
                    "Failed batches are recursively split in half; individual "
                    "records that still fail are logged as mafft_alignment_failed. "
                    "Five consecutive singleton failures stop the run."
                ),
                "batch_size": batch_size,
                "threads": threads,
                "reference_anchor_preserved": True,
            },
            "calling_policy": {
                "minimum_input_sequence_length_nt": min_sequence_length,
                "minimum_aligned_canonical_bases": min_aligned_bases,
                "minimum_pairwise_identity_percent": min_identity,
                "partial_sequences": "Terminal unobserved coverage is not called as deletion.",
                "ambiguous_bases": "IUPAC ambiguity symbols are emitted as uncertain, not definitive, calls.",
                "insertions": "Between-base coordinate N^N+1; terminal insertions are not called.",
                "deletions": "Only alignment gaps within observed reference coverage are called.",
                "orientation": "Input orientation is retained; no automatic reverse-complement search.",
                "circular_origin": "Linear reference coordinates; origin-spanning alignment is a limitation.",
            },
            "counts": {
                "input_sequences": counts["input_sequences"],
                "reference_records": counts["reference_records"],
                "eligible_sequences": counts["eligible_sequences"],
                "analyzed_sequences": counts["analyzed_sequences"],
                "skipped_sequences": counts["skipped_sequences"],
                "sequences_with_variants": len(samples_with_variants),
                "alignment_batches": batch_count,
                "mutation_rows": counts["mutation_rows"],
                "definitive_rows": counts["definitive_rows"],
                "uncertain_rows": counts["uncertain_rows"],
                "variants_by_type": dict(sorted(variants_by_type.items())),
                "skipped_by_reason": {
                    key.split(":", 1)[1]: value
                    for key, value in sorted(counts.items())
                    if key.startswith("skipped:")
                },
            },
            "outputs": {
                "mutations_csv": str(output_csv),
                "skipped_csv": str(skipped_csv),
            },
        }
        report_temp.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        os.replace(output_temp, output_csv)
        os.replace(skipped_temp, skipped_csv)
        os.replace(report_temp, report_json)
        return report
    except Exception:
        for path in (output_temp, skipped_temp, report_temp):
            path.unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Align Task 2 HBV sequences to the NC_003977.2 reference with MAFFT "
            "and export reference-coordinate mutation events."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Task 2 FASTA input")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA, help="Task 2 metadata CSV")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Mutation CSV output")
    parser.add_argument("--skipped", type=Path, default=DEFAULT_SKIPPED, help="Skipped-sequences CSV output")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT, help="JSON run report output")
    parser.add_argument("--reference-accession", default=DEFAULT_REFERENCE, help="Reference FASTA accession")
    parser.add_argument("--mafft", default="mafft", help="MAFFT executable (default: mafft)")
    parser.add_argument("--min-sequence-length", type=int, default=100,
                        help="Skip shorter input sequences (default: 100 nt)")
    parser.add_argument("--min-aligned-bases", type=int, default=80,
                        help="Minimum concrete aligned base pairs required (default: 80)")
    parser.add_argument("--min-identity", type=float, default=70.0,
                        help="Minimum pairwise identity in percent (default: 70)")
    parser.add_argument("--batch-size", type=int, default=500,
                        help="Maximum fragments per MAFFT invocation (default: 500)")
    parser.add_argument("--threads", type=int, default=1,
                        help="MAFFT threads per batch (default: 1)")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        report = analyze_dataset(
            input_fasta=args.input,
            metadata_csv=args.metadata,
            output_csv=args.output,
            skipped_csv=args.skipped,
            report_json=args.report,
            reference_accession=args.reference_accession,
            mafft=args.mafft,
            min_sequence_length=args.min_sequence_length,
            min_aligned_bases=args.min_aligned_bases,
            min_identity=args.min_identity,
            batch_size=args.batch_size,
            threads=args.threads,
        )
    except (ValueError, RuntimeError) as exc:
        print(f"[Task 3] ERROR: {exc}", file=sys.stderr)
        return 2

    counts = report["counts"]
    assert isinstance(counts, dict)
    print("[Task 3] Sequence comparison completed successfully")
    print(f"[Task 3] Reference: {args.reference_accession}")
    print(f"[Task 3] Input records: {counts['input_sequences']}")
    print(f"[Task 3] Analyzed sequences: {counts['analyzed_sequences']}")
    print(f"[Task 3] Skipped sequences: {counts['skipped_sequences']}")
    print(f"[Task 3] Mutation rows: {counts['mutation_rows']}")
    print(f"[Task 3] Variants by type: {counts['variants_by_type']}")
    print(f"[Task 3] Mutation table: {args.output}")
    print(f"[Task 3] Run report: {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
