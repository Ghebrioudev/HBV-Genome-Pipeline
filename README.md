# HBV-Genome-Pipeline

A reproducible bioinformatics pipeline for analysing **Hepatitis B Virus (HBV)**
genomic data using Shell and Python.

## Status

| Phase | Scope | State |
|---|---|---|
| 1 | Data collection (NCBI API + local file) | ✅ implemented |
| 2 | Preprocessing | not started |
| 3 | Sequence comparison / mutation detection | not started |
| 4 | Statistics and figures | not started |
| 5 | Report | not started |

## Requirements

- Python ≥ 3.8 (tested on 3.12.3) — **standard library only**
- bash ≥ 4

No third-party Python packages are required for Phase 1.

## Quick start

```bash
# Download the HBV/Hepadnaviridae dataset from NCBI (default source)
./download.sh

# Hepatitis B virus only
./download.sh --taxon 10407

# Use a local FASTA / .gz / .zip instead
./download.sh --file path/to/hbv.fasta
```

Outputs land in `data-collection/data/`, and the full provenance is documented in
[`data-collection/dataset_documentation.md`](data-collection/dataset_documentation.md).

## Repository layout

```
download.sh                  entry point (delegates to data-collection/)
data-collection/            Phase 1 — data acquisition
  download.sh               orchestrator (NCBI default, --file alternative)
  dataset_documentation.md  provenance, scope, checksums
  README.md                 detailed setup and troubleshooting
  scripts/validate_dataset.sh
  src/ncbi_client.py        NCBI Datasets REST client
  src/file_loader.py        local + NCBI normalisation
  src/dataset_validator.py  integrity / format / size checks
  data/{raw,processed,checksums}/
preprocess.py                Phase 2 (empty)
analyze_mutations.py         Phase 3 (empty)
statistics.py                Phase 4 (empty)
report/report.md             Phase 5 (empty)
```

See [`AGENTS.md`](AGENTS.md) for conventions used by contributors and AI agents.

## Dataset

Default query: NCBI Virus, *Hepadnaviridae* (taxid 10404), nucleotide records.
Last retrieval: **142,254 sequences · 131,955,960 nt · 405.55 MiB** extracted
(31.5 MB archive). 99.0 % are *Hepatitis B virus*; ~1 % are other hepadnaviruses.
