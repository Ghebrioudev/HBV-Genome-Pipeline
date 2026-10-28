# Data collection (Phase 1)

Acquires the HBV nucleotide dataset for the pipeline, from **NCBI by default**,
with a **local-file** pathway as an alternative. Both pathways produce the same
normalised outputs, so downstream stages never need to know where the data came
from.

## Quick start

```bash
# From the repository root
./download.sh                 # NCBI, Hepadnaviridae (taxon 10404)
./download.sh --taxon 10407   # HBV only
./download.sh --file my.fasta # local dataset
```

## Dependencies

| Requirement | Version | Notes |
|---|---|---|
| Python | ≥ 3.8 (tested 3.12.3) | **standard library only** — no pip installs |
| bash | ≥ 4 | `set -euo pipefail`, arrays |
| `zip`/`unzip` | any | only needed for the tests, not the pipeline |

No `pandas`, `requests`, `samtools` or `datasets` CLI is required. The official
`datasets` CLI is *optional* and is useful only for inspecting the data package.

## What it does

```
                 ┌── NCBI Datasets REST API ──┐
                 │   (default, taxon 10404)   │
input source ────┤                            ├──> normalise ──> verify ──> checksums
                 │   local FASTA/.gz/.zip     │
                 └────────────────────────────┘
```

1. **Retrieve** — `src/ncbi_client.py` POSTs to
   `https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download`, streams the
   ZIP to `<name>.part`, validates it, then promotes it atomically.
2. **Validate** — `src/dataset_validator.py` checks ZIP integrity, FASTA parsing,
   IUPAC characters, record counts and the 100 MB – 1 GB size window.
3. **Normalise** — `src/file_loader.py` writes `sequences.fasta`, `metadata.csv`
   and `manifest.json`. The processed directory is staged and only swapped in
   after every step succeeds.
4. **Checksum** — SHA-256 for the source and every processed artefact.

## Options

| Flag | Default | Meaning |
|---|---|---|
| `--taxon ID` | `10404` | NCBI taxonomy id or name |
| `--include LIST` | `genome` | `genome,cds,protein,annotation,biosample` |
| `--file PATH` | — | Use a local dataset instead of the API |
| `--out-dir DIR` | `data/raw` | Raw download directory |
| `--min-mb N` | `100` | Minimum dataset size |
| `--max-mb N` | `1024` | Maximum dataset size |
| `--api-key KEY` | `$NCBI_API_KEY` | Optional NCBI API key |
| `--force` | off | Re-download even if a valid package is cached |

The local pathway always uses the **warn** size policy (a small test file should
not be fatal); the NCBI pathway uses the **hard** policy.

## Outputs

```
data/
  raw/
    ncbi_dataset.zip            raw package, as downloaded
    retrieval.json              request body + timestamp + archive SHA-256
    extracted/                  unpacked package (genomic.fna, data_report.jsonl, ...)
  processed/
    sequences.fasta             normalised sequences  (downstream input)
    metadata.csv                normalised metadata
    manifest.json               counts + SHA-256 of the outputs
  checksums/
    dataset.sha256              SHA-256 of source + processed files
    ncbi_validation.json        validation report
    ncbi_extracted.sha256       per-file checksums of the extracted package
```

`data/raw/` and `data/processed/` are git-ignored; `data/checksums/` is kept so
integrity values are versioned.

## Standalone validation

```bash
./scripts/validate_dataset.sh --dataset data/raw/extracted --strict
./scripts/validate_dataset.sh --dataset ~/upload/hbv.fasta          # warn only
```

## Tests / manual verification

```bash
bash -n download.sh scripts/validate_dataset.sh
python3 -m py_compile src/*.py

# small dataset smoke test (fast, ~12 KB)
./download.sh --taxon 10406 --min-mb 0

# local pathways
./download.sh --file data-collection/tests/data/sample.fasta
./download.sh --file data-collection/tests/data/sample.fasta.gz --min-mb 0
./download.sh --file data-collection/tests/data/sample.zip       --min-mb 0
```

## Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `HTTP 429` | Rate limit (5 req/min). The client backs off automatically; wait and retry. |
| `HTTP 404` | Wrong endpoint — only `POST /v2/virus/genome/download` exists. |
| `below the required minimum of 100 MB` | Query returned a small subset (e.g. `--complete-only`). Raise scope or lower `--min-mb`. |
| `not a valid ZIP archive` | Interrupted download. Re-run with `--force`. |
| `no FASTA records found` | The local file is not FASTA, or a ZIP has no FASTA members. |
| Stale dataset after changing `--taxon` | Handled automatically: the cached request is compared and re-downloaded when it differs. |

## Dataset documentation

See [`dataset_documentation.md`](dataset_documentation.md) for provenance, scope,
size, metadata fields and checksums.
