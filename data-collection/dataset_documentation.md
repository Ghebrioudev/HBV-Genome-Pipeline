# Dataset documentation — HBV nucleotide dataset (Phase 1)

This file records the provenance, scope, size and integrity of the dataset used by
the pipeline. Values below come from an actual retrieval, not estimates. The
authoritative machine-readable copy of the last retrieval is:

- `data/raw/retrieval.json` — request payload, timestamp, archive checksum
- `data/checksums/ncbi_validation.json` — record counts, sizes, per-file SHA-256
- `data/processed/manifest.json` — checksums of the normalised outputs

## 1. Dataset name and description

**HBV genome dataset (Hepadnaviridae, NCBI Virus / NCBI Datasets).**

Nucleotide sequences and associated metadata for Hepatitis B Virus (HBV) and the
wider *Hepadnaviridae* family, retrieved programmatically from the NCBI Virus
service through the official NCBI Datasets REST API. The dataset is the primary
input for every downstream stage of the project (preprocessing, mutation
detection, statistics).

## 2. Original NCBI Virus query

The dataset was requested to match this NCBI Virus web page:

```
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Hepadnaviridae%20(hepatitis%20B-type%20viruses),%20taxid:10404
```

- `SeqType_s=Nucleotide` — nucleotide records only
- `VirusLineage_ss=Hepadnaviridae ... taxid:10404` — family-level filter

## 3. API used

| Item | Value |
|---|---|
| Service | NCBI Datasets — Virus (v2) |
| Endpoint | `POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download` |
| Equivalent CLI | `datasets download virus genome taxon 10404 --include genome` |
| REST docs | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/ |
| Virus how-to | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/virus-download/ |
| Virus docs | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/help/ |
| Data package spec | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/ |

The `GET /datasets/v2/virus/genome/taxon/{id}/download` form does **not** exist
(it returns HTTP 404); the download is a single streaming `POST` with a JSON
body. There is no pagination — the whole package is returned in one response.

### Exact request body

```json
{
  "annotated_only": false,
  "complete_only": false,
  "format": "tsv",
  "geo_location": "",
  "host": "",
  "include_sequence": ["GENOME"],
  "pangolin_classification": "",
  "refseq_only": false,
  "released_since": "0001-01-01T00:00:00Z",
  "taxons": ["10404"],
  "updated_since": "0001-01-01T00:00:00Z",
  "usa_state": ""
}
```

### Rate limits

- `X-Ratelimit-Limit: 5` — **5 requests per minute** per IP.
- The client backs off at least 15 s between retries and honours `Retry-After`.
- An API key (`--api-key` / `NCBI_API_KEY`) is optional for this endpoint.

## 4. Taxonomic scope requested vs. records actually returned

**Requested:** *Hepadnaviridae* (taxid **10404**) — the family, exactly as the web
query specifies.

**Returned:** 142,254 records, of which **140,866 (99.0 %)** have *Hepatitis B
virus* in their lineage and **1,388 (1.0 %)** are other hepadnaviruses.

| Group | Records | Share |
|---|---:|---:|
| Hepatitis B virus (incl. genotypes) | 140,866 | 99.02 % |
| Woodchuck hepatitis virus | 490 | 0.34 % |
| Duck hepatitis B virus | 275 | 0.19 % |
| Domestic cat hepadnavirus | 86 | 0.06 % |
| Equine / bat / other hepadnaviruses | 537 | 0.38 % |

> **Important:** querying the family is *not* equivalent to querying HBV alone.
> Downstream stages that require HBV only must filter on `lineage` /
> `organism_name` in `metadata.csv`. Using taxid **10407** (`--taxon 10407`)
> retrieves *Hepatitis B virus* only.

## 5. Retrieval details

| Field | Value |
|---|---|
| Retrieval date (UTC) | 2026-10-28T19:15Z |
| Tool | `data-collection/download.sh` (stdlib Python 3.12 client) |
| Sequence count | **142,254** |
| Total nucleotides | **131,955,960** |
| Complete records | 14,100 |
| Partial records | 128,154 |
| Longest record | 10,123 nt |
| Archive size | 31,512,362 bytes (30.1 MiB) |
| Extracted package size | **405.55 MiB** |
| Normalised output size | 190.6 MiB |

## 6. File formats

| File | Format | Description |
|---|---|---|
| `data/raw/ncbi_dataset.zip` | ZIP | Raw NCBI data package (as downloaded) |
| `data/raw/extracted/ncbi_dataset/data/genomic.fna` | FASTA | Nucleotide sequences |
| `data/raw/extracted/ncbi_dataset/data/data_report.jsonl` | JSON Lines | Full per-record metadata (kept raw) |
| `data/raw/extracted/ncbi_dataset/data/dataset_catalog.json` | JSON | Package file inventory |
| `data/raw/extracted/md5sum.txt` | text | NCBI-provided MD5 checksums |
| `data/processed/sequences.fasta` | FASTA | Normalised sequences (downstream input) |
| `data/processed/metadata.csv` | CSV | Normalised metadata (derived view) |
| `data/processed/manifest.json` | JSON | Output checksums and counts |

The raw `data_report.jsonl` is preserved untouched; `metadata.csv` is a derived,
portable view and never replaces it.

## 7. Metadata fields

`metadata.csv` columns:

| Column | Meaning |
|---|---|
| `accession` | GenBank nucleotide accession (versioned) |
| `length` | Sequence length in nucleotides |
| `completeness` | `COMPLETE` or `PARTIAL` |
| `is_annotated` | Whether the record carries annotation |
| `gene_count`, `protein_count` | Number of annotated genes / proteins |
| `source_database` | `GenBank` or `RefSeq` |
| `release_date`, `update_date` | NCBI release / update timestamps |
| `organism_name` | Organism as assigned by NCBI |
| `tax_id` | NCBI taxonomy id of the record |
| `isolate` | Isolate name, when available |
| `lineage` | Full taxonomic lineage (`;`-separated) |
| `submitter` | Submitter name(s), when available |

Fields that NCBI does not provide for a record are left empty. The source of truth
for all fields is `data_report.jsonl`.

## 8. Accession numbers

Records are identified by versioned GenBank nucleotide accessions (e.g.
`NC_003977.2`, the HBV reference genome). The complete accession list is the
`accession` column of `metadata.csv`; 142,254 unique accessions are present.

## 9. Integrity verification

SHA-256 is calculated by the pipeline for the archive, the extracted sequences and
the normalised outputs.

| Artefact | SHA-256 |
|---|---|
| `ncbi_dataset.zip` | `299a16e549f9bb6a5159c9f88d87b017ed21211d69ad37e9bf1c36da06e111e4` |
| `genomic.fna` | `b39fe94375f8e8f3f9dfcf02e275a8f2d4c8472110a1b0e41a6a5086603eeb3a` |
| `sequences.fasta` | `b39fe94375f8e8f3f9dfcf02e275a8f2d4c8472110a1b0e41a6a5086603eeb3a` |
| `metadata.csv` | `9abb45ae4caa3fa5001ba775fe2613220a551410adc3a5791c9f586820874bff` |

Full checksum lists: `data/checksums/dataset.sha256`,
`data/checksums/ncbi_extracted.sha256`.

Verification steps applied: HTTP success, ZIP integrity (`ZipFile.testzip`),
required members present, FASTA parse, IUPAC character check, record count, and
the 100 MB – 1 GB size gate.

## 10. Dataset size requirement

The project requires a dataset between **100 MB and 1 GB**.

- Measured extracted size: **405.55 MiB** — inside the window.
- The archive itself is only **30.1 MiB** (highly compressible), so the size gate
  is evaluated on **extracted bytes**, not on the compressed download.
- `complete_only=true` cannot satisfy the 100 MB floor (complete genomes alone are
  ~44 MB); partial sequences are deliberately included.

## 11. Known limitations

- **Family-wide scope:** ~1 % of records are non-HBV hepadnaviruses (see §4).
- **Mostly partial sequences:** 90 % of records are flagged `PARTIAL`; sequence
  lengths range from 1 nt to 10,123 nt.
- **Heterogeneous metadata:** fields such as `isolate` and `submitter` are absent
  for many records; there is no guaranteed country/collection-date column in the
  package report.
- **Metadata is large:** `data_report.jsonl` is ~258 MB of the 405 MB package.
- **No pagination:** a single POST returns everything; interrupted downloads are
  retried from scratch (the client retries up to 4 times).
- **Serotype/genotype** are not normalised into a dedicated column; they must be
  derived from `organism_name` / `lineage`.
- The dataset reflects NCBI content **as of the retrieval date** and is not
  version-pinned upstream.

## 12. How to reproduce

```bash
# Default: NCBI API, taxon 10404
./download.sh

# HBV only
./download.sh --taxon 10407

# Force a fresh download (ignore cached package)
./download.sh --force

# Local dataset instead of the API
./download.sh --file path/to/hbv.fasta
```

Outputs:

```
data-collection/data/raw/ncbi_dataset.zip
data-collection/data/raw/extracted/ncbi_dataset/data/{genomic.fna,data_report.jsonl}
data-collection/data/processed/{sequences.fasta,metadata.csv,manifest.json}
data-collection/data/checksums/{dataset.sha256,ncbi_validation.json}
```

A full run of the default query took ~64 s and ~400 MB of disk.
