# Sources and provenance

Records every external data source used by the project: where the data came
from, when it was retrieved, and how it is identified.

## Data sources

### 1. NCBI Virus — *Hepadnaviridae* nucleotide records (primary)

| Field | Value |
|---|---|
| Database | NCBI Virus (nucleotide) / NCBI Datasets |
| Web query | https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Hepadnaviridae%20(hepatitis%20B-type%20viruses),%20taxid:10404 |
| API | `POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download` |
| Taxon | 10404 — *Hepadnaviridae* (family) |
| Retrieved (UTC) | 2026-10-28T19:15Z |
| Records | 142,254 |
| Nucleotides | 131,955,960 |
| Archive | `data-collection/data/raw/ncbi_dataset.zip` (31,512,362 bytes) |
| SHA-256 | `299a16e549f9bb6a5159c9f88d87b017ed21211d69ad37e9bf1c36da06e111e4` |
| Accessions | GenBank nucleotide accessions, listed in `data-collection/data/processed/metadata.csv` |
| Documentation | [`data-collection/dataset_documentation.md`](data-collection/dataset_documentation.md) |

### 2. HBV reference genome

| Field | Value |
|---|---|
| Accession | `NC_003977.2` |
| Organism | Hepatitis B virus (taxid 10407) |
| Included in | the NCBI dataset above (RefSeq record) |
| Used by | Phase 3 (mutation detection), as the coordinate reference |

## Official documentation consulted

- NCBI Virus help — https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/
- NCBI Datasets REST API — https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
- NCBI Datasets virus how-to — https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/
- Virus data package reference — https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/
- NCBI Datasets API keys — https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/

## Notes

- The web query filters on the *family* (taxid 10404). The API returns the same
  family scope, so ~1 % of records are non-HBV hepadnaviruses. See
  `dataset_documentation.md` §4 for the measured breakdown.
- NCBI content changes over time; the counts above are a snapshot and are not
  version-pinned upstream.
- Rate limit on the virus download endpoint: 5 requests per minute per IP.
