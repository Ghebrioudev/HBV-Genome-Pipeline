# Sources and provenance

## 1. Dataset

| Field | Value |
|---|---|
| Name | NCBI Virus — *Hepadnaviridae* nucleotide genome data package (NCBI Datasets "virus genome" package) |
| Identifier | NCBI Taxonomy ID **10404** (*Hepadnaviridae*); HBV records (taxon **10407**, *Hepatitis B virus*) are selected by `preprocess.py` |
| Provider | National Center for Biotechnology Information (NCBI), NCBI Virus / NCBI Datasets v2 |
| Download address | `POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download` with body `{"taxons": ["10404"], "include_sequence": ["GENOME"]}` |
| Equivalent web page | <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=taxid:10404> |
| Access | public, free, no account; optional API key (`NCBI_API_KEY`) only raises the rate limit; no personal or controlled-access data |
| Download method | automatic: `./run.sh` → `data-collection/download.sh` → `data-collection/src/ncbi_client.py` (Python `urllib`, retries, ZIP validation) |
| Download date of the reference snapshot | **2026-11-14 12:36:39 UTC** (first retrieval); re-verified 2026-11-15 00:18 UTC with the identical archive SHA-256 |
| Format | ZIP containing `genomic.fna` (FASTA), `data_report.jsonl` (metadata, JSON Lines), `dataset_catalog.json`, `md5sum.txt` |
| Size — downloaded ZIP | 32,898,096 bytes (31.4 MiB, compressed) |
| Size — dataset (extracted package) | **425,368,669 bytes = 405.7 MiB** (`genomic.fna` 146,731,669 B + `data_report.jsonl` 278,634,873 B + small files) — inside the 100 MB – 1 GB range required, checked automatically by `dataset_validator.py --size-policy hard` |
| Content | 142,288 nucleotide records, 132,005,421 nt (140,889 HBV records after filtering) |

### Integrity

| Check | Where | Value |
|---|---|---|
| NCBI MD5 of every file in the package | `data-collection/data/raw/extracted/md5sum.txt`, verified by `md5sum -c` in `run.sh` | `genomic.fna` = `63b752be24d2064e39088f5960bad3f2` |
| SHA-256 of the downloaded ZIP | `data-collection/data/checksums/dataset.sha256`, `data-collection/data/raw/retrieval.json` | `968ce5165734b61df7164451d57de429c0ac22317997b3d18e88ff9d42bd324b` |
| SHA-256 of the extracted files | `data-collection/data/checksums/ncbi_extracted.sha256` | `genomic.fna` = `b02f7a2ebf188467b92ef2be03e40556d458bed05b655c3db8aeb0c1c9f3b3d1` |
| SHA-256 of the normalised Task 1 files | `data-collection/data/checksums/dataset.sha256` | `sequences.fasta` = `b02f7a2e…`, `metadata.csv` = `e90499ba688d6ca37f00617ec68151425b1b71242a9f25686da9194b6aa5e202` |
| SHA-256 of the cleaned HBV files | `data/processed/manifest.json` (written and checked by `preprocess.py`) | `hbv_sequences.fasta` = `a2b5ccad…` |
| SHA-256 of all results | `results/checksums.sha256`, written then verified by `run.sh` | — |

All lists are recomputed at each download and **verified at each run**
(`sha256sum -c`, `md5sum -c`); a mismatch stops the pipeline with exit code 5.
After verification the raw package is made read-only.

### Reproducibility note

NCBI Virus is updated continuously. A new download (`./run.sh -f`) therefore
returns more records than the reference snapshot, with a new date and new
checksums, all recorded automatically in `data-collection/data/raw/retrieval.json`
(`retrieved_at` = actual download time, `verified_at` = last check) and in the
checksum files. The numbers in `README.md` and `rapport.pdf` describe the
reference snapshot above.

## 2. Reference genome

| Field | Value |
|---|---|
| Accession | `NC_003977.2` (NCBI RefSeq) |
| Organism | *Hepatitis B virus* (taxon 10407), strain ayw |
| Length | 3,182 nt |
| Address | <https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2> |
| How it is obtained | it is one of the records of the downloaded package; `analyze_mutations.py` extracts it by accession |

## 3. Software

| Tool | Version | Use | Address |
|---|---|---|---|
| MAFFT | ≥ 7.487 (reference run: v7.525; clean-machine test: v7.505) | `--6merpair --addfragments` alignment | <https://mafft.cbrc.jp/alignment/software/> |
| Python | ≥ 3.11 | all analysis scripts | <https://www.python.org> |
| pandas | 3.1.1 (`requirements.txt`) | statistics | <https://pandas.pydata.org> |
| matplotlib | 3.11.2 (`requirements.txt`) | figures | <https://matplotlib.org> |
| GNU coreutils, grep, sed, awk, tar | system | shell stages, checksums, archive | — |

The exact MAFFT version and parameters used are stored in
`results/mutation_analysis.json`.

## 4. Documentation consulted

- NCBI Virus help: <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/>
- NCBI Datasets REST API: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/>
- NCBI Datasets virus how-to: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/>
- Virus genome data package: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/>
- MAFFT `--addfragments`: <https://mafft.cbrc.jp/alignment/software/addsequences.html>
