# Sources and provenance

This file records the external data source and software used by the bundled
HBV Genome Pipeline run.

## NCBI nucleotide data

| Field | Value |
|---|---|
| Provider | NCBI Virus / NCBI Datasets Virus API (v2) |
| Endpoint | `POST https://api.ncbi.nlm.nih.gov/datasets/v2/virus/genome/download` |
| Taxon requested | 10404 — *Hepadnaviridae* family |
| Retrieval timestamp | 2026-11-14 12:36:39 UTC |
| Input records | 142,288 |
| Input nucleotide positions | 132,005,421 |
| Raw ZIP size | 32,898,096 bytes |
| Raw ZIP SHA-256 | `968ce5165734b61df7164451d57de429c0ac22317997b3d18e88ff9d42bd324b` |
| Task 1 processed FASTA SHA-256 | `b02f7a2ebf188467b92ef2be03e40556d458bed05b655c3db8aeb0c1c9f3b3d1` |
| Task 1 processed metadata SHA-256 | `e90499ba688d6ca37f00617ec68151425b1b71242a9f25686da9194b6aa5e202` |

The NCBI family-level request is intentionally broader than HBV alone. Task 2
filters the NCBI metadata by HBV taxon/lineage, excludes 1,399 non-HBV records,
and retains 140,889 HBV records. The Task 1 and Task 2 manifests and checksum
files provide the detailed file-level provenance.

Official references:

- NCBI Virus help: <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/>
- NCBI Datasets REST API: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/>
- NCBI Datasets Virus guide: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/>
- Virus genome data package: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/virus-genome/>
- NCBI Datasets API keys: <https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/>

## Reference genome

| Field | Value |
|---|---|
| Accession | `NC_003977.2` |
| Organism | *Hepatitis B virus* (taxon 10407), strain ayw |
| Length | 3,182 nt |
| Used by | Task 3 mutation detection |
| Record link | <https://www.ncbi.nlm.nih.gov/nuccore/NC_003977.2> |

## Alignment software

| Field | Value |
|---|---|
| Tool | MAFFT |
| Version recorded in run | v7.525 (2024/Mar/13) |
| Supported minimum | 7.487 |
| Method | `--6merpair --addfragments` |
| Default batch size | 500 sample sequences |
| Default thread count | 1 |
| Software documentation | <https://mafft.cbrc.jp/alignment/software/addsequences.html> |
| Installation guidance | <https://mafft.cbrc.jp/alignment/software/linux.html> |

The exact MAFFT version and parameters are stored in
`results/mutation_analysis.json`.

## Python dependencies

Task 1–3 use the Python standard library. Task 4 uses `pandas` and `matplotlib`,
with the versions pinned in `requirements.txt`.
