# Baylor Genetics Technical Assessment

A two-part pipeline that extracts gene mentions from a scientific paper, enriches them with metadata from public databases, and models the results in a normalized SQL database.

Paper: PMID 38790019 - Diagnostic yield of exome and genome sequencing after non-diagnostic multi-gene panels in patients with single-system diseases


## File Reference

 `drpt_genes.py`:Part 1 — full pipeline 
 `gene_output.csv`: Output from Part 1 
 `gene_output.json`: Same data in JSON format 
 `gene_schema.sql`:Part 2 — schema
 `README.md`: This file 

---

## Quick Start

```bash
pip install requests beautifulsoup4 lxml
python3 drpt_genes.py
sqlite3 gene.db < gene_schema.sql
```

---

## Part 1 — Gene Extraction Pipeline (`drpt_genes.py`)

```
PMC full text → regex extraction → HGNC validation → Ensembl coords → disease resolution → CSV + JSON
```

1. Fetch — Downloads the paper XML from NCBI E-utilities. Bibliography nodes are stripped before text extraction to avoid false positives from cited genes.
2. Extract — A regex matches HGNC-format tokens (`[A-Z][A-Z0-9]{1,5}`). A stop-token list removes common non-gene acronyms (DNA, PCR, FDA, etc.). Remaining candidates are ranked by frequency of mention so the most prominent genes are processed first.
3. Validate — Each candidate is looked up via the HGNC REST API. Any symbol not found in HGNC is discarded. This is the primary false-positive filter.
4. Coordinates — Genomic coordinates are fetched from Ensembl for both GRCh38 (hg38) and GRCh37 (hg19) using separate endpoints, formatted as `chr:start-end`.
5. Disease — The Ensembl gene ID from the HGNC record is used to query Open Targets (aggregates OMIM, Orphanet, ClinVar). Returned disease names are cross-referenced against the paper text, only diseases actually mentioned in the paper are kept. If none are confirmed, the pipeline falls back to extracting disease-relevant sentences directly from the paper.

### Output fields

- HGNC ID | HGNC REST API |
- HGNC Gene Name | HGNC REST API |
- Gene Aliases | HGNC (`alias_symbol` + `prev_symbol`) |
- hg38 Coordinates | Ensembl GRCh38 |
- hg19 Coordinates | Ensembl GRCh37 |
- Disease | Open Targets, confirmed against paper text |

### Output files

- `gene_output.csv`
- `gene_output.json`

---

## Part 2 — SQL Schema (`gene_schema.sql`)

### What it does

Models the CSV output as a normalized relational database and provides two queries against it.

### Schema (4 tables)

```
gene          — one row per confirmed gene
gene_alias    — one row per alias, FK → gene       (one-to-many)
disease       — one row per unique disease name
gene_disease  — gene ↔ disease bridge table        (many-to-many)
```

The flat CSV stores aliases and diseases as pipe-separated strings (`SRN1 | PDCN`). The schema normalizes these into separate tables so each value is individually queryable, joinable, and free of string-parsing logic. Diseases get their own table because multiple genes share the same disease (e.g., both NPHS2 and APOL1 link to focal segmental glomerulosclerosis), storing the name once and referencing it via a bridge table avoids duplication and update anomalies.

