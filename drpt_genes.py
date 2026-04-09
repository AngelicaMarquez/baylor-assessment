"""
Baylor Genetics Technical Assessment - Part 1
Fetches PMID 38790019 (PMC11127317), extracts HGNC gene symbols,
and enriches each with metadata from HGNC, Ensembl (hg38 + hg19).

Disease strategy — combined approach:
  1. Query MyGene.info for curated disease names associated with the gene
     (aggregated from OMIM, HPO, DisGeNET — no API key required).
  2. Cross-reference each disease name against the paper text to confirm
     it is actually mentioned in this specific paper.
  3. Fall back to sentence-level keyword extraction if no curated name
     is confirmed in the paper.
"""

import re
import json
import time
import logging
import csv
from dataclasses import dataclass, field, asdict
from typing import Optional

import requests
from bs4 import BeautifulSoup

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PMC_FETCH_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=pmc&id=11127317&rettype=full&retmode=xml"
)
HGNC_SEARCH_URL  = "https://rest.genenames.org/search/symbol/{symbol}"
HGNC_FETCH_URL   = "https://rest.genenames.org/fetch/symbol/{symbol}"
ENSEMBL38_URL    = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
ENSEMBL37_URL    = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json"
OPENTARGETS_URL   = "https://api.platform.opentargets.org/api/v4/graphql"

# HGNC gene symbols: 2–6 uppercase letters optionally followed by digits/hyphens.
# We require word boundaries to avoid matching abbreviations inside words.
GENE_SYMBOL_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,5}(?:-[A-Z0-9]+)?)\b")

# Common non-gene uppercase tokens to skip (acronyms, units, section headers…)
STOP_TOKENS: set[str] = {
    "DNA", "RNA", "NGS", "PCR", "SNP", "CNV", "VUS", "LOF", "NMD",
    "USA", "UK", "EU", "FDA", "ACMG", "IRB", "CI", "OR", "HR",
    "ES", "GS", "WES", "WGS", "IGV", "OMIM", "UCSC", "NCBI",
    "CDS", "UTR", "INDEL", "SV", "MT", "ATP", "ADP", "BP", "KB",
    "MB", "GB", "NA", "SD", "IQR", "N", "P", "HGNC", "EGBP",
    "PRaUD", "Mayo", "MRI", "CT", "QC", "SOP", "DECIPHER", "CLINVAR",
    "HPO", "GO", "KEGG", "GWAS", "ENCODE", "GTEx", "TCGA",
}

REQUEST_DELAY = 0.4  # seconds between API calls (be polite)


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class GeneRecord:
    hgnc_id: str = ""
    hgnc_gene_name: str = ""
    gene_aliases: str = ""          # pipe-separated list
    hg38_coordinates: str = ""      # chr:start-end
    hg19_coordinates: str = ""      # chr:start-end
    disease: str = ""               # curated name (MyGene/OMIM) confirmed in paper, or sentence fallback


# ---------------------------------------------------------------------------
# Step 1: Fetch full-text XML from PMC
# ---------------------------------------------------------------------------

def fetch_paper_text(url: str = PMC_FETCH_URL) -> str:
    """
    Download the full-text PMC XML and return all plain text extracted from it.
    Raises requests.HTTPError on non-2xx responses.
    """
    log.info("Fetching paper from PMC …")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    soup = BeautifulSoup(resp.text, "lxml-xml")

    # Remove reference/bibliography nodes – they inflate gene-hit counts
    for tag in soup.find_all(["ref-list", "back"]):
        tag.decompose()

    text = soup.get_text(separator=" ")
    log.info("Paper fetched – %d characters of text", len(text))
    return text


# ---------------------------------------------------------------------------
# Step 2: Extract candidate gene symbols
# ---------------------------------------------------------------------------

def extract_gene_candidates(text: str) -> list[str]:
    """
    Return a deduplicated, ordered list of candidate HGNC-format gene symbols
    found in *text*, filtering out common non-gene uppercase tokens.
    """
    seen: dict[str, int] = {}
    for match in GENE_SYMBOL_RE.finditer(text):
        symbol = match.group(1)
        if symbol not in STOP_TOKENS:
            seen[symbol] = seen.get(symbol, 0) + 1

    # Sort by frequency descending so the most-mentioned genes come first
    ranked = sorted(seen.items(), key=lambda kv: kv[1], reverse=True)
    log.info("Found %d candidate tokens; top 20: %s",
             len(ranked), [s for s, _ in ranked[:20]])
    return [s for s, _ in ranked]


# ---------------------------------------------------------------------------
# Step 3: Validate symbols via HGNC API
# ---------------------------------------------------------------------------

def validate_with_hgnc(symbol: str) -> Optional[dict]:
    """
    Query the HGNC REST API for *symbol*.
    Returns the HGNC document dict if found, else None.
    """
    url = HGNC_FETCH_URL.format(symbol=symbol)
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        docs = data.get("response", {}).get("docs", [])
        if docs:
            return docs[0]
    except Exception as exc:
        log.warning("HGNC lookup failed for %s: %s", symbol, exc)
    finally:
        time.sleep(REQUEST_DELAY)
    return None


def parse_hgnc_metadata(doc: dict) -> tuple[str, str, str, str]:
    """
    Extract (hgnc_id, approved_name, pipe-separated aliases, ensembl_gene_id) from a HGNC doc.
    """
    hgnc_id     = doc.get("hgnc_id", "")
    gene_name   = doc.get("name", "")
    ensembl_id  = doc.get("ensembl_gene_id", "")

    aliases: list[str] = []
    aliases += doc.get("alias_symbol", [])
    aliases += doc.get("prev_symbol",  [])
    aliases_str = " | ".join(aliases) if aliases else ""

    return hgnc_id, gene_name, aliases_str, ensembl_id


# ---------------------------------------------------------------------------
# Step 4: Fetch genomic coordinates from Ensembl
# ---------------------------------------------------------------------------

def fetch_ensembl_coords(symbol: str, grch37: bool = False) -> str:
    """
    Return 'chr:start-end' for *symbol* using Ensembl GRCh38 (default)
    or GRCh37 (if grch37=True).
    Returns empty string on failure.
    """
    url = (ENSEMBL37_URL if grch37 else ENSEMBL38_URL).format(symbol=symbol)
    try:
        resp = requests.get(url, headers={"Content-Type": "application/json"}, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            chrom = data.get("seq_region_name", "")
            start = data.get("start", "")
            end   = data.get("end", "")
            if chrom and start and end:
                prefix = "" if chrom.startswith("chr") else "chr"
                return f"{prefix}{chrom}:{start}-{end}"
    except Exception as exc:
        log.warning("Ensembl lookup failed for %s (grch37=%s): %s", symbol, grch37, exc)
    finally:
        time.sleep(REQUEST_DELAY)
    return ""


# ---------------------------------------------------------------------------
# Step 5: Extract disease context from paper text
# ---------------------------------------------------------------------------

def extract_disease_context(symbol: str, text: str) -> str:
    """
    Scan sentences containing *symbol* for disease-related keywords and
    return a concise context string (up to 3 unique snippets joined by ' | ').
    """
    disease_keywords = re.compile(
        r"\b(disease|disorder|syndrome|mutation|variant|pathogenic|phenotype|"
        r"hypomagnesemia|nephropathy|nephritis|cardiomyopathy|neuropathy|"
        r"myopathy|retinopathy|epilepsy|cancer|tumor|carcinoma|leukemia|"
        r"dysplasia|ataxia|dystrophy|deficiency|insufficiency|failure|"
        r"hematuria|proteinuria|glomerulopathy|glomerulosclerosis)\b",
        re.IGNORECASE,
    )

    # Split on sentence-ending punctuation (simple heuristic)
    sentences = re.split(r"(?<=[.!?])\s+", text)

    snippets: list[str] = []
    seen_snippets: set[str] = set()

    for sent in sentences:
        if re.search(rf"\b{re.escape(symbol)}\b", sent) and disease_keywords.search(sent):
            # Trim to a manageable length
            snippet = sent.strip()[:200]
            if snippet not in seen_snippets:
                seen_snippets.add(snippet)
                snippets.append(snippet)
            if len(snippets) >= 3:
                break

    return " | ".join(snippets)


# ---------------------------------------------------------------------------
# Step 5b: Fetch curated disease names from Open Targets
# ---------------------------------------------------------------------------

def fetch_diseases_from_opentargets(ensembl_id: str) -> list[str]:
    """
    Query the Open Targets GraphQL API for diseases associated with *ensembl_id*.

    Open Targets aggregates evidence from OMIM, Orphanet, ClinVar, and other
    sources without requiring an API key.  Returns disease name strings sorted
    by association score (highest first), e.g. ["nephrotic syndrome", ...].
    Returns an empty list on failure or when no associations are found.
    """
    if not ensembl_id:
        return []

    query = """
    {
      target(ensemblId: "%s") {
        associatedDiseases(page: {index: 0, size: 10}) {
          rows {
            disease { name }
            score
          }
        }
      }
    }
    """ % ensembl_id

    try:
        resp = requests.post(
            OPENTARGETS_URL,
            json={"query": query},
            headers={"Content-Type": "application/json"},
            timeout=15,
        )
        resp.raise_for_status()
        rows = (
            resp.json()
            .get("data", {})
            .get("target", {})
            .get("associatedDiseases", {})
            .get("rows", [])
        )
        return [row["disease"]["name"] for row in rows if row.get("disease", {}).get("name")]

    except Exception as exc:
        log.warning("Open Targets lookup failed for %s: %s", ensembl_id, exc)
        return []
    finally:
        time.sleep(REQUEST_DELAY)


def resolve_disease(symbol: str, paper_text: str, ensembl_id: str = "") -> str:
    """
    Combined disease resolution strategy:

    1. Fetch curated disease names from Open Targets (structured, from OMIM/Orphanet/ClinVar).
    2. Cross-reference each name against *paper_text* — keep only those
       that are actually mentioned in this paper (case-insensitive substring match).
    3. If any confirmed names remain, return them joined by ' | '.
    4. Fall back to sentence-level keyword extraction if no curated name
       is confirmed in the paper text.

    This gives structured disease names (e.g. "nephrotic syndrome") rather than
    raw sentence fragments, while ensuring the result is relevant to the paper.
    """
    curated = fetch_diseases_from_opentargets(ensembl_id)
    log.info("  Open Targets diseases for %s: %s", symbol, curated[:5])

    paper_lower = paper_text.lower()

    # Keep diseases whose name (or a meaningful stem of it) appears in the paper.
    # We match on the first two significant words to tolerate minor phrasing differences
    # e.g. "Alport Syndrome, Autosomal Recessive" → check for "alport syndrome".
    confirmed: list[str] = []
    for disease in curated:
        # OMIM names often include qualifiers after a comma or semicolon,
        # e.g. "Alport Syndrome, Autosomal Recessive".
        # Strip the qualifier first so the search key is clean, then try
        # matching with progressively fewer words (2 → 1) to maximise recall.
        base = re.split(r"[,;]", disease)[0].strip().lower()
        words = base.split()
        matched = any(
            " ".join(words[:n]) in paper_lower
            for n in range(min(2, len(words)), 0, -1)
        )
        if matched:
            confirmed.append(disease)

    if confirmed:
        log.info("  Confirmed in paper: %s", confirmed)
        return " | ".join(confirmed)

    # Fallback: sentence-level keyword extraction
    log.info("  No curated disease confirmed in paper for %s — falling back to context extraction", symbol)
    return extract_disease_context(symbol, paper_text)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def build_gene_records(
    text: str,
    min_genes: int = 5,
    max_genes: int = 20,
) -> list[GeneRecord]:
    """
    Full pipeline: candidates → HGNC validation → Ensembl coordinates
    → combined disease resolution (MyGene.info + paper cross-reference).

    Stops once *max_genes* valid records are collected (or all candidates
    are exhausted), guaranteeing at least *min_genes* if the paper contains them.
    """
    candidates = extract_gene_candidates(text)
    records: list[GeneRecord] = []

    for symbol in candidates:
        if len(records) >= max_genes:
            break

        log.info("Validating %s …", symbol)
        hgnc_doc = validate_with_hgnc(symbol)
        if hgnc_doc is None:
            log.debug("  %s not found in HGNC – skipping", symbol)
            continue

        hgnc_id, gene_name, aliases, ensembl_id = parse_hgnc_metadata(hgnc_doc)

        log.info("  Fetching coordinates for %s …", symbol)
        hg38 = fetch_ensembl_coords(symbol, grch37=False)
        hg19 = fetch_ensembl_coords(symbol, grch37=True)

        log.info("  Resolving disease for %s …", symbol)
        disease = resolve_disease(symbol, text, ensembl_id)

        records.append(GeneRecord(
            hgnc_id=hgnc_id,
            hgnc_gene_name=gene_name,
            gene_aliases=aliases,
            hg38_coordinates=hg38,
            hg19_coordinates=hg19,
            disease=disease,
        ))
        log.info("  ✓ %s → %s | disease: %s", symbol, gene_name, disease[:60])

    if len(records) < min_genes:
        log.warning(
            "Only %d valid genes found (minimum requested: %d). "
            "Consider expanding STOP_TOKENS or lowering min_genes.",
            len(records), min_genes,
        )

    return records


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def write_csv(records: list[GeneRecord], path: str = "gene_output.csv") -> None:
    fieldnames = list(asdict(records[0]).keys()) if records else []
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in records:
            writer.writerow(asdict(r))
    log.info("CSV written → %s", path)


def write_json(records: list[GeneRecord], path: str = "gene_output.json") -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump([asdict(r) for r in records], fh, indent=2)
    log.info("JSON written → %s", path)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    paper_text = fetch_paper_text()
    records    = build_gene_records(paper_text)

    write_csv(records,  "gene_output.csv")
    write_json(records, "gene_output.json")

    # Pretty-print summary
    print(f"\n{'='*70}")
    print(f"{'HGNC ID':<14} {'Symbol/Name':<20} {'hg38':<28} {'Disease'}")
    print(f"{'-'*70}")
    for r in records:
        name_col = r.hgnc_gene_name[:19]
        hg38_col = r.hg38_coordinates[:27]
        dis_col  = (r.disease[:40] + "…") if len(r.disease) > 40 else r.disease
        print(f"{r.hgnc_id:<14} {name_col:<20} {hg38_col:<28} {dis_col}")
    print(f"{'='*70}")
    print(f"\n{len(records)} gene records written to gene_output.csv / gene_output.json")


if __name__ == "__main__":
    main()