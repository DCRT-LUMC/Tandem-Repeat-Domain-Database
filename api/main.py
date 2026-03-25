from fastapi import FastAPI, Query
import json
import csv
import math
from pathlib import Path

app = FastAPI(title="TReXomeDB API")

# -----------------------------
# Load source files
# -----------------------------
BASE_DIR = Path(__file__).resolve().parent.parent

DATA_PATH = BASE_DIR / "docs" / "all_annotated_repeats.json"
CSV_PATH = BASE_DIR / "docs" / "DCRT_Gene_list_August_2024.csv"

with open(DATA_PATH, encoding="utf-8") as f:
    repeats_data = json.load(f)

gene_eligibility = {}
with open(CSV_PATH, encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gene = (row.get("Entity Name") or "").strip()
        eligibility = (row.get("Eligibility") or "").strip()
        if gene:
            gene_eligibility[gene] = eligibility if eligibility else "N/A"


# -----------------------------
# Helpers
# -----------------------------
def safe_lower(value):
    if value is None:
        return ""
    return str(value).strip().lower()


def parse_bool(value):
    if value is None:
        return False
    return str(value).strip().lower() in {"true", "1", "yes"}


def flatten_repeat_exon_matches():
    """
    Create one row per repeat–exon match.
    """
    results = []

    for item in repeats_data:
        gene = item.get("geneName", "Unknown")
        exon_info = item.get("ensembl_exon_info", {})
        transcripts = exon_info.get("transcripts", [])

        for transcript in transcripts:
            transcript_id = transcript.get("transcript_id")
            containing_exons = transcript.get("containing_exons", [])

            for exon in containing_exons:
                results.append({
                    "gene": gene,
                    "uniprot_id": item.get("uniProtId"),
                    "repeat_type": item.get("repeatType"),
                    "status": item.get("status"),
                    "chrom": item.get("chrom"),
                    "protein_start": item.get("protein_start"),
                    "protein_end": item.get("protein_end"),
                    "block_count": item.get("blockCount"),
                    "eligibility": gene_eligibility.get(gene, "N/A"),
                    "transcript_id": transcript_id,
                    "exon_id": exon.get("exon_id"),
                    "exon_number": exon.get("exon_number"),
                    "frame_status": exon.get("frame_status"),
                    "coding_status": exon.get("coding_status"),
                    "exon_start": exon.get("exon_start"),
                    "exon_end": exon.get("exon_end"),
                    "source": item
                })

    return results


def build_single_repeat_exon_keys():
    """
    Returns a set of (transcript_id, exon_id) pairs for exons that are:
    - in-frame
    - fully coding
    - contain exactly one repeat
    """
    exon_counts = {}

    for item in repeats_data:
        exon_info = item.get("ensembl_exon_info", {})
        transcripts = exon_info.get("transcripts", [])

        for transcript in transcripts:
            transcript_id = transcript.get("transcript_id")
            containing_exons = transcript.get("containing_exons", [])

            for exon in containing_exons:
                exon_id = exon.get("exon_id")
                key = (transcript_id, exon_id)

                if key not in exon_counts:
                    exon_counts[key] = {
                        "count": 0,
                        "frame_status": exon.get("frame_status"),
                        "coding_status": exon.get("coding_status"),
                    }

                exon_counts[key]["count"] += 1

    valid_keys = set()
    for key, info in exon_counts.items():
        if (
            info["count"] == 1
            and info["frame_status"] == "in_frame"
            and info["coding_status"] == "fully_coding"
        ):
            valid_keys.add(key)

    return valid_keys


ALL_REPEAT_EXON_MATCHES = flatten_repeat_exon_matches()
SINGLE_REPEAT_EXON_KEYS = build_single_repeat_exon_keys()


# -----------------------------
# Endpoints
# -----------------------------
@app.get("/")
def root():
    return {"message": "TReXomeDB API is running"}


@app.get("/api/repeats/raw")
def get_raw_repeats():
    return {
        "count": len(repeats_data),
        "results": repeats_data[:10]
    }


@app.get("/api/repeats")
def get_repeats(
    gene: str | None = None,
    uniprot_id: str | None = None,
    transcript_id: str | None = None,
    exon_id: str | None = None,
    exon_number: str | None = None,
    single_repeat_exon: bool = False,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=200),
):
    results = []

    for row in ALL_REPEAT_EXON_MATCHES:
        if gene and safe_lower(row["gene"]) != safe_lower(gene):
            continue

        if uniprot_id and safe_lower(row["uniprot_id"]) != safe_lower(uniprot_id):
            continue

        if transcript_id and safe_lower(row["transcript_id"]) != safe_lower(transcript_id):
            continue

        if exon_id and safe_lower(row["exon_id"]) != safe_lower(exon_id):
            continue

        if exon_number is not None and str(row["exon_number"]) != str(exon_number):
            continue

        if single_repeat_exon:
            key = (row["transcript_id"], row["exon_id"])
            if key not in SINGLE_REPEAT_EXON_KEYS:
                continue

        results.append(row)

    total_results = len(results)
    total_pages = math.ceil(total_results / page_size) if total_results > 0 else 0

    start_idx = (page - 1) * page_size
    end_idx = start_idx + page_size
    paginated_results = results[start_idx:end_idx]

    return {
        "page": page,
        "page_size": page_size,
        "total_results": total_results,
        "total_pages": total_pages,
        "results": paginated_results
    }
