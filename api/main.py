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
                "chrom_start": item.get("chromStart"),
                "chrom_end": item.get("chromEnd"),
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
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
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

        if chrom and safe_lower(row["chrom"]) != safe_lower(chrom):
            continue

        # Overlap logic for genomic interval queries
        row_start = row.get("chrom_start")
        row_end = row.get("chrom_end")

        if start is not None and row_end is not None and row_end < start:
            continue

        if end is not None and row_start is not None and row_start > end:
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

@app.get("/api/proteins")
def get_proteins(
    gene: str | None = None,
    uniprot_id: str | None = None,
    repeat_type: str | None = None,
    status: str | None = None,
    eligibility: str | None = None,
    single_repeat_exon: bool = False,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=200),
):
    proteins = {}

    for row in ALL_REPEAT_EXON_MATCHES:
        key = row["uniprot_id"]
        if not key:
            continue

        if key not in proteins:
            proteins[key] = {
                "gene": row["gene"],
                "uniprot_id": row["uniprot_id"],
                "repeat_type": row["repeat_type"],
                "status": row["status"],
                "chrom": row["chrom"],
                "eligibility": row["eligibility"],
                "repeat_count": 0,
                "transcript_ids": set(),
                "exon_ids": set(),
                "exon_numbers": set(),
                "has_in_frame_exons": False,
                "has_fully_coding_exons": False,
                "has_single_repeat_exon": False,
                "source_repeats": []
            }

        protein = proteins[key]
        protein["transcript_ids"].add(row["transcript_id"])
        protein["exon_ids"].add(row["exon_id"])
        protein["exon_numbers"].add(row["exon_number"])

        if row["frame_status"] == "in_frame":
            protein["has_in_frame_exons"] = True

        if row["coding_status"] == "fully_coding":
            protein["has_fully_coding_exons"] = True

        exon_key = (row["transcript_id"], row["exon_id"])
        if exon_key in SINGLE_REPEAT_EXON_KEYS:
            protein["has_single_repeat_exon"] = True

        source_repeat = row["source"]
        repeat_identifier = (
            source_repeat.get("uniProtId"),
            source_repeat.get("protein_start"),
            source_repeat.get("protein_end"),
            source_repeat.get("position"),
        )

        if "_seen_repeats" not in protein:
            protein["_seen_repeats"] = set()

        if repeat_identifier not in protein["_seen_repeats"]:
            protein["_seen_repeats"].add(repeat_identifier)
            protein["repeat_count"] += 1
            protein["source_repeats"].append(source_repeat)

    results = []
    for protein in proteins.values():
        protein.pop("_seen_repeats", None)

        protein["transcript_ids"] = sorted([x for x in protein["transcript_ids"] if x is not None])
        protein["exon_ids"] = sorted([x for x in protein["exon_ids"] if x is not None])
        protein["exon_numbers"] = sorted([x for x in protein["exon_numbers"] if x is not None])

        if gene and safe_lower(protein["gene"]) != safe_lower(gene):
            continue

        if uniprot_id and safe_lower(protein["uniprot_id"]) != safe_lower(uniprot_id):
            continue

        if repeat_type and safe_lower(protein["repeat_type"]) != safe_lower(repeat_type):
            continue

        if status and safe_lower(protein["status"]) != safe_lower(status):
            continue

        if eligibility and safe_lower(protein["eligibility"]) != safe_lower(eligibility):
            continue

        if single_repeat_exon and not protein["has_single_repeat_exon"]:
            continue

        results.append(protein)

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



@app.get("/api/proteins")
def get_proteins(
    gene: str | None = None,
    uniprot_id: str | None = None,
    repeat_type: str | None = None,
    status: str | None = None,
    eligibility: str | None = None,
    single_repeat_exon: bool = False,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=200),
):
    proteins = {}

    for row in ALL_REPEAT_EXON_MATCHES:
        key = row["uniprot_id"]
        if not key:
            continue

        if key not in proteins:
            proteins[key] = {
                "gene": row["gene"],
                "uniprot_id": row["uniprot_id"],
                "repeat_type": row["repeat_type"],
                "status": row["status"],
                "chrom": row["chrom"],
                "eligibility": row["eligibility"],
                "repeat_count": 0,
                "transcript_ids": set(),
                "exon_ids": set(),
                "exon_numbers": set(),
                "has_in_frame_exons": False,
                "has_fully_coding_exons": False,
                "has_single_repeat_exon": False,
                "source_repeats": []
            }

        protein = proteins[key]
        protein["transcript_ids"].add(row["transcript_id"])
        protein["exon_ids"].add(row["exon_id"])
        protein["exon_numbers"].add(row["exon_number"])

        if row["frame_status"] == "in_frame":
            protein["has_in_frame_exons"] = True

        if row["coding_status"] == "fully_coding":
            protein["has_fully_coding_exons"] = True

        exon_key = (row["transcript_id"], row["exon_id"])
        if exon_key in SINGLE_REPEAT_EXON_KEYS:
            protein["has_single_repeat_exon"] = True

        source_repeat = row["source"]
        repeat_identifier = (
            source_repeat.get("uniProtId"),
            source_repeat.get("protein_start"),
            source_repeat.get("protein_end"),
            source_repeat.get("position"),
        )

        if "_seen_repeats" not in protein:
            protein["_seen_repeats"] = set()

        if repeat_identifier not in protein["_seen_repeats"]:
            protein["_seen_repeats"].add(repeat_identifier)
            protein["repeat_count"] += 1
            protein["source_repeats"].append(source_repeat)

    results = []
    for protein in proteins.values():
        protein.pop("_seen_repeats", None)

        protein["transcript_ids"] = sorted([x for x in protein["transcript_ids"] if x is not None])
        protein["exon_ids"] = sorted([x for x in protein["exon_ids"] if x is not None])
        protein["exon_numbers"] = sorted([x for x in protein["exon_numbers"] if x is not None])

        if gene and safe_lower(protein["gene"]) != safe_lower(gene):
            continue

        if uniprot_id and safe_lower(protein["uniprot_id"]) != safe_lower(uniprot_id):
            continue

        if repeat_type and safe_lower(protein["repeat_type"]) != safe_lower(repeat_type):
            continue

        if status and safe_lower(protein["status"]) != safe_lower(status):
            continue

        if eligibility and safe_lower(protein["eligibility"]) != safe_lower(eligibility):
            continue

        if single_repeat_exon and not protein["has_single_repeat_exon"]:
            continue

        results.append(protein)

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
