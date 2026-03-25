from fastapi import FastAPI
import json
import csv
from pathlib import Path

app = FastAPI(title="TReXomeDB API")

# Load the JSON data
DATA_PATH = Path(__file__).resolve().parent.parent / "docs" / "all_annotated_repeats.json"

with open(DATA_PATH) as f:
    repeats_data = json.load(f)

# Load gene eligibility CSV
CSV_PATH = Path(__file__).resolve().parent.parent / "docs" / "DCRT_Gene_list_August_2024.csv"

gene_eligibility = {}

with open(CSV_PATH) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gene = row.get("Entity Name")
        eligibility = row.get("Eligibility")

        if gene and eligibility:
            gene_eligibility[gene] = eligibility


@app.get("/")
def root():
    return {"message": "TReXomeDB API is running"}


@app.get("/api/repeats/raw")
def get_raw_repeats():
    return {
        "count": len(repeats_data),
        "results": repeats_data[:10]
    }
