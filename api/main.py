from fastapi import FastAPI
import json
from pathlib import Path

app = FastAPI(title="TReXomeDB API")

# Load the JSON data
DATA_PATH = Path(__file__).resolve().parent.parent / "docs" / "all_annotated_repeats.json"

with open(DATA_PATH) as f:
    repeats_data = json.load(f)


@app.get("/")
def root():
    return {"message": "TReXomeDB API is running"}


@app.get("/api/repeats/raw")
def get_raw_repeats():
    return {
        "count": len(repeats_data),
        "results": repeats_data[:10]  # return only first 10 for testing
    }
