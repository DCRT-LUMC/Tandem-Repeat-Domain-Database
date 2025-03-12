import json
import time
import os
import requests
from requests.adapters import HTTPAdapter, Retry
from tqdm import tqdm

# Constants
API_URL = "https://rest.uniprot.org"

# Setup session with proper retry handling
retries = Retry(
    total=5,
    backoff_factor=0.25,
    status_forcelist=[429, 500, 502, 503, 504],
    respect_retry_after_header=True
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    """Check if the API response is valid"""
    try:
        response.raise_for_status()
        return True
    except requests.HTTPError as e:
        print(f"API Error: {e}")
        print(f"Response: {response.text}")
        return False


def get_uniprot_gene_info(uniprot_id):
    """Fetch gene name and aliases from UniProt API for a given UniProt ID"""
    if not uniprot_id:
        return None, None
    
    url = f"{API_URL}/uniprotkb/{uniprot_id}.json"
    
    try:
        response = session.get(url)
        if not check_response(response):
            return None, None
        
        data = response.json()
        
        # Extract gene names
        gene_info = {}
        if "genes" in data:
            for gene in data["genes"]:
                if "geneName" in gene:
                    gene_info["primary"] = gene["geneName"].get("value", "")
                
                if "synonyms" in gene:
                    gene_info["synonyms"] = [s["value"] for s in gene["synonyms"]]
        
        primary_name = gene_info.get("primary", "")
        synonyms = gene_info.get("synonyms", [])
        
        return primary_name, synonyms
    
    except Exception as e:
        print(f"Exception when fetching {uniprot_id}: {str(e)}")
        return None, None


def process_entry(entry):
    """Process a single JSON entry to update gene names"""
    uniprot_id = entry.get("uniProtId", "")
    if uniprot_id:
        primary_name, synonyms = get_uniprot_gene_info(uniprot_id)
        
        # Update entry with gene information
        if primary_name:
            entry["geneName"] = primary_name
        
        # Store synonyms as an array
        if synonyms:
            entry["geneName2"] = synonyms
    
    return entry


def main():
    input_file = os.path.join("merge_test_small", "data", "hg38_repeats_100.json")
    output_file = os.path.join("merge_test_small", "data", "gname_hg38_repeats.json")
    
    # Ensure the input file exists
    if not os.path.exists(input_file):
        print(f"Input file {input_file} not found.")
        return
    
    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    total_entries = len(data)
    print(f"Loaded {total_entries} entries. Starting API queries...")
    
    # Process entries sequentially with a progress bar
    updated_data = []
    progress_bar = tqdm(data, desc="Processing entries")
    
    for entry in progress_bar:
        updated_entry = process_entry(entry)
        updated_data.append(updated_entry)
        # Small delay to avoid overwhelming the API
        time.sleep(0.1)
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    print(f"Updated data written to {output_file}")


if __name__ == "__main__":
    main()