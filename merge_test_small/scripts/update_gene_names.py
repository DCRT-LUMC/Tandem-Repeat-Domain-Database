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

# Cache for UniProt gene information
uniprot_cache = {}
cache_hits = 0
api_calls = 0

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
    global cache_hits, api_calls
    
    if not uniprot_id:
        return None, None
    
    # Check if we have this UniProt ID in cache
    if uniprot_id in uniprot_cache:
        cache_hits += 1
        return uniprot_cache[uniprot_id]
    
    url = f"{API_URL}/uniprotkb/{uniprot_id}.json"
    api_calls += 1
    
    try:
        response = session.get(url)
        if not check_response(response):
            uniprot_cache[uniprot_id] = (None, None)
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
        
        # Store in cache for future use
        result = (primary_name, synonyms)
        uniprot_cache[uniprot_id] = result
        return result
    
    except Exception as e:
        print(f"Exception when fetching {uniprot_id}: {str(e)}")
        uniprot_cache[uniprot_id] = (None, None)
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
            entry["aliases"] = synonyms
    
    return entry


def main():
    input_file = os.path.join("merge_test_small", "data", "1000_length_filtered_hg38_repeats.json")
    output_file = os.path.join("merge_test_small", "data", "1000_gname_hg38_repeats.json")
    
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
        # Small delay only when making actual API calls
        if api_calls > 0 and api_calls % 10 == 0:
            time.sleep(0.1)
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    # Report cache statistics
    unique_proteins = len(uniprot_cache)
    print(f"\nCache statistics:")
    print(f"Total entries processed: {total_entries}")
    print(f"Unique UniProt IDs: {unique_proteins}")
    print(f"API calls made: {api_calls}")
    print(f"Cache hits: {cache_hits}")
    if api_calls > 0:
        print(f"API call reduction: {cache_hits/(api_calls+cache_hits)*100:.2f}%")
    
    print(f"\nUpdated data written to {output_file}")


if __name__ == "__main__":
    main()