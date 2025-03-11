import json
import time
import os
import requests
from requests.adapters import HTTPAdapter, Retry
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# Constants
API_URL = "https://rest.uniprot.org"
DEFAULT_BATCH_SIZE = 100
MAX_WORKERS = 5

# Setup session with proper retry handling
retries = Retry(
    total=5,
    backoff_factor=0.5,
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


def process_batch(batch):
    """Process a batch of entries with proper rate limiting"""
    results = []
    # Add tqdm progress bar for entries within a batch
    for entry in tqdm(batch, desc="Processing entries", leave=False):
        results.append(process_entry(entry))
        # Small delay to avoid overwhelming the API
        time.sleep(0.1)
    return results


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
    
    # Split data into batches for parallel processing
    batch_size = min(DEFAULT_BATCH_SIZE, total_entries)
    batches = [data[i:i+batch_size] for i in range(0, total_entries, batch_size)]
    
    # Process batches in parallel
    updated_data = []
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_batch = {executor.submit(process_batch, batch): i for i, batch in enumerate(batches)}
        
        # Create progress bar for overall batch processing
        progress_bar = tqdm(total=len(batches), desc="Processing batches")
        
        for future in as_completed(future_to_batch):
            batch_index = future_to_batch[future]
            try:
                batch_results = future.result()
                updated_data.extend(batch_results)
                
                # Update progress bar
                progress_bar.update(1)
                progress_bar.set_postfix(batch=f"{batch_index+1}/{len(batches)}", entries=f"{len(updated_data)}/{total_entries}")
                
            except Exception as e:
                print(f"Batch {batch_index} generated an exception: {e}")
        
        progress_bar.close()
    
    # Ensure the output is in the same order as the input
    if len(updated_data) != total_entries:
        print(f"Warning: Output size ({len(updated_data)}) doesn't match input size ({total_entries})")
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    print(f"Updated data written to {output_file}")


if __name__ == "__main__":
    main()