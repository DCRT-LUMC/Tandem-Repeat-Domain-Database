import json
import time
import os
import requests
from requests.adapters import HTTPAdapter, Retry
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging
import datetime
import threading

# Constants
API_URL = "https://rest.ensembl.org"
DEFAULT_BATCH_SIZE = 100
MAX_WORKERS = 15  # Increased from 5 to maximize throughput
REQUESTS_PER_SEC = 15  # Ensembl API rate limit

# Setup logging
logs_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "logs")
os.makedirs(logs_dir, exist_ok=True)
log_filename = f"ensembl_gene_names_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
log_filepath = os.path.join(logs_dir, log_filename)

file_handler = logging.FileHandler(log_filepath)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

console_handler = logging.StreamHandler()
console_handler.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))  # No timestamp

logging.basicConfig(
    level=logging.INFO,
    handlers=[file_handler, console_handler]
)

logging.info(f"Logging to file: {log_filepath}")

# Global stats dictionary for tracking API usage
api_stats = {
    "requests": 0,
    "rate_limits": 0,
    "errors": 0
}

# Create a cache to avoid redundant API calls
query_cache = {}

class EnsemblRestClient(object):
    """
    Client for the Ensembl REST API with proper rate limiting
    - Implements the strategy from the ensembl_example_client.py
    """
    def __init__(self, server=API_URL, reqs_per_sec=REQUESTS_PER_SEC):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0
        self.req_lock = threading.RLock()  # Add thread lock for synchronization

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        # Build URL with parameters
        url = self.server + endpoint
        
        data = None

        # Use lock to synchronize across threads
        with self.req_lock:
            # Check if we need to rate limit ourselves
            if self.req_count >= self.reqs_per_sec:
                delta = time.time() - self.last_req
                if delta < 1:
                    wait_time = 1 - delta
                    api_stats["rate_limits"] += 1
                    logging.debug(f"Rate limited: waiting {wait_time:.2f}s before requesting {endpoint}")
                    time.sleep(wait_time)
                self.last_req = time.time()
                self.req_count = 0
            
            try:
                # Log the request
                api_stats["requests"] += 1
                logging.debug(f"API request: {endpoint}")
                
                # Use requests with timeout
                response = requests.get(url, headers=hdrs, params=params, timeout=15)
                
                if response.status_code == 200:
                    data = response.json()
                    self.req_count += 1
                    
                # Check if we are being rate limited by the server
                elif response.status_code == 429:
                    if 'Retry-After' in response.headers:
                        retry = response.headers['Retry-After']
                        api_stats["rate_limits"] += 1
                        logging.warning(f"Server rate limit hit: waiting {retry}s before retrying {endpoint}")
                        time.sleep(float(retry))
                        return self.perform_rest_action(endpoint, hdrs, params)
                else:
                    api_stats["errors"] += 1
                    logging.error(f"Request failed: {endpoint} (Status {response.status_code})")
                    
            except Exception as e:
                api_stats["errors"] += 1
                logging.error(f"Request error: {endpoint} - {str(e)}")
            
        return data

# Create a single global client instance
ensembl_client = EnsemblRestClient()

def get_ensembl_gene_info(chrom, start, end, uniprot_id=None, species="human"):
    """
    Get gene name and synonyms using the Ensembl API.
    First checks by genomic location, and if that doesn't work, tries by UniProt ID.
    """
    # Format chromosome correctly (Ensembl doesn't use "chr" prefix)
    chrom_id = chrom.replace("chr", "")
    
    # Create cache key
    cache_key = f"{chrom_id}:{start}-{end}"
    if cache_key in query_cache:
        return query_cache[cache_key]
    
    # Use the global client instance instead of creating a new one
    headers = {"Content-Type": "application/json"}
    
    result = {"gene_name": None, "synonyms": []}
    
    # Get genes overlapping the region
    genes = ensembl_client.perform_rest_action(
        endpoint=f"/overlap/region/{species}/{chrom_id}:{start}-{end}",
        hdrs=headers,
        params={'feature': 'gene'}
    )
    
    # If we found genes, get the gene name and synonyms
    if genes and len(genes) > 0:
        # Sort genes by biotype priority (protein_coding takes precedence)
        genes.sort(key=lambda g: 0 if g.get("biotype", "") == "protein_coding" else 1)
        
        # Get the first gene (highest priority)
        gene_id = genes[0].get("id")
        
        if gene_id:
            # Get detailed gene information
            gene_details = ensembl_client.perform_rest_action(
                endpoint=f"/lookup/id/{gene_id}",
                hdrs=headers,
                params={'expand': 1}
            )
            
            if gene_details:
                result["gene_name"] = gene_details.get("display_name", "")
                
                # Get gene synonyms
                xrefs = ensembl_client.perform_rest_action(
                    endpoint=f"/xrefs/id/{gene_id}",
                    hdrs=headers,
                )
                
                if xrefs:
                    synonyms = []
                    for xref in xrefs:
                        # Get all possible synonym sources
                        if "display_id" in xref and xref["display_id"] != result["gene_name"]:
                            synonyms.append(xref["display_id"])
                        if "synonym" in xref and isinstance(xref["synonym"], list):
                            synonyms.extend(xref["synonym"])
                        elif "synonym" in xref and xref["synonym"]:
                            synonyms.append(xref["synonym"])
                    result["synonyms"] = list(set(synonyms))  # Remove duplicates
    
    # If we still don't have a gene name and we have a UniProt ID, try looking up by UniProt
    if (not result["gene_name"] or not result["synonyms"]) and uniprot_id:
        try:
            # Try to get Ensembl gene ID from UniProt ID
            uniprot_lookup = ensembl_client.perform_rest_action(
                endpoint=f"/xrefs/symbol/{species}/{uniprot_id}",
                hdrs=headers
            )
            
            if not uniprot_lookup or len(uniprot_lookup) == 0:
                # If direct lookup fails, try searching by external references
                uniprot_lookup = ensembl_client.perform_rest_action(
                    endpoint=f"/xrefs/name/{species}/{uniprot_id}",
                    hdrs=headers
                )
            
            if uniprot_lookup and len(uniprot_lookup) > 0:
                gene_id = uniprot_lookup[0].get("id")
                
                if gene_id:
                    # Get detailed gene information
                    gene_details = ensembl_client.perform_rest_action(
                        endpoint=f"/lookup/id/{gene_id}",
                        hdrs=headers,
                        params={'expand': 1}
                    )
                    
                    if gene_details:
                        result["gene_name"] = result["gene_name"] or gene_details.get("display_name", "")
                        
                        # Get additional gene synonyms
                        xrefs = ensembl_client.perform_rest_action(
                            endpoint=f"/xrefs/id/{gene_id}",
                            hdrs=headers,
                        )
                        
                        if xrefs:
                            new_synonyms = set(result["synonyms"])
                            for xref in xrefs:
                                if "display_id" in xref and xref["display_id"] != result["gene_name"]:
                                    new_synonyms.add(xref["display_id"])
                                if "synonym" in xref and isinstance(xref["synonym"], list):
                                    new_synonyms.update(xref["synonym"])
                                elif "synonym" in xref and xref["synonym"]:
                                    new_synonyms.add(xref["synonym"])
                            result["synonyms"] = list(new_synonyms)
        except Exception as e:
            logging.error(f"Error looking up UniProt ID {uniprot_id}: {str(e)}")
    
    # Store in cache
    query_cache[cache_key] = result
    return result

def process_entry(entry):
    """Process a single JSON entry to update gene names"""
    chrom = entry.get("chrom", "")
    start = int(entry.get("chromStart", 0))
    end = int(entry.get("chromEnd", 0))
    uniprot_id = entry.get("uniProtId", "")
    
    if chrom and start and end:
        gene_info = get_ensembl_gene_info(chrom, start, end, uniprot_id)
        
        # Update entry with gene information
        if gene_info["gene_name"]:
            entry["geneName"] = gene_info["gene_name"]
        
        # Store synonyms as an array
        if gene_info["synonyms"]:
            entry["geneName2"] = gene_info["synonyms"]
    
    return entry

def main():
    input_file = os.path.join("merge_test_small", "data", "hg38_repeats_100.json")
    output_file = os.path.join("merge_test_small", "data", "ensembl_gname_hg38_repeats.json")
    
    # Ensure the input file exists
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return
    
    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    total_entries = len(data)
    logging.info(f"Loaded {total_entries} entries. Starting API queries...")
    
    # Start timing BEFORE processing
    start_time = time.time()
    
    # Process entries directly in parallel without batching
    updated_data = [None] * total_entries  # Pre-allocate list to maintain order
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # Submit each entry as a separate job
        future_to_index = {executor.submit(process_entry, entry): i for i, entry in enumerate(data)}
        
        # Create a single progress bar for all entries
        with tqdm(total=total_entries, desc="Processing entries") as progress_bar:
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    # Store result at the same index to maintain order
                    updated_data[index] = future.result()
                    progress_bar.update(1)
                except Exception as e:
                    logging.error(f"Entry {index} generated an exception: {e}")
                    # Keep original entry on error
                    updated_data[index] = data[index]
                    progress_bar.update(1)
    
    # Ensure all entries were processed
    if None in updated_data:
        logging.warning("Some entries were not processed")
        # Fill any missing entries with original data
        for i, entry in enumerate(updated_data):
            if entry is None:
                updated_data[i] = data[i]
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    # Calculate duration AFTER processing
    duration = time.time() - start_time
    
    logging.info(f"Updated data written to {output_file}")
    logging.info("=== Processing Summary ===")
    logging.info(f"Total API requests: {api_stats['requests']}")
    logging.info(f"Rate limits encountered: {api_stats['rate_limits']}")
    logging.info(f"Errors: {api_stats['errors']}")
    logging.info(f"Total runtime: {duration:.2f} seconds")
    logging.info("========================")

if __name__ == "__main__":
    main()
