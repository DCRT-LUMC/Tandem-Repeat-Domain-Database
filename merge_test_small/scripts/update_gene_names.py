import json
import requests
import time
import os
from concurrent.futures import ThreadPoolExecutor

def get_uniprot_gene_info(uniprot_id):
    """Fetch gene name and aliases from UniProt API for a given UniProt ID"""
    if not uniprot_id:
        return None, None
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Error fetching {uniprot_id}: {response.status_code}")
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
        
        # Store synonyms as an array instead of comma-delimited string
        if synonyms:
            entry["geneName2"] = synonyms
        
        # Add delay to avoid hitting API rate limits
        time.sleep(0.5)
    
    return entry

def main():
    # Path to your input JSON file
    input_file = "merge_test_small/data/hg38_repeats_100.json"
    output_file = "merge_test_small/data/gname_hg38_repeats_100.json"
    
    # Ensure the input file exists
    if not os.path.exists(input_file):
        print(f"Input file {input_file} not found.")
        return
    
    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    print(f"Loaded {len(data)} entries. Starting API queries...")
    
    # Process entries in parallel with a thread pool
    updated_data = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        for i, updated_entry in enumerate(executor.map(process_entry, data)):
            updated_data.append(updated_entry)
            if (i + 1) % 10 == 0:
                print(f"Processed {i + 1}/{len(data)} entries")
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    print(f"Updated data written to {output_file}")

if __name__ == "__main__":
    main()