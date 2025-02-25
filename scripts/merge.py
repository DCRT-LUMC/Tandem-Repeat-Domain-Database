import json
import os
from collections import defaultdict
from tqdm import tqdm

def merge_repeat_and_transcript_data():
    # Define file paths
    repeats_file = 'data/hg38_repeats.json'
    transcripts_file = 'data/ensembl_transcripts.json'
    output_file = 'data/merged_repeats_transcripts.json'
    
    # Load the data
    with open(repeats_file, 'r') as f:
        repeat_data = json.load(f)
    
    with open(transcripts_file, 'r') as f:
        transcript_data = json.load(f)
    
    print(f"Loaded {len(repeat_data)} repeat entries and {len(transcript_data)} transcript entries")
    
    # Create mappings for repeat data by various identifiers
    uniprot_to_repeats = defaultdict(list)
    gene_to_repeats = defaultdict(list)
    
    for entry in tqdm(repeat_data, desc="Processing repeat entries"):
        if "uniProtId" in entry and entry["uniProtId"]:
            uniprot_to_repeats[entry["uniProtId"]].append(entry)
        
        # Also map by geneName if available
        if "geneName" in entry and entry["geneName"]:
            gene_name = entry["geneName"].upper()  # Normalize to uppercase for matching
            gene_to_repeats[gene_name].append(entry)
    
    # Create mappings for transcript data
    gene_to_transcripts = defaultdict(list)
    alias_to_transcript = defaultdict(list)
    
    for entry in tqdm(transcript_data, desc="Processing transcript entries"):
        if "gene" in entry and entry["gene"]:
            gene_name = entry["gene"].upper()  # Normalize to uppercase
            gene_to_transcripts[gene_name].append(entry)
        
        # Also map by aliases
        if "aliases" in entry and entry["aliases"]:
            for alias in entry["aliases"]:
                alias = alias.upper()  # Normalize to uppercase
                alias_to_transcript[alias].append(entry)
    
    # Perform the actual merging
    merged_entries = {}  # Use a dict to avoid duplicates
    
    # First, process all transcript entries and find matching repeat entries
    for transcript_entry in tqdm(transcript_data, desc="Merging entries"):
        gene_name = transcript_entry.get("gene", "").upper()
        merged_key = f"gene_{gene_name}"  # Use gene as key for uniqueness
        
        if gene_name:
            # Start with the transcript entry
            if merged_key not in merged_entries:
                merged_entries[merged_key] = dict(transcript_entry)
                merged_entries[merged_key]["repeats"] = []  # Add a list to hold repeat data
            
            # Try to find matching repeat entries by gene name
            if gene_name in gene_to_repeats:
                for repeat_entry in gene_to_repeats[gene_name]:
                    # Add repeat info to the repeats list
                    repeat_info = {k: v for k, v in repeat_entry.items() if k not in ["gene", "aliases", "transcripts"]}
                    merged_entries[merged_key]["repeats"].append(repeat_info)
            
            # Also check aliases for matches
            if "aliases" in transcript_entry:
                for alias in transcript_entry["aliases"]:
                    alias = alias.upper()
                    if alias in gene_to_repeats and alias != gene_name:
                        for repeat_entry in gene_to_repeats[alias]:
                            repeat_info = {k: v for k, v in repeat_entry.items() if k not in ["gene", "aliases", "transcripts"]}
                            merged_entries[merged_key]["repeats"].append(repeat_info)
    
    # Now find repeat entries that didn't match any transcript
    print("Finding unmatched repeat entries...")
    for gene_name, repeat_list in gene_to_repeats.items():
        merged_key = f"gene_{gene_name}"
        
        # If we haven't already created an entry for this gene
        if merged_key not in merged_entries:
            # Create a new entry
            base_entry = {
                "gene": gene_name,
                "ensembl_id": None,
                "aliases": [],
                "transcripts": [],
                "repeats": []
            }
            
            # Add repeat info to the repeats list
            for repeat_entry in repeat_list:
                repeat_info = {k: v for k, v in repeat_entry.items() if k not in ["gene", "aliases", "transcripts"]}
                base_entry["repeats"].append(repeat_info)
            
            merged_entries[merged_key] = base_entry
    
    # Convert the dictionary to a list for output
    merged_data = list(merged_entries.values())
    
    # Write the merged data to the output file
    with open(output_file, 'w') as f:
        json.dump(merged_data, f, indent=4)
    
    print(f"Merged data written to {output_file}")
    print(f"Total entries: {len(merged_data)}")
    print(f"Entries with repeat data: {sum(1 for entry in merged_data if entry.get('repeats'))}")

# Run the merge function
merge_repeat_and_transcript_data()