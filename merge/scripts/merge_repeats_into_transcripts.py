import json
import os
from collections import defaultdict
from tqdm import tqdm

def merge_repeat_and_transcript_data():
    # Define file paths
    repeats_file = 'data/hg38_repeats.json'
    transcripts_file = 'data/ensembl_transcripts_merged.json'
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
    
    # Keep track of which repeat entries have been matched
    repeat_entry_matched = {i: False for i in range(len(repeat_data))}
    
    for i, entry in enumerate(tqdm(repeat_data, desc="Processing repeat entries")):
        if "position" in entry and "protein" in entry.get("position", ""):
            # Extract UniProt ID from position field if available (e.g., "amino acids 343-389 on protein Q6TDP4")
            parts = entry.get("position", "").split()
            for part in parts:
                if len(part) == 6 and part[0] in "OPQ" and part[1:].isalnum():
                    entry["uniProtId"] = part
                    break
        
        if "uniProtId" in entry and entry["uniProtId"]:
            uniprot_to_repeats[entry["uniProtId"]].append((i, entry))
        
        # Also map by geneName if available
        if "geneName" in entry and entry["geneName"]:
            gene_name = entry["geneName"].upper()  # Normalize to uppercase for matching
            gene_to_repeats[gene_name].append((i, entry))
    
    # Track statistics
    total_repeat_entries = len(repeat_data)
    matches_by_uniprot = 0
    matches_by_gene = 0
    matches_by_alias = 0
    transcripts_with_repeats = 0
    transcripts_without_repeats = 0
    unmatched_repeats = 0
    
    # Perform the actual merging
    merged_entries = {}  # Use a dict to avoid duplicates
    
    # First, process all transcript entries and find matching repeat entries
    for transcript_entry in tqdm(transcript_data, desc="Merging entries"):
        gene_name = transcript_entry.get("gene", "").upper()
        uniprot_id = transcript_entry.get("uniprot_id", "")
        merged_key = f"gene_{gene_name}"  # Use gene as key for uniqueness
        
        if gene_name:
            # Start with the transcript entry
            if merged_key not in merged_entries:
                merged_entries[merged_key] = dict(transcript_entry)
                merged_entries[merged_key]["repeats"] = []  # Add a list to hold repeat data
            
            has_repeats = False
            
            # PRIORITY 1: Try to find matching repeat entries by UniProt ID first
            if uniprot_id and uniprot_id in uniprot_to_repeats:
                has_repeats = True
                matches_by_uniprot += 1
                for idx, repeat_entry in uniprot_to_repeats[uniprot_id]:
                    repeat_info = {k: v for k, v in repeat_entry.items()}
                    merged_entries[merged_key]["repeats"].append(repeat_info)
                    repeat_entry_matched[idx] = True
            
            # PRIORITY 2: Try to find matching repeat entries by gene name
            if gene_name in gene_to_repeats:
                for idx, repeat_entry in gene_to_repeats[gene_name]:
                    if not repeat_entry_matched[idx]:  # Only add if not already added via UniProt
                        has_repeats = True
                        matches_by_gene += 1
                        repeat_info = {k: v for k, v in repeat_entry.items()}
                        merged_entries[merged_key]["repeats"].append(repeat_info)
                        repeat_entry_matched[idx] = True
            
            # PRIORITY 3: Check aliases for matches
            if "aliases" in transcript_entry:
                for alias in transcript_entry["aliases"]:
                    alias = alias.upper()
                    if alias in gene_to_repeats and alias != gene_name:
                        for idx, repeat_entry in gene_to_repeats[alias]:
                            if not repeat_entry_matched[idx]:
                                has_repeats = True
                                matches_by_alias += 1
                                repeat_info = {k: v for k, v in repeat_entry.items()}
                                merged_entries[merged_key]["repeats"].append(repeat_info)
                                repeat_entry_matched[idx] = True
            
            # Count entries with/without repeats
            if has_repeats:
                transcripts_with_repeats += 1
            else:
                transcripts_without_repeats += 1
    
    # Now add any unmatched repeat entries as standalone entries
    print("Adding unmatched repeat entries...")
    
    for i, entry in enumerate(repeat_data):
        if not repeat_entry_matched[i]:
            unmatched_repeats += 1
            # Try to extract a meaningful key
            if "geneName" in entry and entry["geneName"]:
                gene_name = entry["geneName"].upper()
            elif "uniProtId" in entry and entry["uniProtId"]:
                gene_name = f"PROT_{entry['uniProtId']}"
            else:
                gene_name = f"UNNAMED_REPEAT_{i}"
            
            merged_key = f"unmatched_{gene_name}_{i}"
            
            # Create a new entry
            merged_entries[merged_key] = {
                "gene": gene_name,
                "ensembl_id": None,
                "uniprot_id": entry.get("uniProtId", ""),
                "aliases": [],
                "transcripts": [],
                "repeats": [entry]  # Add the unmatched repeat entry
            }
    
    # Convert the dictionary to a list for output
    merged_data = list(merged_entries.values())
    
    # Write the merged data to the output file
    with open(output_file, 'w') as f:
        json.dump(merged_data, f, indent=4)
    
    # Print detailed merge report
    print("\n" + "="*50)
    print("MERGE REPORT")
    print("="*50)
    print(f"Input Data:")
    print(f"- Transcript entries: {len(transcript_data)}")
    print(f"- Repeat entries: {total_repeat_entries}")
    print(f"\nMatching Results:")
    print(f"- Transcripts with repeats: {transcripts_with_repeats} ({transcripts_with_repeats/len(transcript_data)*100:.1f}%)")
    print(f"- Transcripts without repeats: {transcripts_without_repeats} ({transcripts_without_repeats/len(transcript_data)*100:.1f}%)")
    print(f"\nMatching Method:")
    print(f"- Matched by UniProt ID: {matches_by_uniprot}")
    print(f"- Matched by gene name: {matches_by_gene}")
    print(f"- Matched by alias: {matches_by_alias}")
    print(f"- Unmatched repeat entries: {unmatched_repeats} ({unmatched_repeats/total_repeat_entries*100:.1f}%)")
    print(f"\nOutput Data:")
    print(f"- Total merged entries: {len(merged_data)}")
    print(f"- Entries with repeats: {sum(1 for entry in merged_data if entry.get('repeats'))}")
    print("="*50)
    
    print(f"\nMerged data written to {output_file}")

# Run the merge function
if __name__ == "__main__":
    merge_repeat_and_transcript_data()