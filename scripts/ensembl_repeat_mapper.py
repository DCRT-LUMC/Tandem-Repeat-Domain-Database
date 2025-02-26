#!/usr/bin/env python3

import requests
import json
import time
import sys
import os
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from datetime import datetime

# Ensembl REST API base URL
ENSEMBL_API_BASE = "https://rest.ensembl.org"
MAX_RETRIES = 3
RATE_LIMIT_SLEEP = 0.1  # Sleep time between requests to avoid rate limiting

# Statistics tracking
stats = {
    "genes_processed": 0,
    "genes_with_repeats": 0,
    "total_repeats": 0,
    "repeats_with_coordinates": 0,
    "repeats_with_transcripts": 0,
    "total_transcript_overlaps": 0,
    "repeats_in_exons": 0,
    "repeats_in_introns": 0,
    "api_errors": 0,
    "api_rate_limits": 0,
    "start_time": None,
    "end_time": None
}

def get_overlapping_features(species, chrom, start, end, feature_type):
    """
    Get genomic features that overlap with a region using the Ensembl API.
    """
    # Format chromosome name for Ensembl (remove 'chr' prefix if present)
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    
    # Prepare API endpoint
    region = f"{chrom}:{start}-{end}"
    endpoint = f"{ENSEMBL_API_BASE}/overlap/region/{species}/{region}"
    params = {"feature": feature_type, "content-type": "application/json"}
    
    # Try the request with retries
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(endpoint, params=params, headers={"Content-Type": "application/json"})
            
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:  # Too Many Requests
                wait_time = int(response.headers.get('Retry-After', RATE_LIMIT_SLEEP * 10))
                print(f"Rate limited. Waiting {wait_time} seconds.")
                stats["api_rate_limits"] += 1
                time.sleep(wait_time)
            else:
                print(f"Error for {region}: HTTP {response.status_code}")
                stats["api_errors"] += 1
                time.sleep(RATE_LIMIT_SLEEP * (attempt + 1))
        except Exception as e:
            print(f"Request exception for {region}: {str(e)}")
            stats["api_errors"] += 1
            time.sleep(RATE_LIMIT_SLEEP * (attempt + 1))
    
    return []  # Return empty list if all attempts fail

def is_repeat_in_exon(transcript_id, chrom, start, end):
    """
    Determine if a repeat region falls within an exon of the specified transcript.
    """
    endpoint = f"{ENSEMBL_API_BASE}/lookup/id/{transcript_id}"
    params = {"expand": 1, "content-type": "application/json"}
    
    for attempt in range(MAX_RETRIES):
        try:
            response = requests.get(endpoint, params=params, headers={"Content-Type": "application/json"})
            
            if response.status_code == 200:
                transcript_data = response.json()
                
                # Check if the region overlaps with any exons
                if "Exon" in transcript_data:
                    for exon in transcript_data["Exon"]:
                        exon_start = exon["start"]
                        exon_end = exon["end"]
                        
                        # Check for overlap
                        if not (end < exon_start or start > exon_end):
                            return True, {"exon_id": exon["id"], "exon_start": exon_start, "exon_end": exon_end}
                
                return False, {}
            
            elif response.status_code == 429:  # Too Many Requests
                wait_time = int(response.headers.get('Retry-After', RATE_LIMIT_SLEEP * 10))
                print(f"Rate limited. Waiting {wait_time} seconds.")
                stats["api_rate_limits"] += 1
                time.sleep(wait_time)
            else:
                stats["api_errors"] += 1
                time.sleep(RATE_LIMIT_SLEEP * (attempt + 1))
        except Exception as e:
            print(f"Error checking exons for {transcript_id}: {str(e)}")
            stats["api_errors"] += 1
            time.sleep(RATE_LIMIT_SLEEP * (attempt + 1))
    
    return False, {}

def process_repeat(repeat, known_transcripts):
    """Process a single repeat to find overlapping transcripts."""
    if not all(key in repeat for key in ['chrom', 'chromStart', 'chromEnd']):
        return repeat
    
    chrom = repeat['chrom']
    start = repeat['chromStart']
    end = repeat['chromEnd']
    
    # Get overlapping transcripts
    transcripts = get_overlapping_features("human", chrom, start, end, "transcript")
    
    # Process each transcript
    transcript_info = []
    has_transcript_overlaps = False
    for transcript in transcripts:
        transcript_id = transcript["id"]
        is_known = transcript_id in known_transcripts
        
        # Check if the repeat is in an exon of this transcript
        in_exon, exon_details = is_repeat_in_exon(transcript_id, chrom, start, end)
        
        transcript_info.append({
            "transcript_id": transcript_id,
            "gene_id": transcript.get("Parent", ""),
            "biotype": transcript.get("biotype", ""),
            "is_known_transcript": is_known,
            "in_exon": in_exon,
            "exon_details": exon_details if in_exon else {}
        })
        
        has_transcript_overlaps = True
        stats["total_transcript_overlaps"] += 1
        if in_exon:
            stats["repeats_in_exons"] += 1
        else:
            stats["repeats_in_introns"] += 1
        
        time.sleep(RATE_LIMIT_SLEEP)  # Respect API rate limits
    
    # Add transcript information to the repeat
    repeat["transcript_overlaps"] = transcript_info
    if has_transcript_overlaps:
        stats["repeats_with_transcripts"] += 1
    
    return repeat

def process_gene_data(gene_data):
    """Process a single gene entry."""
    gene_name = gene_data.get('gene', 'Unknown')
    
    # Skip genes without repeats
    if 'repeats' not in gene_data or not gene_data['repeats']:
        return gene_data
    
    # Update stats for this gene
    stats["genes_with_repeats"] += 1
    stats["total_repeats"] += len(gene_data['repeats'])
    
    # Get known transcripts for this gene
    known_transcripts = set(gene_data.get('transcripts', []))
    
    # Process repeats that have coordinates
    valid_repeats = [r for r in gene_data['repeats'] if all(key in r for key in ['chrom', 'chromStart', 'chromEnd'])]
    stats["repeats_with_coordinates"] += len(valid_repeats)
    
    # Process repeats with a thread pool for efficiency
    with ThreadPoolExecutor(max_workers=5) as executor:
        processed_repeats = list(executor.map(
            lambda r: process_repeat(r, known_transcripts),
            valid_repeats
        ))
    
    # Update the valid repeats with transcript information
    valid_repeat_idx = 0
    for i, repeat in enumerate(gene_data['repeats']):
        if all(key in repeat for key in ['chrom', 'chromStart', 'chromEnd']):
            gene_data['repeats'][i] = processed_repeats[valid_repeat_idx]
            valid_repeat_idx += 1
    
    return gene_data

def print_summary_log():
    """Print a summary of all operations performed."""
    duration = stats["end_time"] - stats["start_time"]
    hours, remainder = divmod(duration.total_seconds(), 3600)
    minutes, seconds = divmod(remainder, 60)
    
    print("\n" + "="*80)
    print(" ENSEMBL REPEAT MAPPER - EXECUTION SUMMARY ")
    print("="*80)
    print(f"Start time: {stats['start_time'].strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"End time:   {stats['end_time'].strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Duration:   {int(hours)}h {int(minutes)}m {int(seconds)}s")
    print("-"*80)
    print("PROCESSING STATISTICS:")
    print(f"  • Total genes processed:          {stats['genes_processed']}")
    print(f"  • Genes with repeats:             {stats['genes_with_repeats']} ({stats['genes_with_repeats']/stats['genes_processed']*100:.1f}% of total)")
    print(f"  • Total repeats found:            {stats['total_repeats']}")
    print(f"  • Repeats with coordinates:       {stats['repeats_with_coordinates']} ({stats['repeats_with_coordinates']/stats['total_repeats']*100:.1f}% of total)")
    print(f"  • Repeats with transcript data:   {stats['repeats_with_transcripts']} ({stats['repeats_with_transcripts']/stats['repeats_with_coordinates']*100:.1f}% of those with coordinates)")
    print(f"  • Total transcript overlaps:      {stats['total_transcript_overlaps']}")
    print("-"*80)
    print("EXONIC vs INTRONIC REPEATS:")
    total_located = stats['repeats_in_exons'] + stats['repeats_in_introns']
    if total_located > 0:
        print(f"  • Repeats in exonic regions:      {stats['repeats_in_exons']} ({stats['repeats_in_exons']/total_located*100:.1f}%)")
        print(f"  • Repeats in intronic regions:    {stats['repeats_in_introns']} ({stats['repeats_in_introns']/total_located*100:.1f}%)")
    print("-"*80)
    print("API INTERACTIONS:")
    print(f"  • API rate limit hits:            {stats['api_rate_limits']}")
    print(f"  • API errors encountered:         {stats['api_errors']}")
    print("="*80)

def main(input_file, output_file, limit=None):
    """
    Main function to process repeats and find overlapping transcripts.
    
    Args:
        input_file: Path to input JSON file
        output_file: Path to output JSON file
        limit: Optional limit on number of genes to process (for testing)
    """
    stats["start_time"] = datetime.now()
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found.")
        sys.exit(1)
    
    # Load data
    try:
        with open(input_file, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Could not parse {input_file} as JSON.")
        sys.exit(1)
    
    # Apply limit if specified for testing
    if limit and isinstance(limit, int) and limit > 0:
        data = data[:limit]
        print(f"TESTING MODE: Processing only the first {limit} entries")
    
    # Process each gene one at a time (to avoid overwhelming the API)
    processed_data = []
    
    # Create progress bar
    with tqdm(total=len(data), desc="Processing genes", unit="gene") as pbar:
        for gene_idx, gene_data in enumerate(data):
            gene_name = gene_data.get('gene', f"Gene #{gene_idx+1}")
            pbar.set_description(f"Processing {gene_name}")
            
            processed_gene = process_gene_data(gene_data)
            processed_data.append(processed_gene)
            stats["genes_processed"] += 1
            
            # Save progress periodically
            if (gene_idx + 1) % 10 == 0:
                with open(f"{output_file}.partial", 'w') as f:
                    json.dump(processed_data, f, indent=2)
                pbar.set_postfix(saved="partial")
            
            pbar.update(1)
    
    # Write final output
    with open(output_file, 'w') as f:
        json.dump(processed_data, f, indent=2)
    
    stats["end_time"] = datetime.now()
    print_summary_log()
    print(f"Processing complete. Results saved to {output_file}")

if __name__ == "__main__":
    # Check arguments - now accepts an optional limit parameter
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python ensembl_repeat_mapper.py input_file.json output_file.json [limit]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Process optional limit parameter
    limit = None
    if len(sys.argv) == 4:
        try:
            limit = int(sys.argv[3])
            if limit <= 0:
                print("Error: limit must be a positive integer")
                sys.exit(1)
        except ValueError:
            print("Error: limit must be an integer")
            sys.exit(1)
    
    main(input_file, output_file, limit)