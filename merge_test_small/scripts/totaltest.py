#!/usr/bin/env python
# filepath: /C:/Users/ojfab/Documents/GitHub/Tandem-Repeat-Domain-Database/merge_test_small/scripts/repeat_database.py
import json
import requests
import time
import os
import sys
import logging
import datetime
import argparse
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor

# Set up logging
def setup_logging():
    """Set up logging configuration"""
    # Determine project root directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    # Create logs directory if it doesn't exist
    logs_dir = os.path.join(project_root, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    
    # Setup log file
    log_filename = f"repeat_db_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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
    return log_filepath

# Global stats dictionary for tracking API usage
api_stats = {
    "requests": 0,
    "rate_limits": 0,
    "errors": 0
}

# Create a cache to avoid redundant API calls
query_cache = {}

######################
# GENE NAME FUNCTIONS
######################

def get_uniprot_gene_info(uniprot_id):
    """Fetch gene name and aliases from UniProt API for a given UniProt ID"""
    if not uniprot_id:
        return None, None
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    try:
        logging.debug(f"Querying UniProt for {uniprot_id}")
        response = requests.get(url)
        if response.status_code != 200:
            logging.warning(f"Error fetching {uniprot_id}: {response.status_code}")
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
        logging.error(f"Exception when fetching {uniprot_id}: {str(e)}")
        return None, None

def process_entry_gn(entry):
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

def update_gene_names(input_file, output_file):
    """Update gene names in the repeat data using UniProt API"""
    logging.info("Starting gene name update process")
    
    # Ensure the input file exists
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} not found.")
        return False
    
    # Load JSON data
    with open(input_file, 'r') as f:
        data = json.load(f)
    
    logging.info(f"Loaded {len(data)} entries. Starting UniProt API queries...")
    
    # Process entries in parallel with a thread pool
    updated_data = []
    with ThreadPoolExecutor(max_workers=5) as executor:
        for i, updated_entry in enumerate(executor.map(process_entry_gn, data)):
            updated_data.append(updated_entry)
            if (i + 1) % 10 == 0:
                logging.info(f"Processed {i + 1}/{len(data)} entries")
    
    # Write updated JSON data
    with open(output_file, 'w') as f:
        json.dump(updated_data, f, indent=2)
    
    logging.info(f"Gene names updated and saved to {output_file}")
    return True

######################
# REPEAT LENGTH FUNCTIONS
######################

def add_repeat_length_and_filter(input_file, output_file, min_length=60):
    """
    Adds a repeatLength field to each repeat entry and filters out repeats 
    shorter than the minimum length.
    """
    logging.info(f"Adding repeat length and filtering repeats shorter than {min_length} bp")
    
    # Load the JSON data
    with open(input_file, 'r') as f:
        repeats = json.load(f)
    
    logging.info(f"Processing {len(repeats)} repeat entries...")
    
    # Create a new list for filtered repeats
    filtered_repeats = []
    excluded_count = 0
    
    # Loop through each repeat entry
    for repeat in tqdm(repeats):
        # Calculate repeat length
        chrom_start = repeat.get("chromStart", 0)
        chrom_end = repeat.get("chromEnd", 0)
        repeat_length = chrom_end - chrom_start
        
        # Skip repeats shorter than min_length
        if repeat_length < min_length:
            excluded_count += 1
            continue
        
        # Create a new dictionary with the repeatLength field inserted after "comments"
        new_repeat = {}
        for key, value in repeat.items():
            new_repeat[key] = value
            # Insert repeatLength right after comments
            if key == "comments":
                new_repeat["repeatLength"] = repeat_length
        
        # If "comments" wasn't found, add it at the end
        if "repeatLength" not in new_repeat:
            new_repeat["repeatLength"] = repeat_length
        
        # Add the modified repeat to our filtered list
        filtered_repeats.append(new_repeat)
    
    # Save the modified and filtered JSON
    with open(output_file, 'w') as f:
        json.dump(filtered_repeats, f, indent=2)
    
    logging.info(f"Added repeatLength field to entries and excluded {excluded_count} repeats shorter than {min_length} bp.")
    logging.info(f"Saved {len(filtered_repeats)} repeats to {output_file}")
    return True

######################
# EXON INFO FUNCTIONS
######################

class EnsemblRestClient(object):
    """
    Client for the Ensembl REST API with proper rate limiting
    """
    def __init__(self, server='https://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        # Build URL with parameters
        url = self.server + endpoint
        
        data = None

        # Check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                wait_time = 1 - delta
                api_stats["rate_limits"] += 1
                logging.info(f"Rate limited: waiting {wait_time:.2f}s before requesting {endpoint}")
                time.sleep(wait_time)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            # Log the request
            api_stats["requests"] += 1
            logging.debug(f"API request: {endpoint}")
            
            # Use requests instead of urllib for consistency with the rest of your code
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

def get_ensembl_info(chrom, start, end, species="human"):
    """
    Get transcript and exon information using the Ensembl API.
    """
    # Format chromosome correctly (Ensembl doesn't use "chr" prefix)
    chrom_id = chrom.replace("chr", "")
    
    # Create cache key
    cache_key = f"{chrom_id}:{start}-{end}"
    if cache_key in query_cache:
        return query_cache[cache_key]
    
    # Initialize Ensembl REST client
    client = EnsemblRestClient()
    headers = {"Content-Type": "application/json"}
    
    result = {"transcripts": [], "exons": []}
    
    # Get transcripts overlapping the region
    transcripts = client.perform_rest_action(
        endpoint=f"/overlap/region/{species}/{chrom_id}:{start}-{end}",
        hdrs=headers,
        params={'feature': 'transcript'}
    )
    
    if transcripts:
        result["transcripts"] = transcripts
    
    # Get exons overlapping the region
    exons = client.perform_rest_action(
        endpoint=f"/overlap/region/{species}/{chrom_id}:{start}-{end}",
        hdrs=headers,
        params={'feature': 'exon'}
    )
    
    if exons:
        result["exons"] = exons
    
    # For each transcript, get detailed info including all its exons
    transcript_details = {}
    
    for transcript in result["transcripts"]:
        try:
            transcript_id = transcript["id"]
            
            # Get detailed transcript information with all exons
            transcript_detail = client.perform_rest_action(
                endpoint=f"/lookup/id/{transcript_id}",
                hdrs=headers,
                params={'expand': 1}
            )
            
            if transcript_detail:
                transcript_details[transcript_id] = transcript_detail
                
        except Exception as e:
            logging.error(f"Error querying Ensembl transcript details for {transcript_id}: {str(e)}")
    
    result["transcript_details"] = transcript_details
    
    # Store in cache
    query_cache[cache_key] = result
    return result

def is_canonical_transcript(transcript, all_transcripts):
    """Determine if transcript is canonical based on MANE Select or longest CDS"""
    # Check if transcript has MANE Select tag - this indicates the canonical transcript
    if "Tags" in transcript and "MANE_Select" in transcript.get("Tags", []):
        return True
        
    # If no MANE Select, check if it's marked as canonical
    if transcript.get("is_canonical", 0) == 1:
        return True
        
    # If still can't determine, use longest CDS as fallback
    gene_id = transcript.get("Parent")
    if not gene_id:
        return False
        
    gene_transcripts = [t for t in all_transcripts if t.get("Parent") == gene_id]
    
    # If there's only one transcript for this gene, it's canonical by default
    if len(gene_transcripts) == 1:
        return True
        
    # Find transcript with longest translation
    try:
        # Use transcript length as proxy for CDS length if translation_length not available
        translation_lengths = []
        for t in gene_transcripts:
            if "Translation" in t:
                translation_lengths.append((t["id"], len(t.get("Translation", {}).get("seq", ""))))
            else:
                translation_lengths.append((t["id"], t.get("end", 0) - t.get("start", 0)))
                
        if not translation_lengths:
            return False
            
        longest_transcript = max(translation_lengths, key=lambda x: x[1])
        return transcript["id"] == longest_transcript[0]
    except Exception:
        return False

def classify_repeat_location(repeat, transcript):
    """Classify if the repeat is exonic, intronic, or outside the transcript"""
    repeat_start = int(repeat["chromStart"])
    repeat_end = int(repeat["chromEnd"])
    
    tx_start = int(transcript.get("start", 0))
    tx_end = int(transcript.get("end", 0))
    
    # Check if completely outside the transcript
    if repeat_end <= tx_start or repeat_start >= tx_end:
        return "outside"
    
    # Get exons from the transcript
    exons = transcript.get("Exon", [])
    
    # Check if overlapping any exon
    for exon in exons:
        exon_start = int(exon.get("start", 0))
        exon_end = int(exon.get("end", 0))
        
        if max(repeat_start, exon_start) < min(repeat_end, exon_end):
            return "exonic"
    
    # If not exonic but within transcript bounds, must be intronic
    return "intronic"

def get_coding_status(exon, transcript, repeat_start, repeat_end):
    """Determine if exon is coding, non-coding, or partial and if it's in UTR"""
    exon_start = int(exon.get("start", 0))
    exon_end = int(exon.get("end", 0))
    
    # Get coding region if available
    translation = transcript.get("Translation", {})
    has_translation = bool(translation)
    
    if not has_translation:
        return "non_coding", "non_coding_transcript", 0
    
    # Get CDS start/end
    cds_start = int(translation.get("start", 0))
    cds_end = int(translation.get("end", 0))
    
    # Completely non-coding exon
    if exon_end <= cds_start or exon_start >= cds_end:
        # Determine if it's 5' or 3' UTR based on strand
        strand = 1 if transcript.get("strand") == 1 else -1
        if exon_end <= cds_start:
            utr_type = "5'_UTR" if strand == 1 else "3'_UTR"
        else:
            utr_type = "3'_UTR" if strand == 1 else "5'_UTR"
        return "non_coding", utr_type, 0
    
    # Partially coding exon
    if exon_start < cds_start or exon_end > cds_end:
        utr_regions = []
        strand = 1 if transcript.get("strand") == 1 else -1
        if exon_start < cds_start:
            utr_regions.append("5'_UTR" if strand == 1 else "3'_UTR")
        if exon_end > cds_end:
            utr_regions.append("3'_UTR" if strand == 1 else "5'_UTR")
        
        # Calculate coding percentage
        exon_length = exon_end - exon_start
        coding_start = max(exon_start, cds_start)
        coding_end = min(exon_end, cds_end)
        coding_length = coding_end - coding_start
        coding_percentage = (coding_length / exon_length) * 100
        
        return "partial_coding", "+".join(utr_regions), round(coding_percentage, 2)
    
    # Fully coding exon
    return "fully_coding", "none", 100

def process_exon_info(input_file, output_file, limit=None):
    """
    Process the repeat data JSON and add exon information using Ensembl API.
    """
    logging.info("Starting exon information annotation")
    
    # Load repeat data
    with open(input_file, 'r') as f:
        repeats = json.load(f)
    
    # Filter out entries that don't have proper coordinate data
    valid_repeats = [r for r in repeats if "chrom" in r and "chromStart" in r and "chromEnd" in r]
    
    # Apply limit if specified
    if limit and isinstance(limit, int) and limit > 0:
        valid_repeats = valid_repeats[:limit]
        logging.info(f"Processing first {limit} out of {len(repeats)} repeats...")
    else:
        logging.info(f"Processing {len(valid_repeats)} out of {len(repeats)} repeats with valid coordinates...")
    
    # Process each repeat
    for repeat_idx, repeat in enumerate(tqdm(valid_repeats)):
        # Save intermediate results every 10 repeats to avoid losing progress
        if repeat_idx > 0 and repeat_idx % 10 == 0:
            with open(output_file + ".temp", 'w') as f:
                json.dump(repeats, f, indent=2)
        
        chrom = repeat["chrom"]
        start = int(repeat["chromStart"])
        end = int(repeat["chromEnd"])

        # Get the repeat's strand
        repeat_strand = repeat.get("strand", "")
        # Convert to Ensembl format for comparison
        expected_ensembl_strand = 1 if repeat_strand == "+" else -1 if repeat_strand == "-" else None
        
        # Get transcript and exon information from Ensembl
        api_data = get_ensembl_info(chrom, start, end)
        
        if not api_data or not api_data["transcripts"]:
            repeat["ensembl_exon_info"] = {
                "transcripts_count": 0,
                "has_canonical_transcript": False,
                "location_summary": "unknown",
                "transcripts": []
            }
            continue
        
        all_transcripts = []
        transcript_details = api_data["transcript_details"]
        
        # Use transcript details for comprehensive information
        for transcript_id, transcript in transcript_details.items():
            all_transcripts.append(transcript)
            
        transcript_info = []
        locations = set()
        
        for transcript in all_transcripts:
            try:
                # Skip transcripts with different strand if repeat strand is specified
                if expected_ensembl_strand is not None and transcript.get("strand") != expected_ensembl_strand:
                    continue

                # Check if this is likely the canonical transcript
                is_canonical = is_canonical_transcript(transcript, all_transcripts)
                
                # Classify location (exonic, intronic, outside)
                location = classify_repeat_location(repeat, transcript)
                locations.add(location)
                
                # Get basic transcript info
                transcript_id = transcript["id"]
                
                # Find exons in this transcript
                exons = transcript.get("Exon", [])
                exon_count = len(exons)
                
                # Sort exons by genomic order
                strand = transcript.get("strand")
                if strand == 1:  # "+" strand
                    # For + strand, exons are ordered from 5' to 3'
                    exons = sorted(exons, key=lambda e: e.get("start", 0))
                else:  # "-" strand
                    # For - strand, exons are ordered from 3' to 5'
                    exons = sorted(exons, key=lambda e: e.get("start", 0), reverse=True)
                
                containing_exons = []
                for i, exon in enumerate(exons):
                    exon_start = int(exon.get("start", 0))
                    exon_end = int(exon.get("end", 0))
                    
                    # Check if the repeat overlaps this exon
                    if max(start, exon_start) < min(end, exon_end):  # Overlap
                        exon_number = i + 1  # 1-based exon numbering
                        
                        overlap_start = max(start, exon_start)
                        overlap_end = min(end, exon_end)
                        overlap_length = overlap_end - overlap_start
                        exon_length = exon_end - exon_start
                        overlap_percentage = (overlap_length / exon_length) * 100
                        
                        # Determine position in transcript
                        if exon_count == 1:
                            position = "single_exon"
                        elif i == 0:
                            position = "first_exon"
                        elif i == exon_count - 1:
                            position = "last_exon"
                        else:
                            position = f"middle_exon_{exon_number}"
                        
                        # Get coding status
                        coding_status, utr_status, coding_percentage = get_coding_status(exon, transcript, start, end)
                        
                        exon_info = {
                            "exon_id": exon.get("id", ""),
                            "exon_number": exon_number,
                            "position": position,
                            "coding_status": coding_status,
                            "overlap_bp": overlap_length,
                            "overlap_percentage": round(overlap_percentage, 2),
                            "utr_status": utr_status,
                            "coding_percentage": coding_percentage
                        }
                        
                        containing_exons.append(exon_info)
                
                # Get transcript biotype
                biotype = transcript.get("biotype", "unknown")
                
                transcript_info.append({
                    "transcript_id": transcript_id,
                    "transcript_name": transcript.get("display_name", ""),
                    "is_canonical": is_canonical,
                    "biotype": biotype,
                    "location": location,
                    "exon_count": exon_count,
                    "containing_exons": containing_exons
                })
            except Exception as e:
                logging.error(f"Error processing transcript {transcript.get('id', 'unknown')}: {str(e)}")
                continue
        
        # Summarize location (prioritize exonic > intronic > outside/unknown)
        location_summary = "unknown"
        if "exonic" in locations:
            location_summary = "exonic"
        elif "intronic" in locations:
            location_summary = "intronic"
        elif "outside" in locations:
            location_summary = "intergenic"
        
        has_canonical = any(t["is_canonical"] for t in transcript_info) if transcript_info else False
        
        # Add exon information to the repeat
        repeat["ensembl_exon_info"] = {
            "transcripts_count": len(transcript_info),
            "has_canonical_transcript": has_canonical,
            "location_summary": location_summary,
            "transcripts": transcript_info
        }
    
    # Save updated repeat data
    with open(output_file, 'w') as f:
        json.dump(repeats, f, indent=2)
    
    # Remove temp file if exists
    if os.path.exists(output_file + ".temp"):
        os.remove(output_file + ".temp")
    
    logging.info(f"Updated repeat data with exon information saved to {output_file}")
    return True

######################
# MAIN PROGRAM
######################

def main():
    print(f"Command line arguments: {sys.argv}")
    """Main entry point for the command line tool"""
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory (project root)
    project_root = os.path.dirname(script_dir)
    # Default data directory
    data_dir = os.path.join(project_root, "data")
    
    parser = argparse.ArgumentParser(
        description="Tandem Repeat Domain Database processing tool",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # Make input and output optional with defaults in data directory
    parser.add_argument("-i", "--input", dest="input_file",
                        default=os.path.join(data_dir, "hg38_repeats_100.json"),
                        help="Input JSON file with repeat data (default: data/hg38_repeats_100.json)")
    
    parser.add_argument("-o", "--output", dest="output_file",
                        default=None,  # We'll auto-generate this based on operations
                        help="Output file to save processed data (default: auto-generated in data folder)")
    
    parser.add_argument("-gn", "--gene-names", action="store_true", 
                        help="Update gene names using UniProt API")
    
    parser.add_argument("-rl", "--repeat-length", type=int, metavar="MIN_LENGTH",
                        help="Add repeatLength field and filter repeats shorter than MIN_LENGTH")
    
    parser.add_argument("-e", "--exon-info", action="store_true",
                        help="Add exon information using Ensembl API")
    
    parser.add_argument("-l", "--limit", type=int, default=None,
                        help="Limit the number of repeats to process (useful for testing)")
    
    args = parser.parse_args()
    
    # Check if input file exists, if not, try to find it in the data directory
    if not os.path.exists(args.input_file):
        alternative_path = os.path.join(data_dir, os.path.basename(args.input_file))
        if os.path.exists(alternative_path):
            logging.info(f"Input file not found at {args.input_file}, using {alternative_path} instead")
            args.input_file = alternative_path
    
    # Generate sensible output filename if none provided
    if args.output_file is None:
        base_name = os.path.splitext(os.path.basename(args.input_file))[0]
        parts = [base_name]
        
        if args.gene_names:
            parts.append("gname")
            
        if args.repeat_length:
            parts.append(f"length{args.repeat_length}")
            
        if args.exon_info:
            parts.append("exons")
            
        if args.limit:
            parts.append(f"limit{args.limit}")
            
        output_filename = "_".join(parts) + ".json"
        args.output_file = os.path.join(data_dir, output_filename)
        
    logging.info(f"Input file: {args.input_file}")
    logging.info(f"Output will be saved to: {args.output_file}")
    
    # Create pipeline based on requested operations
    pipeline = []
    temp_files = []
    current_input = args.input_file
    
    if args.gene_names:
        pipeline.append(("Update gene names", update_gene_names))
        temp_files.append(os.path.join(data_dir, "temp_gene_names.json"))
    
    if args.repeat_length:
        pipeline.append((f"Add repeat length (min {args.repeat_length} bp)", 
                        lambda i, o: add_repeat_length_and_filter(i, o, args.repeat_length)))
        temp_files.append("temp_repeat_length.json")
    
    if args.exon_info:
        pipeline.append(("Add exon information", 
                        lambda i, o: process_exon_info(i, o, args.limit)))
        temp_files.append("temp_exon_info.json")
    
    # Execute pipeline
    for i, (description, operation) in enumerate(pipeline):
        logging.info(f"Step {i+1}/{len(pipeline)}: {description}")
        output = args.output_file if i == len(pipeline)-1 else temp_files[i]
        success = operation(current_input, output)
        if not success:
            logging.error(f"Step {i+1} failed. Stopping pipeline.")
            return
        current_input = output

if __name__ == "__main__":
    # Setup logging
    log_file = setup_logging()
    
    try:
        main()
    except Exception as e:
        logging.exception(f"Unhandled exception: {e}")
    finally:
        # Print summary statistics if any API calls were made
        if api_stats["requests"] > 0:
            logging.info("=== Processing Summary ===")
            logging.info(f"Total API requests: {api_stats['requests']}")
            logging.info(f"Rate limits encountered: {api_stats['rate_limits']}")
            logging.info(f"Errors: {api_stats['errors']}")
            logging.info("========================")