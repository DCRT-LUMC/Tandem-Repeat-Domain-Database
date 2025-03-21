import json
import requests
import time
from tqdm import tqdm
import os
import sys
import logging
import datetime

# First, determine project root directory
script_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(script_dir)

# Create logs directory if it doesn't exist
logs_dir = os.path.join(project_root, "logs")
os.makedirs(logs_dir, exist_ok=True)

# Simple logging setup with logs in the logs directory
log_filename = f"ensembl_api_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
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

def convert_to_zero_based(data):
    """
    Convert Ensembl 1-based coordinates to 0-based coordinates.
    This function can be applied to a single object or a list of objects.
    """
    if isinstance(data, list):
        for item in data:
            convert_to_zero_based(item)
        return data
    
    # Skip if not a dictionary
    if not isinstance(data, dict):
        return data
    
    # Convert start coordinate from 1-based to 0-based
    if "start" in data:
        data["start"] = data["start"] - 1
    
    # Process nested objects
    for key, value in data.items():
        if isinstance(value, (dict, list)):
            convert_to_zero_based(value)
    
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
        # Convert coordinates to 0-based
        result["transcripts"] = convert_to_zero_based(transcripts)
    
    # Get exons overlapping the region
    exons = client.perform_rest_action(
        endpoint=f"/overlap/region/{species}/{chrom_id}:{start}-{end}",
        hdrs=headers,
        params={'feature': 'exon'}
    )
    
    if exons:
        # Convert coordinates to 0-based
        result["exons"] = convert_to_zero_based(exons)
    
    # For each transcript, get detailed info including all its exons
    transcript_details = {}
    
    for transcript in result["transcripts"]:
        try:
            transcript_id = transcript["id"]
            
            # Get detailed transcript information with all exons and version information
            transcript_detail = client.perform_rest_action(
                endpoint=f"/lookup/id/{transcript_id}",
                hdrs=headers,
                params={'expand': 1, 'format': 'full'}
            )
            
            # Convert all coordinates to 0-based
            if transcript_detail:
                convert_to_zero_based(transcript_detail)
            
            # Construct the versioned ID
            if transcript_detail and "version" in transcript_detail:
                versioned_transcript_id = f"{transcript_id}.{transcript_detail['version']}"
            else:
                versioned_transcript_id = transcript_id
            
            if transcript_detail:
                transcript_details[transcript_id] = transcript_detail
                
        except Exception as e:
            print(f"Error querying Ensembl transcript details for {transcript_id}: {e}")
    
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

def process_repeat_data(repeat_data_file, output_file, limit=None):
    """
    Process the repeat data JSON and add exon information using Ensembl API.
    
    Parameters:
        repeat_data_file: Input JSON file with repeat data
        output_file: Output file to save the updated data
        limit: Optional. If set, process only this many entries
    """
    
    # Load repeat data
    with open(repeat_data_file, 'r') as f:
        repeats = json.load(f)
    
    # Filter out entries that don't have proper coordinate data
    valid_repeats = [r for r in repeats if "chrom" in r and "chromStart" in r and "chromEnd" in r]
    
    # Apply limit if specified
    if limit and isinstance(limit, int) and limit > 0:
        valid_repeats = valid_repeats[:limit]
        print(f"Processing first {limit} out of {len(repeats)} repeats...")
    else:
        print(f"Processing {len(valid_repeats)} out of {len(repeats)} repeats with valid coordinates...")
    
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
        
        # Create a dictionary of exon IDs to exon objects from the overlap endpoint results
        # This gives us access to the phase information
        exon_phase_map = {}
        if "exons" in api_data and api_data["exons"]:
            for exon in api_data["exons"]:
                if "id" in exon:
                    exon_phase_map[exon["id"]] = exon
        
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
                gene_name = transcript.get("display_name", "").split('-')[0]
                strand = "+" if transcript.get("strand") == 1 else "-"
                
                # Find exons in this transcript
                exons = transcript.get("Exon", [])
                exon_count = len(exons)
                
                # Sort exons by genomic order
                if strand == "+":
                    # For + strand, exons are ordered from 5' to 3'
                    exons = sorted(exons, key=lambda e: e.get("start", 0))
                else:
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

                        # Add this adjustment for BED format
                        if exon_end == end:
                            overlap_end += 1  # Adjust for BED's exclusive end coordinate

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
                        
                        # Look up the exon in our phase map from the overlap endpoint
                        exon_id = exon.get("id", "")
                        overlap_exon = exon_phase_map.get(exon_id, {})
                        
                        # Get phase and end_phase from the overlap endpoint data
                        phase = overlap_exon.get("ensembl_phase", -1)
                        end_phase = overlap_exon.get("ensembl_end_phase", -1)
                        
                        frame_status = "non_coding"
                        if coding_status != "non_coding":
                            if phase == end_phase and phase != -1:
                                frame_status = "in_frame"  # Exon contains complete codons or maintains reading frame
                            elif phase == -1 or end_phase == -1:
                                frame_status = "non_coding"  # Non-coding exon
                            else:
                                frame_status = "out_of_frame"  # Exon contains partial codons
                        
                        exon_info = {
                            "exon_number": exon_number,
                            "exon_id": exon.get("id", ""),
                            "overlap_bp": overlap_length,
                            "position": position,
                            "overlap_percentage": round(overlap_percentage, 2),
                            "coding_status": coding_status,
                            "utr_status": utr_status,
                            "coding_percentage": coding_percentage,
                            "phase": phase,  # Store the direct value
                            "end_phase": end_phase,  # Store the direct value
                            "frame_status": frame_status
                        }
                        
                        containing_exons.append(exon_info)
                
                # Get transcript biotype
                biotype = transcript.get("biotype", "unknown")
                
                # Create versioned transcript ID
                versioned_transcript_id = transcript_id
                if "version" in transcript:
                    versioned_transcript_id = f"{transcript_id}.{transcript['version']}"
                
                transcript_info.append({
                    "transcript_id": transcript_id,  # Keep the unversioned ID for API queries
                    "versioned_transcript_id": versioned_transcript_id,  # Add this new field
                    "transcript_name": transcript.get("display_name", ""),
                    "is_canonical": is_canonical,
                    "biotype": biotype,
                    "location": location,
                    "exon_count": exon_count,
                    "containing_exons": containing_exons
                })
            except Exception as e:
                print(f"Error processing transcript {transcript.get('id', 'unknown')}: {e}")
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
    print(f"Updated repeat data saved to {output_file}")
    
    # Remove temp file if exists
    if os.path.exists(output_file + ".temp"):
        os.remove(output_file + ".temp")

if __name__ == "__main__":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory (project root)
    project_root = os.path.dirname(script_dir)
    
    # Set up argument parser for more user-friendly command line options
    import argparse
    parser = argparse.ArgumentParser(description="Process repeats and add Ensembl exon information.")
    parser.add_argument("--input", "-i", default=os.path.join(project_root, "data", "1000_gname_hg38_repeats.json"),
                        help="Input JSON file containing repeat data")
    parser.add_argument("--output", "-o", default=os.path.join(project_root, "data", "10_test_exons_hg38_repeats.json"),
                        help="Output JSON file to save results")
    parser.add_argument("--limit", "-l", type=int, default=None,
                        help="Limit processing to first N entries (e.g., 10, 100)")
    args = parser.parse_args()
    
    input_file = args.input
    output_file = args.output
    limit = args.limit
    
    # Display limit information
    if limit:
        print(f"Limiting processing to first {limit} entries in the input file")
    else:
        print("Processing all entries in the input file")
        
    try:
        # Start timing BEFORE processing
        start_time = time.time()
        
        # Run the processing with the specified limit
        process_repeat_data(input_file, output_file, limit=limit)
        
        # Calculate duration AFTER processing
        duration = time.time() - start_time
        
        # Log summary statistics
        logging.info("=== Processing Summary ===")
        logging.info(f"Total API requests: {api_stats['requests']}")
        logging.info(f"Rate limits encountered: {api_stats['rate_limits']}")
        logging.info(f"Errors: {api_stats['errors']}")
        logging.info(f"Total runtime: {duration:.2f} seconds")
        logging.info("========================")
    except Exception as e:
        logging.exception("Script failed with error:")