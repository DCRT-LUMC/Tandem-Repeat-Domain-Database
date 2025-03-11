import json
import requests
import time
import pandas as pd
from tqdm import tqdm
import os
from collections import defaultdict
import sys

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
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            # Use requests instead of urllib for consistency with the rest of your code
            response = requests.get(url, headers=hdrs, params=params, timeout=15)
            
            if response.status_code == 200:
                data = response.json()
                self.req_count += 1
                
            # Check if we are being rate limited by the server
            elif response.status_code == 429:
                if 'Retry-After' in response.headers:
                    retry = response.headers['Retry-After']
                    time.sleep(float(retry))
                    return self.perform_rest_action(endpoint, hdrs, params)
            else:
                print(f"Request failed for {endpoint}: Status code: {response.status_code}")
                
        except Exception as e:
            print(f"Error performing request to {endpoint}: {e}")
            
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
    
    # Remove temp file if exists
    if os.path.exists(output_file + ".temp"):
        os.remove(output_file + ".temp")
    
    print(f"Updated repeat data saved to {output_file}")

if __name__ == "__main__":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory (project root)
    project_root = os.path.dirname(script_dir)
    
    # Default paths relative to project root
    input_file = os.path.join(project_root, "data", "gname_hg38_repeats_100.json")
    output_file = os.path.join(project_root, "data", "ensembl_exons_hg38_repeats_10.json")
    
    # Allow command-line arguments to override default file paths
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
        
    # Process only the first 10 repeats
    process_repeat_data(input_file, output_file, limit=10)