import json
import requests
import time
from tqdm import tqdm
import os
import sys
import logging
import datetime
import sqlite3

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
            
            # Get detailed transcript information with all exons and version information
            transcript_detail = client.perform_rest_action(
                endpoint=f"/lookup/id/{transcript_id}",
                hdrs=headers,
                params={'expand': 1, 'format': 'full'}  # Add format=full parameter
            )
            
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

def initialize_database(db_path):
    """
    Initialize the SQLite database and create necessary tables if they don't exist
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Create tables
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS repeats (
        repeat_id INTEGER PRIMARY KEY,
        chrom TEXT,
        start INTEGER,
        end INTEGER,
        sequence TEXT,
        motif TEXT,
        num_copies REAL,
        strand TEXT,
        location_summary TEXT
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS transcripts (
        transcript_id TEXT PRIMARY KEY,
        versioned_id TEXT,
        transcript_name TEXT,
        gene_name TEXT,
        biotype TEXT,
        is_canonical INTEGER
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS repeat_transcript (
        repeat_id INTEGER,
        transcript_id TEXT,
        location TEXT,
        PRIMARY KEY (repeat_id, transcript_id),
        FOREIGN KEY (repeat_id) REFERENCES repeats (repeat_id),
        FOREIGN KEY (transcript_id) REFERENCES transcripts (transcript_id)
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS exons (
        exon_id TEXT PRIMARY KEY,
        transcript_id TEXT,
        exon_number INTEGER,
        start INTEGER,
        end INTEGER,
        phase INTEGER,
        end_phase INTEGER,
        coding_status TEXT,
        FOREIGN KEY (transcript_id) REFERENCES transcripts (transcript_id)
    )
    ''')
    
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS repeat_exon (
        repeat_id INTEGER,
        exon_id TEXT,
        overlap_bp INTEGER,
        overlap_percentage REAL,
        position TEXT,
        utr_status TEXT,
        coding_percentage REAL,
        frame_status TEXT,
        PRIMARY KEY (repeat_id, exon_id),
        FOREIGN KEY (repeat_id) REFERENCES repeats (repeat_id),
        FOREIGN KEY (exon_id) REFERENCES exons (exon_id)
    )
    ''')
    
    conn.commit()
    return conn

def process_repeat_data_to_db(repeat_data_file, db_path, limit=None):
    """
    Process the repeat data JSON and add exon information using Ensembl API.
    Insert data directly into SQLite database.
    
    Parameters:
        repeat_data_file: Input JSON file with repeat data
        db_path: Path to SQLite database
        limit: Optional. If set, process only this many entries
    """
    
    # Initialize database
    conn = initialize_database(db_path)
    cursor = conn.cursor()
    
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
    
    # Start a transaction - commit after every 50 repeats
    conn.execute("BEGIN TRANSACTION")
    commit_counter = 0
    
    # Process each repeat
    for repeat_idx, repeat in enumerate(tqdm(valid_repeats)):
        # Commit transaction every 50 repeats
        commit_counter += 1
        if commit_counter >= 50:
            conn.commit()
            conn.execute("BEGIN TRANSACTION")
            commit_counter = 0
            logging.info(f"Committed batch at repeat {repeat_idx}")
        
        chrom = repeat["chrom"]
        start = int(repeat["chromStart"])
        end = int(repeat["chromEnd"])
        sequence = repeat.get("sequence", "")
        motif = repeat.get("repeatUnit", "")
        num_copies = repeat.get("copyNumber", 0)
        strand = repeat.get("strand", "")

        # Insert repeat into the database
        cursor.execute(
            """INSERT INTO repeats (chrom, start, end, sequence, motif, num_copies, strand, location_summary) 
               VALUES (?, ?, ?, ?, ?, ?, ?, 'unknown')
               RETURNING repeat_id""", 
            (chrom, start, end, sequence, motif, num_copies, strand)
        )
        repeat_id = cursor.fetchone()[0]
        
        # Convert strand to Ensembl format for comparison
        expected_ensembl_strand = 1 if strand == "+" else -1 if strand == "-" else None
        
        # Get transcript and exon information from Ensembl
        api_data = get_ensembl_info(chrom, start, end)
        
        if not api_data or not api_data["transcripts"]:
            # Update repeat with location_summary
            cursor.execute(
                "UPDATE repeats SET location_summary = 'unknown' WHERE repeat_id = ?",
                (repeat_id,)
            )
            continue
        
        all_transcripts = []
        transcript_details = api_data["transcript_details"]
        
        # Use transcript details for comprehensive information
        for transcript_id, transcript in transcript_details.items():
            all_transcripts.append(transcript)
        
        locations = set()
        has_canonical = False
        
        for transcript in all_transcripts:
            try:
                # Skip transcripts with different strand if repeat strand is specified
                if expected_ensembl_strand is not None and transcript.get("strand") != expected_ensembl_strand:
                    continue
                
                # Get basic transcript info
                transcript_id = transcript["id"]
                gene_name = transcript.get("display_name", "").split('-')[0]
                transcript_name = transcript.get("display_name", "")
                biotype = transcript.get("biotype", "unknown")
                
                # Check if this is likely the canonical transcript
                is_canonical = is_canonical_transcript(transcript, all_transcripts)
                if is_canonical:
                    has_canonical = True
                
                # Create versioned transcript ID
                versioned_transcript_id = transcript_id
                if "version" in transcript:
                    versioned_transcript_id = f"{transcript_id}.{transcript['version']}"
                
                # Insert transcript
                cursor.execute(
                    """INSERT OR IGNORE INTO transcripts 
                       (transcript_id, versioned_id, transcript_name, gene_name, biotype, is_canonical) 
                       VALUES (?, ?, ?, ?, ?, ?)""",
                    (transcript_id, versioned_transcript_id, transcript_name, gene_name, biotype, 1 if is_canonical else 0)
                )
                
                # Classify location (exonic, intronic, outside)
                location = classify_repeat_location(repeat, transcript)
                locations.add(location)
                
                # Insert repeat-transcript relationship
                cursor.execute(
                    """INSERT INTO repeat_transcript (repeat_id, transcript_id, location) 
                       VALUES (?, ?, ?)""",
                    (repeat_id, transcript_id, location)
                )
                
                # Find exons in this transcript
                exons = transcript.get("Exon", [])
                strand_value = transcript.get("strand")
                strand = "+" if strand_value == 1 else "-"
                
                # Sort exons by genomic order
                if strand == "+":
                    # For + strand, exons are ordered from 5' to 3'
                    exons = sorted(exons, key=lambda e: e.get("start", 0))
                else:
                    # For - strand, exons are ordered from 3' to 5'
                    exons = sorted(exons, key=lambda e: e.get("start", 0), reverse=True)
                
                exon_count = len(exons)
                
                for i, exon in enumerate(exons):
                    exon_id = exon.get("id", "")
                    exon_start = int(exon.get("start", 0))
                    exon_end = int(exon.get("end", 0))
                    exon_number = i + 1  # 1-based exon numbering
                    phase = exon.get("phase", -1)
                    end_phase = exon.get("end_phase", -1)
                    
                    # Insert exon
                    coding_status, _, _ = get_coding_status(exon, transcript, start, end)
                    cursor.execute(
                        """INSERT OR IGNORE INTO exons 
                           (exon_id, transcript_id, exon_number, start, end, phase, end_phase, coding_status) 
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                        (exon_id, transcript_id, exon_number, exon_start, exon_end, phase, end_phase, coding_status)
                    )
                    
                    # Check if the repeat overlaps this exon
                    if max(start, exon_start) < min(end, exon_end):  # Overlap
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
                        _, utr_status, coding_percentage = get_coding_status(exon, transcript, start, end)
                        
                        # Determine frame status
                        frame_status = "non_coding"
                        if coding_status != "non_coding":
                            if phase == end_phase and phase != -1 and end_phase != -1:
                                frame_status = "in_frame"
                            elif phase == -1 or end_phase == -1:
                                frame_status = "non_coding"
                            else:
                                frame_status = "out_of_frame"
                                
                        # Insert repeat-exon relationship
                        cursor.execute(
                            """INSERT INTO repeat_exon 
                               (repeat_id, exon_id, overlap_bp, overlap_percentage, position, utr_status, 
                                coding_percentage, frame_status) 
                               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                            (repeat_id, exon_id, overlap_length, round(overlap_percentage, 2), 
                             position, utr_status, coding_percentage, frame_status)
                        )
                
            except Exception as e:
                logging.error(f"Error processing transcript {transcript.get('id', 'unknown')}: {e}")
                continue
        
        # Summarize location (prioritize exonic > intronic > outside/unknown)
        location_summary = "unknown"
        if "exonic" in locations:
            location_summary = "exonic"
        elif "intronic" in locations:
            location_summary = "intronic"
        elif "outside" in locations:
            location_summary = "intergenic"
            
        # Update repeat with location_summary
        cursor.execute(
            "UPDATE repeats SET location_summary = ? WHERE repeat_id = ?",
            (location_summary, repeat_id)
        )
    
    # Commit final changes
    conn.commit()
    
    # Add indexes for better query performance
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_repeat_location ON repeats (location_summary)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_transcript_canonical ON transcripts (is_canonical)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_exon_transcript ON exons (transcript_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_repeat_exon_overlap ON repeat_exon (overlap_percentage)")
    
    conn.commit()
    conn.close()
    print(f"Database populated at {db_path}")

if __name__ == "__main__":
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Get the parent directory (project root)
    project_root = os.path.dirname(script_dir)
    
    # Default paths relative to project root
    input_file = os.path.join(project_root, "data", "1000_gname_hg38_repeats.json")
    db_file = os.path.join(project_root, "data", "repeats_database.db")
    
    # Allow command-line arguments to override default file paths
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        db_file = sys.argv[2]
        
    try:
        # Start timing BEFORE processing
        start_time = time.time()
        
        # Run the processing directly to database
        process_repeat_data_to_db(input_file, db_file, limit=None)
        
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