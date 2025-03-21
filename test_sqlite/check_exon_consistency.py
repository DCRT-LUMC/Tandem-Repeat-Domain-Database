import json
import os
import sqlite3
from collections import defaultdict

def connect_to_database(db_path):
    """Connect to the database"""
    if not os.path.exists(db_path):
        print(f"Database file not found: {db_path}")
        return None
    
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # Return rows as dictionaries
    return conn

def load_original_json(json_file):
    """Load and parse the original JSON data"""
    if not os.path.exists(json_file):
        print(f"JSON file not found: {json_file}")
        return None
    
    try:
        with open(json_file, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON file: {e}")
        return None

def check_exon_consistency(db_path, json_file=None):
    """Check exon consistency between database and original JSON"""
    conn = connect_to_database(db_path)
    if not conn:
        return
    
    cursor = conn.cursor()
    
    # First, check for exons with multiple positions in the database
    print("Checking for exons with multiple positions in the database...")
    
    cursor.execute("""
    SELECT 
        e.ensembl_exon_id, 
        COUNT(DISTINCT json_extract(e.ensembl_info, '$.position')) as position_count,
        GROUP_CONCAT(DISTINCT json_extract(e.ensembl_info, '$.position')) as positions,
        GROUP_CONCAT(DISTINCT json_extract(e.ensembl_info, '$.transcript_id')) as transcripts
    FROM exons e
    WHERE e.ensembl_info IS NOT NULL
    GROUP BY e.ensembl_exon_id
    HAVING position_count > 1
    ORDER BY position_count DESC
    LIMIT 20
    """)
    
    results = cursor.fetchall()
    
    if not results:
        print("No exons found with multiple positions in the database.")
    else:
        print(f"\nFound {len(results)} exons with multiple positions in the database:")
        for row in results:
            print(f"\nExon ID: {row['ensembl_exon_id']}")
            print(f"Position Count: {row['position_count']}")
            print(f"Positions: {row['positions']}")
            print(f"Transcripts: {row['transcripts']}")
            
            # Get detailed info for this exon
            cursor.execute("""
            SELECT 
                e.exon_id,
                e.ensembl_exon_id,
                e.ensembl_info,
                r.repeat_id,
                g.gene_name
            FROM exons e
            JOIN repeats r ON e.repeat_id = r.repeat_id
            JOIN proteins p ON r.protein_id = p.protein_id
            JOIN genes g ON p.gene_id = g.gene_id
            WHERE e.ensembl_exon_id = ?
            """, (row['ensembl_exon_id'],))
            
            exon_details = cursor.fetchall()
            if exon_details:
                print(f"Gene: {exon_details[0]['gene_name']}")
                print("Detailed position data:")
                for detail in exon_details:
                    try:
                        info = json.loads(detail['ensembl_info'])
                        transcript_id = info.get('transcript_id', 'Unknown')
                        transcript_name = info.get('transcript_name', 'Unknown')
                        position = info.get('position', 'Unknown')
                        exon_number = info.get('exon_number', 'Unknown')
                        print(f"  Transcript {transcript_id} ({transcript_name}):")
                        print(f"    Position: {position}")
                        print(f"    Exon Number: {exon_number}")
                    except json.JSONDecodeError:
                        print(f"  Error parsing JSON: {detail['ensembl_info']}")
    
    # If JSON file is provided, check consistency with original data
    if json_file and os.path.exists(json_file):
        original_data = load_original_json(json_file)
        if not original_data:
            return
        
        print("\nChecking consistency with original JSON data...")
        
        # Create a mapping of exon positions from original data
        original_positions = {}
        for item in original_data:
            if 'ensembl_exon_info' in item and 'transcripts' in item['ensembl_exon_info']:
                for transcript in item['ensembl_exon_info']['transcripts']:
                    transcript_id = transcript.get('transcript_id')
                    if not transcript_id:
                        continue
                    
                    if 'containing_exons' in transcript:
                        for exon in transcript['containing_exons']:
                            exon_id = exon.get('exon_id')
                            if not exon_id:
                                continue
                            
                            if exon_id not in original_positions:
                                original_positions[exon_id] = {}
                            
                            original_positions[exon_id][transcript_id] = {
                                'position': exon.get('position'),
                                'exon_number': exon.get('exon_number'),
                                'transcript_name': transcript.get('transcript_name')
                            }
        
        # Check a sample of exons against the original data
        cursor.execute("""
        SELECT 
            e.ensembl_exon_id,
            e.ensembl_info,
            g.gene_name
        FROM exons e
        JOIN repeats r ON e.repeat_id = r.repeat_id
        JOIN proteins p ON r.protein_id = p.protein_id
        JOIN genes g ON p.gene_id = g.gene_id
        WHERE e.ensembl_info IS NOT NULL
        LIMIT 100
        """)
        
        sample_exons = cursor.fetchall()
        inconsistencies = []
        
        print(f"Checking {len(sample_exons)} sample exons against original JSON...")
        
        for exon in sample_exons:
            exon_id = exon['ensembl_exon_id']
            if exon_id in original_positions:
                original = original_positions[exon_id]
                
                try:
                    db_info = json.loads(exon['ensembl_info'])
                    db_transcript_id = db_info.get('transcript_id')
                    
                    if db_transcript_id and db_transcript_id in original:
                        db_position = db_info.get('position')
                        original_position = original[db_transcript_id].get('position')
                        
                        if db_position and original_position and db_position != original_position:
                            inconsistencies.append({
                                'exon_id': exon_id,
                                'gene': exon['gene_name'],
                                'transcript_id': db_transcript_id,
                                'db_position': db_position,
                                'original_position': original_position
                            })
                except json.JSONDecodeError:
                    pass
        
        if not inconsistencies:
            print("No inconsistencies found in the sample!")
        else:
            print(f"\nFound {len(inconsistencies)} inconsistencies between database and original JSON:")
            for inc in inconsistencies:
                print(f"\nExon ID: {inc['exon_id']} in gene {inc['gene']}")
                print(f"Transcript: {inc['transcript_id']}")
                print(f"Database Position: {inc['db_position']}")
                print(f"Original Position: {inc['original_position']}")
    
    conn.close()

if __name__ == "__main__":
    db_path = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\repeats.db'
    json_file = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\1000_test_exons_hg38_repeats.json'
    
    check_exon_consistency(db_path, json_file)
