import sqlite3
import json
import os

def connect_to_database(db_path):
    """Connect to the database"""
    if not os.path.exists(db_path):
        print(f"Database file not found: {db_path}")
        return None
    
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row  # Return rows as dictionaries
    return conn

def display_database_stats(conn):
    """Display basic statistics about the database"""
    cursor = conn.cursor()
    
    # Count records in each table
    tables = ['genes', 'gene_aliases', 'proteins', 'repeats', 'transcripts', 'repeat_transcripts', 'exons']
    print("Database Statistics:")
    print("===================")
    
    for table in tables:
        cursor.execute(f"SELECT COUNT(*) as count FROM {table}")
        count = cursor.fetchone()['count']
        print(f"{table.capitalize()}: {count} records")
    
    # Count unique repeat types
    cursor.execute("SELECT repeat_type, COUNT(*) as count FROM repeats GROUP BY repeat_type ORDER BY count DESC")
    repeat_types = cursor.fetchall()
    
    print("\nRepeat Types:")
    for rt in repeat_types:
        print(f"  {rt['repeat_type']}: {rt['count']} repeats")

def find_exon_skipping_subjects(conn, min_repeats=5, min_overlap_percentage=70):
    """Find genes with exon skipping potential
    
    This finds genes that have:
    1. More than 'min_repeats' of the same repeat type
    2. At least one exon with 'min_overlap_percentage' or more overlap
    3. Exons that are in frame
    4. Repeats with block_count = 1 (single-block repeats)eats)
    """
    cursor = conn.cursor()
    
    query = """
    -- Step 1: Find genes with multiple repeats of the same type
    WITH gene_repeat_counts AS (
        SELECT g.gene_id, g.gene_name, r.repeat_type, COUNT(DISTINCT r.repeat_id) as repeat_count
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        -- WHERE r.block_count = 1  -- Only include repeats with block_count = 1
        GROUP BY g.gene_id, r.repeat_type
        HAVING COUNT(DISTINCT r.repeat_id) > ?
    )
    
    -- Step 2: Find exons with significant overlap and in-frame status
    SELECT 
        grc.gene_name,
        grc.repeat_type,
        grc.repeat_count,
        e.ensembl_exon_id,
        e.frame_status,
        r.repeat_id,
        r.position,
        te.overlap_percentage,
        te.overlap_bp,
        te.exon_position_in_transcript,
        te.transcript_id,
        t.transcript_name,
        te.exon_number,
        r.block_count,
        r.repeat_length
    FROM gene_repeat_counts grc
    JOIN genes g ON grc.gene_id = g.gene_id
    JOIN proteins p ON g.gene_id = p.gene_id
    JOIN repeats r ON p.protein_id = r.protein_id AND r.repeat_type = grc.repeat_type
    JOIN repeat_exons re ON r.repeat_id = re.repeat_id
    JOIN exons e ON re.exon_id = e.exon_id
    JOIN transcript_exons te ON e.exon_id = te.exon_id
    JOIN transcripts t ON te.transcript_id = t.transcript_id
    WHERE 
        e.frame_status = 'in_frame'
        AND te.overlap_percentage >= ?
        AND r.block_count = 1  -- Only include repeats with block_count = 1ude repeats with block_count = 1
    ORDER BY grc.gene_name, grc.repeat_count DESC, te.overlap_percentage DESC
    """
    
    cursor.execute(query, (min_repeats, min_overlap_percentage))
    results = cursor.fetchall()
    
    # Organize results by gene
    gene_results = {}
    for row in results:
        gene_name = row['gene_name']
        
        if gene_name not in gene_results:
            gene_results[gene_name] = {
                'gene_name': gene_name,
                'repeat_type': row['repeat_type'],
                'repeat_count': row['repeat_count'],
                'exons': []
            }
        
        gene_results[gene_name]['exons'].append({
            'ensembl_exon_id': row['ensembl_exon_id'],
            'frame_status': row['frame_status'],
            'repeat_id': row['repeat_id'],
            'position': row['position'],
            'overlap_percentage': row['overlap_percentage'],
            'overlap_bp': row['overlap_bp'],
            'exon_position': row['exon_position_in_transcript'],
            'transcript_id': row['transcript_id'],
            'transcript_name': row['transcript_name'],
            'exon_number': row['exon_number'],
            'block_count': row['block_count'],
            'repeat_length': row['repeat_length']
        })
    
    # Display results
    print(f"\nExon Skipping Subjects (>{min_repeats} repeats, >{min_overlap_percentage}% overlap, block_count = 1):")
    print("=====================================================================")
    
    if not gene_results:
        print("No matching genes found.")
        return
    
    for gene_name, data in gene_results.items():
        print(f"\nGene: {gene_name}")
        print(f"Repeat Type: {data['repeat_type']}")
        print(f"Repeat Count: {data['repeat_count']}")
        print(f"In-frame Exons with high overlap: {len(data['exons'])}")
        
        print("\nExon Details:")
        print("------------")
        for idx, exon in enumerate(data['exons'][:5], 1):  # Show first 5 exons
            print(f"{idx}. ID: {exon['ensembl_exon_id']}")
            # Include transcript context in the position display
            transcript_info = f" (in transcript {exon['transcript_id']})"
            print(f"   Position: {exon['exon_position']}{transcript_info}")
            print(f"   Exon Number: {exon['exon_number']} in {exon['transcript_name']}")
            print(f"   Overlap: {exon['overlap_percentage']}% ({exon['overlap_bp']} bp)")
            print(f"   Repeat Position: {exon['position']}")
            print(f"   Block Count: {exon['block_count']} (matches repeat length: {exon['repeat_length']})")
        
        if len(data['exons']) > 5:
            print(f"   ... and {len(data['exons']) - 5} more exons")

def find_specific_gene_repeats(conn, gene_name):
    """Find repeats for a specific gene"""
    cursor = conn.cursor()
    
    query = """
    SELECT 
        g.gene_name,
        r.repeat_type,
        r.repeat_id,
        r.position,
        r.chrom,
        r.chrom_start,
        r.chrom_end,
        r.strand,
        r.repeat_length,
        p.protein_id,
        p.status
    FROM genes g
    JOIN proteins p ON g.gene_id = p.gene_id
    JOIN repeats r ON p.protein_id = r.protein_id
    WHERE g.gene_name = ?
    ORDER BY r.repeat_type, r.chrom_start
    """
    
    cursor.execute(query, (gene_name,))
    results = cursor.fetchall()
    
    print(f"\nRepeats for gene '{gene_name}':")
    print("=============================")
    
    if not results:
        print(f"No repeats found for gene '{gene_name}'.")
        return
    
    # Group by repeat type
    repeats_by_type = {}
    for row in results:
        repeat_type = row['repeat_type']
        
        if repeat_type not in repeats_by_type:
            repeats_by_type[repeat_type] = []
        
        repeats_by_type[repeat_type].append(row)
    
    for repeat_type, repeats in repeats_by_type.items():
        print(f"\nRepeat Type: {repeat_type} ({len(repeats)} repeats)")
        print("----------------------------")
        
        for idx, repeat in enumerate(repeats, 1):
            print(f"{idx}. {repeat['position']}")
            print(f"   Location: {repeat['chrom']}:{repeat['chrom_start']}-{repeat['chrom_end']} ({repeat['strand']})")
            print(f"   Length: {repeat['repeat_length']}")
            print(f"   Protein: {repeat['protein_id']} ({repeat['status']})")

def debug_exon_info(conn, exon_id=None):
    """Debug function to examine specific exon information"""
    cursor = conn.cursor()
    
    if exon_id:
        # Query for a specific exon
        query = """
        SELECT 
            e.exon_id,
            e.ensembl_exon_id,
            e.frame_status,
            e.ensembl_info,
            r.repeat_id,
            g.gene_name,
            r.repeat_type,
            r.position AS repeat_position,
            r.chrom_start,
            r.chrom_end,
            t.transcript_id,
            t.ensembl_transcript_id
        FROM exons e
        JOIN repeats r ON e.repeat_id = r.repeat_id
        JOIN proteins p ON r.protein_id = p.protein_id
        JOIN genes g ON p.gene_id = g.gene_id
        LEFT JOIN repeat_transcripts rt ON r.repeat_id = rt.repeat_id
        LEFT JOIN transcripts t ON rt.transcript_id = t.transcript_id
        WHERE e.ensembl_exon_id = ?
        ORDER BY r.repeat_id
        """
        cursor.execute(query, (exon_id,))
    else:
        # Get a sample of exons
        query = """
        SELECT 
            e.exon_id,
            e.ensembl_exon_id,
            e.frame_status,
            e.ensembl_info,
            r.repeat_id,
            g.gene_name,
            r.repeat_type,
            r.position AS repeat_position,
            r.chrom_start,
            r.chrom_end,
            t.transcript_id,
            t.ensembl_transcript_id
        FROM exons e
        JOIN repeats r ON e.repeat_id = r.repeat_id
        JOIN proteins p ON r.protein_id = p.protein_id
        JOIN genes g ON p.gene_id = g.gene_id
        LEFT JOIN repeat_transcripts rt ON r.repeat_id = rt.repeat_id
        LEFT JOIN transcripts t ON rt.transcript_id = t.transcript_id
        WHERE e.ensembl_info IS NOT NULL
        LIMIT 5
        """
        cursor.execute(query)
    
    results = cursor.fetchall()
    
    print("\nExon Debug Information:")
    print("======================")
    
    if not results:
        print(f"No exon found" + (f" with ID {exon_id}" if exon_id else ""))
        return
    
    # Get original JSON data for verification
    json_file = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\1000_test_exons_hg38_repeats.json'
    original_data_by_transcript = {}
    
    if os.path.exists(json_file):
        try:
            with open(json_file, 'r') as f:
                original_data = json.load(f)
            
            # Extract original positions for the specific exon ID, organized by transcript
            for item in original_data:
                if 'ensembl_exon_info' in item and 'transcripts' in item['ensembl_exon_info']:
                    for transcript in item['ensembl_exon_info']['transcripts']:
                        transcript_id = transcript.get('transcript_id')
                        if not transcript_id:
                            continue
                            
                        if 'containing_exons' in transcript:
                            for exon in transcript['containing_exons']:
                                if exon_id and exon.get('exon_id') == exon_id:
                                    if transcript_id not in original_data_by_transcript:
                                        original_data_by_transcript[transcript_id] = {}
                                    
                                    original_data_by_transcript[transcript_id] = {
                                        'exon_number': exon.get('exon_number'),
                                        'position': exon.get('position'),
                                        'transcript_name': transcript.get('transcript_name'),
                                        'is_canonical': transcript.get('is_canonical')
                                    }
        except Exception as e:
            print(f"Error reading original JSON: {e}")
    
    # Print original data by transcript
    if exon_id and original_data_by_transcript:
        print(f"\nOriginal data for exon {exon_id} by transcript:")
        for transcript_id, data in original_data_by_transcript.items():
            canonical = " (CANONICAL)" if data.get('is_canonical') else ""
            print(f"  Transcript {transcript_id} - {data.get('transcript_name')}{canonical}:")
            print(f"    Exon Number: {data.get('exon_number')}")
            print(f"    Position: {data.get('position')}")
    
    for row in results:
        print(f"\nExon ID: {row['exon_id']}")
        print(f"Ensembl Exon ID: {row['ensembl_exon_id']}")
        print(f"Gene: {row['gene_name']}")
        print(f"Repeat Type: {row['repeat_type']}")
        print(f"Frame Status: {row['frame_status']}")
        print(f"Repeat Position: {row['repeat_position']}")
        print(f"Genomic Range: {row['chrom_start']}-{row['chrom_end']}")
        
        # Parse the ensembl_info JSON
        if row['ensembl_info']:
            try:
                ensembl_info = json.loads(row['ensembl_info'])
                print("\nEnsembl Info from Database:")
                
                # First, show the primary transcript info
                primary_transcript_id = ensembl_info.get('transcript_id')
                primary_transcript_name = ensembl_info.get('transcript_name')
                
                if primary_transcript_id:
                    print(f"  Primary Transcript: {primary_transcript_id} ({primary_transcript_name})")
                    
                for key, value in ensembl_info.items():
                    if key not in ['transcript_id', 'transcript_name', 'other_transcripts']:
                        print(f"  {key}: {value}")
                
                # Then, show other transcripts' info if available
                if 'other_transcripts' in ensembl_info and ensembl_info['other_transcripts']:
                    print("\n  Other Transcripts Data:")
                    for t_id, t_data in ensembl_info['other_transcripts'].items():
                        print(f"    {t_id} ({t_data.get('transcript_name', 'Unknown')}):")
                        print(f"      Exon Number: {t_data.get('exon_number')}")
                        print(f"      Position: {t_data.get('position')}")
                
                # Check for position mismatch with the original JSON
                if primary_transcript_id in original_data_by_transcript:
                    original = original_data_by_transcript[primary_transcript_id]
                    stored_position = ensembl_info.get('position')
                    original_position = original.get('position')
                    
                    if stored_position and original_position and stored_position != original_position:
                        print(f"\n  POSITION MISMATCH for transcript {primary_transcript_id}:")
                        print(f"    Database has: {stored_position}")
                        print(f"    Original JSON has: {original_position}")
                
            except json.JSONDecodeError:
                print(f"Error parsing ensembl_info: {row['ensembl_info']}")

def main():
    db_path = 'test_sqlite/repeats.db'
    
    conn = connect_to_database(db_path)
    if not conn:
        return
    
    try:
        # Display database statistics
        display_database_stats(conn)
        
        # Debug a specific exon by its Ensembl ID
        # debug_exon_info(conn, "ENSE00000769655")  # Check this specific exon
        
        # Also check a few other exons to see if the issue is widespread
        # print("\nChecking other exons for possible indexing issues:")
        # debug_exon_info(conn, "ENSE00001368267")  # Another exon to check
        
        # Find exon skipping subjects
        find_exon_skipping_subjects(conn, min_repeats=5, min_overlap_percentage=70)
        
    finally:
        conn.close()

if __name__ == "__main__":
    main()
