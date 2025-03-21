import json
import sqlite3
import os
import sys
from collections import defaultdict

def create_tables(conn, schema_file):
    """Create database tables from schema file"""
    with open(schema_file, 'r') as f:
        schema_sql = f.read()
        conn.executescript(schema_sql)
    print("Database tables created successfully")

def insert_or_get_gene(cursor, gene_name, chromosome=None):
    """Insert gene if not exists and return gene_id"""
    cursor.execute(
        "SELECT gene_id FROM genes WHERE gene_name = ?", 
        (gene_name,)
    )
    result = cursor.fetchone()
    
    if result:
        return result[0]
    
    cursor.execute(
        "INSERT INTO genes (gene_name, chromosome) VALUES (?, ?)",
        (gene_name, chromosome)
    )
    return cursor.lastrowid

def insert_gene_aliases(cursor, gene_id, aliases):
    """Insert gene aliases"""
    if not aliases:
        return
        
    if isinstance(aliases, list):
        for alias in aliases:
            if alias:
                cursor.execute(
                    "INSERT OR IGNORE INTO gene_aliases (gene_id, alias_name) VALUES (?, ?)",
                    (gene_id, alias)
                )
    elif isinstance(aliases, str) and aliases:
        cursor.execute(
            "INSERT OR IGNORE INTO gene_aliases (gene_id, alias_name) VALUES (?, ?)",
            (gene_id, aliases)
        )

def insert_or_get_protein(cursor, uniprot_id, gene_id, status=None):
    """Insert protein if not exists and return protein_id"""
    cursor.execute(
        "SELECT protein_id FROM proteins WHERE protein_id = ?", 
        (uniprot_id,)
    )
    result = cursor.fetchone()
    
    if result:
        return result[0]
    
    cursor.execute(
        "INSERT INTO proteins (protein_id, gene_id, status) VALUES (?, ?, ?)",
        (uniprot_id, gene_id, status)
    )
    return uniprot_id

def insert_repeat(cursor, protein_id, repeat_data):
    """Insert repeat and return repeat_id"""
    cursor.execute(
        """
        INSERT INTO repeats 
        (protein_id, repeat_type, chrom, chrom_start, chrom_end, 
        strand, position, repeat_length, reserved) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            protein_id, 
            repeat_data.get('repeatType'), 
            repeat_data.get('chrom'), 
            repeat_data.get('chromStart'), 
            repeat_data.get('chromEnd'),
            repeat_data.get('strand'), 
            repeat_data.get('position'), 
            repeat_data.get('repeatLength'),
            json.dumps(repeat_data.get('reserved', []))
        )
    )
    repeat_id = cursor.lastrowid
    
    # Insert into genomic_coordinates table
    cursor.execute(
        """
        INSERT INTO genomic_coordinates 
        (repeat_id, chrom_start, chrom_end, strand) 
        VALUES (?, ?, ?, ?)
        """,
        (
            repeat_id, 
            repeat_data.get('chromStart'), 
            repeat_data.get('chromEnd'), 
            repeat_data.get('strand')
        )
    )
    
    # Extract start and end positions from position field (e.g., "amino acids 343-389 on protein Q6TDP4")
    position = repeat_data.get('position', '')
    start_pos = None
    end_pos = None
    
    if position and 'amino acids' in position:
        try:
            pos_part = position.split('amino acids ')[1].split(' on')[0]
            if '-' in pos_part:
                start_pos, end_pos = map(int, pos_part.split('-'))
        except (IndexError, ValueError):
            pass
    
    # Update start_pos and end_pos
    if start_pos is not None and end_pos is not None:
        cursor.execute(
            "UPDATE repeats SET start_pos = ?, end_pos = ? WHERE repeat_id = ?",
            (start_pos, end_pos, repeat_id)
        )
    
    # Insert exons
    block_count = repeat_data.get('blockCount', 0)
    block_sizes = repeat_data.get('blockSizes', [])
    block_starts = repeat_data.get('chromStarts', [])
    
    if block_count and block_sizes and block_starts:
        if isinstance(block_sizes, list):
            block_sizes_str = ','.join(map(str, block_sizes))
        else:
            block_sizes_str = str(block_sizes)
            
        if isinstance(block_starts, list):
            block_starts_str = ','.join(map(str, block_starts))
        else:
            block_starts_str = str(block_starts)
            
        cursor.execute(
            """
            INSERT INTO exons 
            (repeat_id, block_count, block_sizes, block_starts) 
            VALUES (?, ?, ?, ?)
            """,
            (repeat_id, block_count, block_sizes_str, block_starts_str)
        )
    
    return repeat_id

def insert_or_get_transcript(cursor, transcript_data, gene_id):
    """Insert transcript if not exists and return transcript_id"""
    transcript_id = transcript_data.get('transcript_id')
    
    cursor.execute(
        "SELECT transcript_id FROM transcripts WHERE transcript_id = ?", 
        (transcript_id,)
    )
    result = cursor.fetchone()
    
    if result:
        return result[0]
    
    cursor.execute(
        """
        INSERT INTO transcripts 
        (transcript_id, gene_id, ensembl_transcript_id) 
        VALUES (?, ?, ?)
        """,
        (
            transcript_id, 
            gene_id, 
            transcript_data.get('versioned_transcript_id')
        )
    )
    return transcript_id

def process_ensembl_info(cursor, repeat_id, ensembl_info, gene_id):
    """Process ensembl exon information and insert related records"""
    if not ensembl_info or not isinstance(ensembl_info, dict):
        return
    
    transcripts = ensembl_info.get('transcripts', [])
    if not transcripts:
        return
    
    for transcript_data in transcripts:
        if not isinstance(transcript_data, dict):
            continue
            
        transcript_id = transcript_data.get('transcript_id')
        if not transcript_id:
            continue
            
        # Insert or get transcript
        transcript_id = insert_or_get_transcript(cursor, transcript_data, gene_id)
        
        # Link repeat to transcript
        cursor.execute(
            """
            INSERT OR IGNORE INTO repeat_transcripts 
            (repeat_id, transcript_id) 
            VALUES (?, ?)
            """,
            (repeat_id, transcript_id)
        )
        
        # Process exons
        containing_exons = transcript_data.get('containing_exons', [])
        if not containing_exons:
            continue
            
        for exon_data in containing_exons:
            if not isinstance(exon_data, dict):
                continue
                
            # Process exon
            ensembl_exon_id = exon_data.get('exon_id')
            if not ensembl_exon_id:
                continue
            
            # Store exon number and position with transcript context
            exon_number = exon_data.get('exon_number')
            exon_position = exon_data.get('position')
            
            # Add exon information to existing exon record or create a new one
            cursor.execute(
                """
                SELECT exon_id, ensembl_info FROM exons 
                WHERE repeat_id = ? AND ensembl_exon_id = ?
                """, 
                (repeat_id, ensembl_exon_id)
            )
            exon_result = cursor.fetchone()
            
            phase = exon_data.get('phase')
            end_phase = exon_data.get('end_phase')
            frame_status = exon_data.get('frame_status')
            
            # Create JSON structure with exon information including transcript context
            ensembl_info_dict = {
                'overlap_bp': exon_data.get('overlap_bp'),
                'position': exon_position,
                'overlap_percentage': exon_data.get('overlap_percentage'),
                'coding_status': exon_data.get('coding_status'),
                'utr_status': exon_data.get('utr_status'),
                'coding_percentage': exon_data.get('coding_percentage'),
                'exon_number': exon_number,
                'transcript_id': transcript_id,  # Add transcript context
                'transcript_name': transcript_data.get('transcript_name')  # Add transcript name
            }
            
            # Convert to JSON string
            ensembl_info_str = json.dumps(ensembl_info_dict)
            
            if exon_result:
                # Exon record exists - check if we need to update or create a new transcript-specific record
                existing_exon_id = exon_result[0]
                existing_info = exon_result[1]
                
                # Process existing info
                if existing_info:
                    try:
                        existing_info_dict = json.loads(existing_info)
                        existing_transcript_id = existing_info_dict.get('transcript_id')
                        
                        # If this is info for a different transcript, store both versions
                        if existing_transcript_id and existing_transcript_id != transcript_id:
                            # Keep track of different transcripts and their exon positions
                            if not existing_info_dict.get('other_transcripts'):
                                existing_info_dict['other_transcripts'] = {}
                            
                            # Add this transcript's info to the other_transcripts field
                            existing_info_dict['other_transcripts'][transcript_id] = {
                                'position': exon_position,
                                'exon_number': exon_number,
                                'transcript_name': transcript_data.get('transcript_name')
                            }
                            
                            # Update JSON string with enriched info
                            ensembl_info_str = json.dumps(existing_info_dict)
                        
                    except json.JSONDecodeError:
                        pass  # Use the new JSON if old one is corrupt
                
                # Update exon record with the new or merged information
                cursor.execute(
                    """
                    UPDATE exons SET 
                    phase = CASE WHEN ? IS NOT NULL THEN ? ELSE phase END, 
                    end_phase = CASE WHEN ? IS NOT NULL THEN ? ELSE end_phase END,
                    frame_status = CASE WHEN ? IS NOT NULL THEN ? ELSE frame_status END,
                    ensembl_info = ?
                    WHERE exon_id = ?
                    """,
                    (phase, phase, end_phase, end_phase, 
                     frame_status, frame_status, ensembl_info_str, existing_exon_id)
                )
            else:
                # Insert new exon record
                cursor.execute(
                    """
                    INSERT INTO exons 
                    (repeat_id, ensembl_exon_id, phase, end_phase, frame_status, ensembl_info) 
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (repeat_id, ensembl_exon_id, phase, end_phase, frame_status, ensembl_info_str)
                )

def populate_database(json_file, db_file, schema_file):
    """Populate the database with data from JSON file"""
    # Check if files exist
    if not os.path.exists(json_file):
        print(f"Error: JSON file not found: {json_file}")
        return False
        
    if not os.path.exists(schema_file):
        print(f"Error: Schema file not found: {schema_file}")
        return False
    
    # Remove existing database if it exists
    if os.path.exists(db_file):
        os.remove(db_file)
        print(f"Removed existing database: {db_file}")
    
    # Read JSON data
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            print(f"Loaded {len(data)} items from JSON file")
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON file: {e}")
        return False
    
    # Connect to database and create tables
    conn = sqlite3.connect(db_file)
    conn.execute("PRAGMA foreign_keys = ON")
    cursor = conn.cursor()
    
    try:
        create_tables(conn, schema_file)
        
        # Track processed genes and proteins to avoid duplicates
        processed_genes = {}
        processed_proteins = {}
        
        # Process each item
        for i, item in enumerate(data):
            if not item:  # Skip empty items
                continue
                
            if i % 100 == 0:
                print(f"Processing item {i+1}/{len(data)}")
            
            # Extract gene information
            gene_name = item.get('geneName')
            if not gene_name:
                continue
                
            # Get chromosome from chrom field
            chromosome = item.get('chrom', '').replace('chr', '')
            
            # Insert or get gene
            if gene_name in processed_genes:
                gene_id = processed_genes[gene_name]
            else:
                gene_id = insert_or_get_gene(cursor, gene_name, chromosome)
                processed_genes[gene_name] = gene_id
                
                # Insert gene aliases
                aliases = item.get('aliases', [])
                insert_gene_aliases(cursor, gene_id, aliases)
            
            # Extract protein information
            uniprot_id = item.get('uniProtId')
            if not uniprot_id:
                continue
                
            status = item.get('status', '')
            
            # Insert or get protein
            if uniprot_id in processed_proteins:
                protein_id = processed_proteins[uniprot_id]
            else:
                protein_id = insert_or_get_protein(cursor, uniprot_id, gene_id, status)
                processed_proteins[uniprot_id] = protein_id
            
            # Insert repeat
            repeat_id = insert_repeat(cursor, protein_id, item)
            
            # Process ensembl exon information
            ensembl_info = item.get('ensembl_exon_info')
            if ensembl_info:
                process_ensembl_info(cursor, repeat_id, ensembl_info, gene_id)
        
        # Commit changes
        conn.commit()
        print(f"Database populated successfully with data from {len(processed_genes)} genes and {len(processed_proteins)} proteins")
        
        # Create example query for exon skipping subjects
        print("\nExample SQL query for exon skipping subjects:")
        example_query = """
        -- Find genes with more than 5 of the same repeat type and exons with significant overlap
        SELECT g.gene_name, r.repeat_type, COUNT(DISTINCT r.repeat_id) as repeat_count,
               e.ensembl_exon_id, e.frame_status, e.ensembl_info
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        JOIN exons e ON r.repeat_id = e.repeat_id
        WHERE e.frame_status = 'in_frame'
          AND JSON_EXTRACT(e.ensembl_info, '$.overlap_percentage') >= 70
        GROUP BY g.gene_id, r.repeat_type
        HAVING COUNT(DISTINCT r.repeat_id) > 5
        ORDER BY g.gene_name, repeat_count DESC;
        """
        print(example_query)
        
        return True
        
    except sqlite3.Error as e:
        print(f"SQLite error: {e}")
        conn.rollback()
        return False
    finally:
        conn.close()

def main():
    json_file = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\1000_test_exons_hg38_repeats.json'
    db_file = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\repeats.db'
    schema_file = 'c:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\test_sqlite\\database_schema.sql'
    
    populate_database(json_file, db_file, schema_file)

if __name__ == "__main__":
    main()