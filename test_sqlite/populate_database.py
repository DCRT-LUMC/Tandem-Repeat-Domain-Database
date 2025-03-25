#!/usr/bin/env python3
# filepath: /home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/test_sqlite/populate_database.py
import json
import sqlite3
import os
import sys
from collections import defaultdict
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_tables(conn, schema_file):
    """Create database tables from schema file"""
    with open(schema_file, 'r') as f:
        schema_sql = f.read()
        conn.executescript(schema_sql)
    logger.info("Database tables created successfully")

def insert_or_get_gene(cursor, gene_name, chromosome=None, location=None, gene_type=None, ensembl_gene_id=None):
    """Insert gene if not exists and return gene_id"""
    cursor.execute(
        "SELECT gene_id FROM genes WHERE gene_name = ?", 
        (gene_name,)
    )
    result = cursor.fetchone()
    
    if result:
        return result[0]
    
    cursor.execute(
        """INSERT INTO genes 
           (gene_name, chromosome, location, gene_type, ensembl_gene_id) 
           VALUES (?, ?, ?, ?, ?)""",
        (gene_name, chromosome, location, gene_type, ensembl_gene_id)
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

def insert_or_get_protein(cursor, uniprot_id, gene_id, length=None, description=None, status=None):
    """Insert protein if not exists and return protein_id"""
    cursor.execute(
        "SELECT protein_id FROM proteins WHERE protein_id = ?", 
        (uniprot_id,)
    )
    result = cursor.fetchone()
    
    if result:
        return result[0]
    
    cursor.execute(
        """INSERT INTO proteins 
           (protein_id, gene_id, length, description, status) 
           VALUES (?, ?, ?, ?, ?)""",
        (uniprot_id, gene_id, length, description, status)
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
    
    # Extract start and end positions from position field
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
    
    # Add block data if available
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
            "UPDATE repeats SET block_count = ?, block_sizes = ?, block_starts = ? WHERE repeat_id = ?",
            (block_count, block_sizes_str, block_starts_str, repeat_id)
        )
    
    return repeat_id

def insert_or_get_transcript(cursor, transcript_data, gene_id):
    """Insert transcript if not exists and return transcript_id"""
    transcript_id = transcript_data.get('transcript_id')
    if not transcript_id:
        return None
        
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
        (transcript_id, gene_id, description, ensembl_transcript_id, 
         versioned_transcript_id, transcript_name, is_canonical, biotype, exon_count) 
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            transcript_id, 
            gene_id, 
            transcript_data.get('description'),
            transcript_data.get('ensembl_transcript_id'),
            transcript_data.get('versioned_transcript_id'),
            transcript_data.get('transcript_name'),
            transcript_data.get('is_canonical', 0),
            transcript_data.get('biotype'),
            transcript_data.get('exon_count')
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
        if not transcript_id:
            continue
        
        # Link repeat to transcript
        cursor.execute(
            """
            INSERT OR IGNORE INTO repeat_transcripts 
            (repeat_id, transcript_id, genomic_start, genomic_end, exon_mapping, location) 
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                repeat_id, 
                transcript_id, 
                transcript_data.get('genomic_start'),
                transcript_data.get('genomic_end'),
                json.dumps(transcript_data.get('exon_mapping', {})),
                transcript_data.get('location')
            )
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
            
            exon_number = exon_data.get('exon_number')
            exon_position = exon_data.get('position')
            
            # Prepare exon JSON info structure
            ensembl_info_dict = {
                'overlap_bp': exon_data.get('overlap_bp'),
                'position': exon_position,
                'overlap_percentage': exon_data.get('overlap_percentage'),
                'coding_status': exon_data.get('coding_status'),
                'utr_status': exon_data.get('utr_status'),
                'coding_percentage': exon_data.get('coding_percentage'),
                'exon_number': exon_number,
                'transcript_id': transcript_id,
                'transcript_name': transcript_data.get('transcript_name')
            }
            
            # Insert or update exon
            phase = exon_data.get('phase')
            end_phase = exon_data.get('end_phase')
            frame_status = exon_data.get('frame_status')
            coding_status = exon_data.get('coding_status')
            utr_status = exon_data.get('utr_status')
            coding_percentage = exon_data.get('coding_percentage')
            
            cursor.execute(
                "SELECT exon_id FROM exons WHERE ensembl_exon_id = ?",
                (ensembl_exon_id,)
            )
            exon_result = cursor.fetchone()
            
            if exon_result:
                # Update existing exon
                exon_id = exon_result[0]
                cursor.execute(
                    """
                    UPDATE exons SET 
                    phase = COALESCE(?, phase),
                    end_phase = COALESCE(?, end_phase),
                    frame_status = COALESCE(?, frame_status),
                    coding_status = COALESCE(?, coding_status),
                    utr_status = COALESCE(?, utr_status),
                    coding_percentage = COALESCE(?, coding_percentage)
                    WHERE exon_id = ?
                    """,
                    (phase, end_phase, frame_status, coding_status, 
                     utr_status, coding_percentage, exon_id)
                )
            else:
                # Insert new exon
                cursor.execute(
                    """
                    INSERT INTO exons 
                    (ensembl_exon_id, phase, end_phase, frame_status, 
                     coding_status, utr_status, coding_percentage) 
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (ensembl_exon_id, phase, end_phase, frame_status,
                     coding_status, utr_status, coding_percentage)
                )
                exon_id = cursor.lastrowid
            
            # Link exon to repeat
            cursor.execute(
                "INSERT OR IGNORE INTO repeat_exons (repeat_id, exon_id) VALUES (?, ?)",
                (repeat_id, exon_id)
            )
            
            # Link exon to transcript with position info
            cursor.execute(
                """
                INSERT OR IGNORE INTO transcript_exons 
                (transcript_id, exon_id, exon_number, overlap_bp, 
                 exon_position_in_transcript, overlap_percentage) 
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (transcript_id, exon_id, exon_number, 
                 ensembl_info_dict.get('overlap_bp'),
                 exon_position,
                 ensembl_info_dict.get('overlap_percentage'))
            )

def populate_database(json_file, db_file, schema_file):
    """Populate the database with data from JSON file"""
    # Check if files exist
    if not os.path.exists(json_file):
        logger.error(f"JSON file not found: {json_file}")
        return False
        
    if not os.path.exists(schema_file):
        logger.error(f"Schema file not found: {schema_file}")
        return False
    
    # Remove existing database if it exists
    if os.path.exists(db_file):
        os.remove(db_file)
        logger.info(f"Removed existing database: {db_file}")
    
    # Read JSON data
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)
            logger.info(f"Loaded {len(data)} items from JSON file")
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing JSON file: {e}")
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
                logger.info(f"Processing item {i+1}/{len(data)}")
            
            # Extract gene information
            gene_name = item.get('geneName')
            if not gene_name:
                continue
                
            # Get chromosome from chrom field
            chromosome = item.get('chrom', '').replace('chr', '')
            location = item.get('location')
            
            # Insert or get gene
            if gene_name in processed_genes:
                gene_id = processed_genes[gene_name]
            else:
                gene_id = insert_or_get_gene(cursor, gene_name, chromosome, location)
                processed_genes[gene_name] = gene_id
                
                # Insert gene aliases
                aliases = item.get('aliases', [])
                insert_gene_aliases(cursor, gene_id, aliases)
            
            # Extract protein information
            uniprot_id = item.get('uniProtId')
            if not uniprot_id:
                continue
                
            status = item.get('status', '')
            protein_length = item.get('length')
            protein_description = item.get('description')
            
            # Insert or get protein
            if uniprot_id in processed_proteins:
                protein_id = processed_proteins[uniprot_id]
            else:
                protein_id = insert_or_get_protein(
                    cursor, uniprot_id, gene_id, 
                    protein_length, protein_description, status
                )
                processed_proteins[uniprot_id] = protein_id
            
            # Insert repeat
            repeat_id = insert_repeat(cursor, protein_id, item)
            
            # Process ensembl exon information
            ensembl_info = item.get('ensembl_exon_info')
            if ensembl_info:
                process_ensembl_info(cursor, repeat_id, ensembl_info, gene_id)
        
        conn.commit()
        logger.info(f"Database populated successfully with data from {len(processed_genes)} genes and {len(processed_proteins)} proteins")
        
        # Create example query for exon skipping subjects
        logger.info("\nExample SQL query for exon skipping subjects:")
        example_query = """
        -- Find genes with more than 5 of the same repeat type and exons with significant overlap
        SELECT g.gene_name, r.repeat_type, COUNT(DISTINCT r.repeat_id) as repeat_count,
               e.ensembl_exon_id, e.frame_status, te.overlap_percentage
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        JOIN repeat_exons re ON r.repeat_id = re.repeat_id
        JOIN exons e ON re.exon_id = e.exon_id
        JOIN transcript_exons te ON e.exon_id = te.exon_id
        WHERE e.frame_status = 'in_frame'
          AND te.overlap_percentage >= 70
        GROUP BY g.gene_id, r.repeat_type
        HAVING COUNT(DISTINCT r.repeat_id) > 5
        ORDER BY g.gene_name, repeat_count DESC;
        """
        logger.info(example_query)
        
        return True
        
    except sqlite3.Error as e:
        logger.error(f"SQLite error: {e}")
        conn.rollback()
        return False
    finally:
        conn.close()

def main():
    json_file = 'test_sqlite/1000_test_exons_hg38_repeats.json'
    db_file = 'test_sqlite/repeats.db'
    schema_file = 'test_sqlite/database_schema.sql'
    
    if len(sys.argv) > 1:
        json_file = sys.argv[1]
    if len(sys.argv) > 2:
        db_file = sys.argv[2]
    if len(sys.argv) > 3:
        schema_file = sys.argv[3]
        
    success = populate_database(json_file, db_file, schema_file)
    if success:
        logger.info("Database population completed successfully")
    else:
        logger.error("Failed to populate database")

if __name__ == "__main__":
    main()