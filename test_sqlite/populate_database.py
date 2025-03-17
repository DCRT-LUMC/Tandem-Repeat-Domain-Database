import json
import sqlite3
import os
from tqdm import tqdm  # For progress bars

def populate_test_database(json_file, db_file):
    """Populate SQLite database from the repeat domain JSON file"""
    print(f"Creating database from {json_file}...")
    
    # Connect to SQLite DB
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Create tables using schema file
    # Get the directory of the script to find the schema file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    schema_file_path = os.path.join(script_dir, "database_schema.sql")
    
    with open(schema_file_path, "r") as schema_file:
        schema = schema_file.read()
        conn.executescript(schema)
    
    # Load JSON data
    with open(json_file, "r") as f:
        data = json.load(f)
    
    # Track processed IDs to avoid duplicates
    processed_genes = {}
    processed_proteins = {}
    processed_transcripts = {}
    processed_exons = {}

    # Process each repeat entry
    for repeat in tqdm(data, desc="Processing repeats"):
        if not isinstance(repeat, dict) or not repeat:
            continue  # Skip empty entries
        
        # Get gene info
        gene_name = repeat.get("geneName", "")
        if not gene_name:
            continue
        
        # Insert gene if not already processed
        if gene_name not in processed_genes:
            aliases = repeat.get("aliases", "")
            if isinstance(aliases, list):
                aliases = ",".join(aliases)
            
            cursor.execute("""
                INSERT INTO genes (gene_name, aliases, chromosome, location)
                VALUES (?, ?, ?, ?)
            """, (
                gene_name,
                aliases,
                repeat.get("chrom", ""),
                f"{repeat.get('chrom', '')}:{repeat.get('chromStart', '')}_{repeat.get('chromEnd', '')}"
            ))
            processed_genes[gene_name] = cursor.lastrowid
        
        gene_id = processed_genes[gene_name]
        
        # Process protein
        uniprot_id = repeat.get("uniProtId", "")
        if uniprot_id and uniprot_id not in processed_proteins:
            status = repeat.get("status", "")
            
            # Extract length from position if possible: "amino acids 343-389 on protein Q6TDP4"
            position = repeat.get("position", "")
            length = 0
            if position and isinstance(position, str):
                try:
                    parts = position.split()
                    if len(parts) >= 3:
                        pos = parts[2].split("-")
                        if len(pos) == 2:
                            length = int(pos[1]) - int(pos[0]) + 1
                except:
                    pass
            
            cursor.execute("""
                INSERT INTO proteins (protein_id, gene_id, length, description, status)
                VALUES (?, ?, ?, ?, ?)
            """, (
                uniprot_id,
                gene_id,
                length,
                f"Protein for {gene_name}",
                status
            ))
            processed_proteins[uniprot_id] = uniprot_id
        
        # Process repeat domain
        amino_start = None
        amino_end = None
        if repeat.get("position") and isinstance(repeat.get("position"), str):
            position = repeat.get("position")
            try:
                # Extract positions from "amino acids 343-389 on protein Q6TDP4"
                parts = position.split()
                if len(parts) >= 3:
                    pos = parts[2].split("-")
                    if len(pos) == 2:
                        amino_start = int(pos[0])
                        amino_end = int(pos[1])
            except:
                pass
        
        # Calculate sequence length from blockSizes
        sequence_length = 0
        block_sizes = repeat.get("blockSizes", [])
        if isinstance(block_sizes, list):
            for size in block_sizes:
                try:
                    sequence_length += int(size)
                except:
                    pass
        
        # Insert repeat
        cursor.execute("""
            INSERT INTO repeats (protein_id, repeat_type, start_pos, end_pos, sequence)
            VALUES (?, ?, ?, ?, ?)
        """, (
            uniprot_id,
            repeat.get("repeatType", ""),
            amino_start,
            amino_end,
            "N" * sequence_length  # Placeholder sequence of Ns
        ))
        repeat_id = cursor.lastrowid
        
        # Process exon information if available
        if "ensembl_exon_info" not in repeat:
            continue
            
        exon_info = repeat.get("ensembl_exon_info", {})
        if not exon_info or "transcripts" not in exon_info:
            continue
            
        # Process each transcript
        for transcript_data in exon_info.get("transcripts", []):
            transcript_id = transcript_data.get("transcript_id")
            if not transcript_id:
                continue
                
            # Insert transcript if not already processed
            if transcript_id not in processed_transcripts:
                cursor.execute("""
                    INSERT INTO transcripts (transcript_id, gene_id, description)
                    VALUES (?, ?, ?)
                """, (
                    transcript_id,
                    gene_id,
                    f"{transcript_data.get('transcript_name', '')} ({transcript_data.get('biotype', '')})"
                ))
                processed_transcripts[transcript_id] = transcript_id
            
            # Create repeat_transcript relationship
            genomic_start = repeat.get("chromStart", 0)
            genomic_end = repeat.get("chromEnd", 0)
            
            # Convert exon mapping to JSON string
            exon_mapping = json.dumps([
                {
                    "exon_id": exon.get("exon_id", ""),
                    "exon_number": exon.get("exon_number", 0),
                    "overlap_bp": exon.get("overlap_bp", 0),
                    "overlap_percentage": exon.get("overlap_percentage", 0),
                    "coding_percentage": exon.get("coding_percentage", 0)
                } for exon in transcript_data.get("containing_exons", [])
            ])
            
            cursor.execute("""
                INSERT INTO repeat_transcripts (
                    repeat_id, transcript_id, genomic_start, genomic_end, exon_mapping
                ) VALUES (?, ?, ?, ?, ?)
            """, (
                repeat_id, 
                transcript_id, 
                genomic_start, 
                genomic_end, 
                exon_mapping
            ))
            
            # Process exons in this transcript
            for exon_data in transcript_data.get("containing_exons", []):
                exon_id = exon_data.get("exon_id")
                if not exon_id or exon_id in processed_exons:
                    continue
                    
                # Estimate exon size from overlap percentage
                exon_size = 0
                if exon_data.get("overlap_percentage") and exon_data.get("overlap_bp"):
                    try:
                        exon_size = int(exon_data.get("overlap_bp") * 100 / exon_data.get("overlap_percentage"))
                    except:
                        pass
                
                # Calculate if skipping would preserve reading frame
                frame_preserving = exon_size % 3 == 0
                
                cursor.execute("""
                    INSERT INTO exons (
                        exon_id, gene_id, length, frame_preserving
                    ) VALUES (?, ?, ?, ?)
                """, (
                    exon_id,
                    gene_id,
                    exon_size,
                    frame_preserving
                ))
                processed_exons[exon_id] = exon_id
                
                # Create transcript_exon relationship
                cursor.execute("""
                    INSERT INTO transcript_exons (
                        transcript_id, exon_id, exon_number
                    ) VALUES (?, ?, ?)
                """, (
                    transcript_id,
                    exon_id,
                    exon_data.get("exon_number", 0)
                ))
                
                # Create repeat_exon relationship
                cursor.execute("""
                    INSERT INTO repeat_exons (
                        repeat_id, exon_id, overlap_bp, overlap_percentage
                    ) VALUES (?, ?, ?, ?)
                """, (
                    repeat_id,
                    exon_id,
                    exon_data.get("overlap_bp", 0),
                    exon_data.get("overlap_percentage", 0)
                ))
    
    # Commit all changes
    conn.commit()
    
    # Print some statistics
    cursor.execute("SELECT COUNT(*) FROM genes")
    gene_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM proteins")
    protein_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM repeats")
    repeat_count = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM exons")
    exon_count = cursor.fetchone()[0]
    
    print(f"\nDatabase populated successfully:")
    print(f"  - {gene_count} genes")
    print(f"  - {protein_count} proteins")
    print(f"  - {repeat_count} repeat domains")
    print(f"  - {exon_count} exons")
    
    conn.close()
    
    print(f"\nDatabase created at: {db_file}")

if __name__ == "__main__":
    # Define file paths - simplified to use the same directory as the script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    json_file = os.path.join(script_dir, "1000_test_exons_hg38_repeats.json")
    db_file = os.path.join(script_dir, "tandem_repeats.db")
    
    # Create the database
    populate_test_database(json_file, db_file)