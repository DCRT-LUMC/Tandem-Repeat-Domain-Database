import sqlite3
import os
import pandas as pd

def run_test_queries():
    """Run a series of test queries on the tandem repeat database"""
    # Connect to the database
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_file = os.path.join(script_dir, "tandem_repeats.db")
    
    print(f"Connecting to database: {db_file}")
    conn = sqlite3.connect(db_file)
    
    # Helper function to run a query and print results
    def run_query(query, title):
        print(f"\n=== {title} ===")
        try:
            result = pd.read_sql_query(query, conn)
            if len(result) > 0:
                print(result)
                print(f"Total rows: {len(result)}")
            else:
                print("No results found.")
        except Exception as e:
            print(f"Error executing query: {str(e)}")
    
    # Basic count queries
    run_query("SELECT COUNT(*) as count FROM genes", "Count of Genes")
    run_query("SELECT COUNT(*) as count FROM proteins", "Count of Proteins")
    run_query("SELECT COUNT(*) as count FROM repeats", "Count of Repeat Domains")
    run_query("SELECT COUNT(*) as count FROM transcripts", "Count of Transcripts")
    run_query("SELECT COUNT(*) as count FROM exons", "Count of Exons")
    
    # Sample data from each table
    run_query("SELECT * FROM genes LIMIT 5", "Sample Genes")
    run_query("SELECT * FROM proteins LIMIT 5", "Sample Proteins")
    run_query("SELECT * FROM repeats LIMIT 5", "Sample Repeats")
    run_query("SELECT * FROM transcripts LIMIT 5", "Sample Transcripts")
    run_query("SELECT * FROM exons LIMIT 5", "Sample Exons")
    
    # Join queries
    run_query("""
        SELECT g.gene_name, p.protein_id, r.repeat_type, r.start_pos, r.end_pos
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        LIMIT 10
    """, "Genes with Proteins and Repeats")
    
    # Find genes with the most repeat domains
    run_query("""
        SELECT g.gene_name, COUNT(r.repeat_id) as repeat_count
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        GROUP BY g.gene_name
        ORDER BY repeat_count DESC
        LIMIT 10
    """, "Genes with Most Repeat Domains")
    
    # Find the most common repeat types
    run_query("""
        SELECT repeat_type, COUNT(*) as count
        FROM repeats
        GROUP BY repeat_type
        ORDER BY count DESC
        LIMIT 10
    """, "Most Common Repeat Types")
    
    # Exons that are associated with repeat domains
    run_query("""
        SELECT e.exon_id, g.gene_name, re.overlap_bp, re.overlap_percentage
        FROM exons e
        JOIN genes g ON e.gene_id = g.gene_id
        JOIN repeat_exons re ON e.exon_id = re.exon_id
        JOIN repeats r ON re.repeat_id = r.repeat_id
        ORDER BY re.overlap_percentage DESC
        LIMIT 10
    """, "Exons with High Overlap to Repeat Domains")
    
    # Transcripts with the most exons
    run_query("""
        SELECT t.transcript_id, g.gene_name, COUNT(te.exon_id) as exon_count
        FROM transcripts t
        JOIN genes g ON t.gene_id = g.gene_id
        JOIN transcript_exons te ON t.transcript_id = te.transcript_id
        GROUP BY t.transcript_id
        ORDER BY exon_count DESC
        LIMIT 10
    """, "Transcripts with Most Exons")
    
    # Frame-preserving exons
    run_query("""
        SELECT g.gene_name, COUNT(e.exon_id) as frame_preserving_count
        FROM genes g
        JOIN exons e ON g.gene_id = e.gene_id
        WHERE e.frame_preserving = 1
        GROUP BY g.gene_name
        ORDER BY frame_preserving_count DESC
        LIMIT 10
    """, "Genes with Most Frame-Preserving Exons")
    
    conn.close()
    print("\nDatabase queries completed.")

def search_gene(gene_name="MTOR"):
    """Search for all information about a specific gene"""
    # Connect to the database
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_file = os.path.join(script_dir, "tandem_repeats.db")
    
    print(f"\nSearching for gene: {gene_name}")
    conn = sqlite3.connect(db_file)
    
    # Helper function to run a query and print results
    def run_query(query, title):
        print(f"\n=== {title} ===")
        try:
            result = pd.read_sql_query(query, conn)
            if len(result) > 0:
                print(result)
                print(f"Total rows: {len(result)}")
            else:
                print("No results found.")
            return result
        except Exception as e:
            print(f"Error executing query: {str(e)}")
            return pd.DataFrame()
    
    # 1. Basic gene information
    gene_info = run_query(f"""
        SELECT * FROM genes WHERE gene_name = '{gene_name}'
    """, f"Gene Information for {gene_name}")
    
    if gene_info.empty:
        print(f"Gene '{gene_name}' not found in database.")
        conn.close()
        return
        
    gene_id = gene_info.iloc[0]['gene_id']
    
    # 2. Protein information
    run_query(f"""
        SELECT * FROM proteins WHERE gene_id = {gene_id}
    """, f"Proteins for {gene_name}")
    
    # 3. Repeat domains in this gene
    run_query(f"""
        SELECT r.* 
        FROM repeats r
        JOIN proteins p ON r.protein_id = p.protein_id
        WHERE p.gene_id = {gene_id}
    """, f"Repeat Domains in {gene_name}")
    
    # 4. Transcripts for this gene
    run_query(f"""
        SELECT * FROM transcripts WHERE gene_id = {gene_id}
    """, f"Transcripts for {gene_name}")
    
    # 5. Exons for this gene
    run_query(f"""
        SELECT * FROM exons WHERE gene_id = {gene_id}
    """, f"Exons for {gene_name}")
    
    # 6. Comprehensive view of repeats, transcripts and exons
    run_query(f"""
        SELECT 
            g.gene_name,
            p.protein_id,
            r.repeat_id, 
            r.repeat_type,
            r.start_pos,
            r.end_pos,
            t.transcript_id,
            e.exon_id,
            re.overlap_bp,
            re.overlap_percentage
        FROM genes g
        JOIN proteins p ON g.gene_id = p.gene_id
        JOIN repeats r ON p.protein_id = r.protein_id
        LEFT JOIN repeat_transcripts rt ON r.repeat_id = rt.repeat_id
        LEFT JOIN transcripts t ON rt.transcript_id = t.transcript_id
        LEFT JOIN repeat_exons re ON r.repeat_id = re.repeat_id
        LEFT JOIN exons e ON re.exon_id = e.exon_id
        WHERE g.gene_name = '{gene_name}'
    """, f"Comprehensive View of {gene_name} Repeat Domains")
    
    conn.close()
    print(f"\nSearch for gene '{gene_name}' completed.")

def find_exon_skipping_candidates():
    """Find exons that are good candidates for skipping experiments and export to file"""
    # Connect to the database
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_file = os.path.join(script_dir, "tandem_repeats.db")
    output_file = os.path.join(script_dir, "exon_skipping_candidates.csv")
    
    conn = sqlite3.connect(db_file)
    
    print("\n=== Finding Exon Skipping Candidates ===")
    
    # First, find genes with transcripts containing 5+ repeats
    gene_transcript_query = """
        SELECT g.gene_id, g.gene_name, t.transcript_id, 
               COUNT(DISTINCT rt.repeat_id) as repeat_count
        FROM genes g
        JOIN transcripts t ON g.gene_id = t.gene_id
        JOIN repeat_transcripts rt ON t.transcript_id = rt.transcript_id
        GROUP BY g.gene_id, t.transcript_id
        HAVING repeat_count >= 5
    """
    
    # Get detailed information about exons in these transcripts
    exon_detail_query = """
        SELECT 
            g.gene_name, 
            g.chromosome,
            g.location,
            t.transcript_id, 
            e.exon_id,
            te.exon_number,
            e.length,
            e.frame_preserving,
            r.repeat_id,
            r.repeat_type,
            r.start_pos as repeat_start,
            r.end_pos as repeat_end,
            re.overlap_bp,
            re.overlap_percentage
        FROM genes g
        JOIN transcripts t ON g.gene_id = t.gene_id
        JOIN transcript_exons te ON t.transcript_id = te.transcript_id
        JOIN exons e ON te.exon_id = e.exon_id
        JOIN repeat_exons re ON e.exon_id = re.exon_id
        JOIN repeats r ON re.repeat_id = r.repeat_id
        WHERE e.frame_preserving = 1 
          AND re.overlap_percentage > 70
          AND t.transcript_id IN (
              SELECT transcript_id FROM (
                  SELECT t.transcript_id, COUNT(DISTINCT rt.repeat_id) as repeat_count
                  FROM transcripts t
                  JOIN repeat_transcripts rt ON t.transcript_id = rt.transcript_id
                  GROUP BY t.transcript_id
                  HAVING repeat_count >= 5
              )
          )
        ORDER BY g.gene_name, t.transcript_id, te.exon_number
    """
    
    try:
        # Get the results
        result = pd.read_sql_query(exon_detail_query, conn)
        
        if len(result) > 0:
            # Save to CSV file
            result.to_csv(output_file, index=False)
            
            # Print summary to console
            print(f"Found {len(result)} candidate exons for skipping experiments")
            print(f"Results saved to: {output_file}")
            
            # Summary statistics
            print("\n=== Summary ===")
            gene_count = result['gene_name'].nunique()
            transcript_count = result['transcript_id'].nunique()
            exon_count = result['exon_id'].nunique()
            print(f"Total genes: {gene_count}")
            print(f"Total transcripts: {transcript_count}")
            print(f"Total candidate exons: {exon_count}")
            
            # Show top 10 rows as preview
            print("\n=== Preview (top 10 rows) ===")
            pd.set_option('display.max_columns', None)
            print(result.head(10))
            
            # Group results by gene for better readability
            print("\n=== Summary by Gene ===")
            gene_groups = result.groupby('gene_name')
            for gene_name, group in gene_groups:
                transcript_count = group['transcript_id'].nunique()
                exon_count = group['exon_id'].nunique()
                max_overlap = group['overlap_percentage'].max()
                print(f"Gene {gene_name}: {transcript_count} transcripts with {exon_count} candidate exons (max overlap: {max_overlap:.1f}%)")
                
        else:
            print("No suitable exon skipping candidates found.")
    except Exception as e:
        print(f"Error executing query: {str(e)}")
    
    conn.close()

if __name__ == "__main__":
    # Uncomment to run
    run_test_queries()
    search_gene("MTOR")
    find_exon_skipping_candidates()
