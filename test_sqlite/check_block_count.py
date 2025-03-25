import sqlite3

def check_block_count_values(db_path="test_sqlite/repeats.db", limit=10):
    """Simple function to display block_count and repeat_length examples"""
    
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    # Query to get sample values
    cursor.execute("""
        SELECT repeat_id, repeat_type, block_count, repeat_length, protein_id
        FROM repeats
        WHERE block_count IS NOT NULL AND repeat_length IS NOT NULL
        LIMIT ?
    """, (limit,))
    
    results = cursor.fetchall()
    
    print("Block Count vs Repeat Length Examples:")
    print("=====================================")
    print(f"{'ID':<6} {'Type':<15} {'Block Count':<12} {'Repeat Len':<12}")
    print("-" * 50)
    
    for row in results:
        print(f"{row['repeat_id']:<6} {row['repeat_type'][:15]:<15} {row['block_count']:<12} {row['repeat_length']:<12}")
    
    # Also check if there are ANY matches at all
    cursor.execute("""
        SELECT COUNT(*) as match_count
        FROM repeats
        WHERE block_count = repeat_length
    """)
    match_count = cursor.fetchone()['match_count']
    
    print(f"\nTotal repeats where block_count = repeat_length: {match_count}")
    conn.close()

if __name__ == "__main__":
    check_block_count_values()