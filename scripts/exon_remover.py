#!/usr/bin/env python3
"""
Exon Removal Script for FASTA files
Removes specified exons and generates sequences with context exons
"""

import argparse
import sys
from pathlib import Path

def parse_fasta(file_path):
    """Parse FASTA file and return list of (header, sequence) tuples"""
    exons = []
    with open(file_path, 'r') as f:
        header = None
        sequence = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    exons.append((header, sequence))
                header = line[1:]  # Remove '>' character
                sequence = ""
            else:
                sequence += line
        
        # Add the last exon
        if header is not None:
            exons.append((header, sequence))
    
    return exons

def extract_exon_number(header):
    """Extract exon number from header (assumes ENSE ID format)"""
    # Look for ENSE followed by numbers
    parts = header.split()
    for part in parts:
        if part.startswith('ENSE'):
            return part
    return None

def generate_sequences(exons, target_exon_num, output_prefix, full_sequence=False):
    """Generate sequences with and without target exon"""
    total_exons = len(exons)
    
    # Convert to 0-based index
    target_index = target_exon_num - 1
    
    if target_index < 0 or target_index >= total_exons:
        print(f"Error: Exon {target_exon_num} not found. File has {total_exons} exons.")
        return
    
    # Calculate range (2 before, target, 2 after)
    start_index = max(0, target_index - 2)
    end_index = min(total_exons - 1, target_index + 2)
    
    # Get exons in range
    context_exons = exons[start_index:end_index + 1]
    
    print(f"\nProcessing exon {target_exon_num} (index {target_index})")
    print(f"Context range: exons {start_index + 1} to {end_index + 1}")
    print(f"Total exons in context: {len(context_exons)}")
    
    # Generate sequence WITH target exon (complete context)
    with_exon_seq = ""
    with_exon_headers = []
    
    for header, seq in context_exons:
        with_exon_seq += seq
        exon_id = extract_exon_number(header)
        if exon_id:
            with_exon_headers.append(exon_id)
    
    # Generate sequence WITHOUT target exon
    without_exon_seq = ""
    without_exon_headers = []
    
    for i, (header, seq) in enumerate(context_exons):
        # Skip the target exon (it's at position target_index - start_index)
        if start_index + i != target_index:
            without_exon_seq += seq
            exon_id = extract_exon_number(header)
            if exon_id:
                without_exon_headers.append(exon_id)
    
    # Generate full sequence WITHOUT target exon if requested
    if full_sequence:
        full_without_exon_seq = ""
        full_without_exon_headers = []
        
        for i, (header, seq) in enumerate(exons):
            # Skip the target exon
            if i != target_index:
                full_without_exon_seq += seq
                exon_id = extract_exon_number(header)
                if exon_id:
                    full_without_exon_headers.append(exon_id)
    
    # Write output files
    gene_name = exons[0][0].split()[0] if exons else "unknown"
    
    # File with exon included
    with_file = f"{output_prefix}_with_exon_{target_exon_num}.fa"
    with open(with_file, 'w') as f:
        f.write(f">{gene_name}_with_exon_{target_exon_num} exons_{'+'.join(with_exon_headers)}\n")
        # Write sequence in 60-character lines
        for i in range(0, len(with_exon_seq), 60):
            f.write(with_exon_seq[i:i+60] + '\n')
    
    # File without exon
    without_file = f"{output_prefix}_without_exon_{target_exon_num}.fa"
    with open(without_file, 'w') as f:
        f.write(f">{gene_name}_without_exon_{target_exon_num} exons_{'+'.join(without_exon_headers)}\n")
        # Write sequence in 60-character lines
        for i in range(0, len(without_exon_seq), 60):
            f.write(without_exon_seq[i:i+60] + '\n')
    
    # Write full sequence file if requested
    if full_sequence:
        full_file = f"{output_prefix}_full_without_exon_{target_exon_num}.fa"
        with open(full_file, 'w') as f:
            f.write(f">{gene_name}_full_without_exon_{target_exon_num} all_exons_except_{target_exon_num}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(full_without_exon_seq), 60):
                f.write(full_without_exon_seq[i:i+60] + '\n')
    
    print(f"\nOutput files generated:")
    print(f"  With exon {target_exon_num}: {with_file}")
    print(f"  Without exon {target_exon_num}: {without_file}")
    if full_sequence:
        print(f"  Full sequence without exon {target_exon_num}: {full_file}")
    
    print(f"\nSequence lengths:")
    print(f"  With exon: {len(with_exon_seq)} bp")
    print(f"  Without exon: {len(without_exon_seq)} bp")
    if full_sequence:
        print(f"  Full without exon: {len(full_without_exon_seq)} bp")
    
    print(f"\nExon IDs included:")
    print(f"  With exon: {', '.join(with_exon_headers)}")
    print(f"  Without exon: {', '.join(without_exon_headers)}")
    if full_sequence:
        print(f"  Full without exon: {len(full_without_exon_headers)} exons total")
    
    # Get the removed exon sequence
    removed_exon_header, removed_exon_seq = exons[target_index]
    removed_exon_id = extract_exon_number(removed_exon_header)
    
    # Print the actual FASTA sequences to console
    print(f"\n{'='*60}")
    print("REMOVED EXON SEQUENCE:")
    print(f"{'='*60}")
    print(f">{gene_name}_exon_{target_exon_num} {removed_exon_id}")
    for i in range(0, len(removed_exon_seq), 60):
        print(removed_exon_seq[i:i+60])
    
    print(f"\n{'='*60}")
    print("FASTA SEQUENCE WITH TARGET EXON:")
    print(f"{'='*60}")
    print(f">{gene_name}_with_exon_{target_exon_num} exons_{'+'.join(with_exon_headers)}")
    for i in range(0, len(with_exon_seq), 60):
        print(with_exon_seq[i:i+60])
    
    print(f"\n{'='*60}")
    print("FASTA SEQUENCE WITHOUT TARGET EXON:")
    print(f"{'='*60}")
    print(f">{gene_name}_without_exon_{target_exon_num} exons_{'+'.join(without_exon_headers)}")
    for i in range(0, len(without_exon_seq), 60):
        print(without_exon_seq[i:i+60])
    
    # Print full sequence if requested
    if full_sequence:
        print(f"\n{'='*60}")
        print("FULL SEQUENCE WITHOUT TARGET EXON:")
        print(f"{'='*60}")
        print(f">{gene_name}_full_without_exon_{target_exon_num} all_exons_except_{target_exon_num}")
        for i in range(0, len(full_without_exon_seq), 60):
            print(full_without_exon_seq[i:i+60])
    
    print(f"{'='*60}")

def main():
    parser = argparse.ArgumentParser(
        description="Remove specific exons from FASTA files and generate context sequences",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python exon_remover.py input.fa 7
  python exon_remover.py input.fa 7 --output custom_name
  python exon_remover.py input.fa 7 --list-exons
  python exon_remover.py input.fa 7 --full
        """
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file path')
    parser.add_argument('exon_number', type=int, help='Exon number to remove (1-based)')
    parser.add_argument('--output', '-o', help='Output file prefix (default: input filename without extension)')
    parser.add_argument('--list-exons', '-l', action='store_true', help='List all exons in the file')
    parser.add_argument('--full', '-f', action='store_true', help='Generate full sequence with target exon removed')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.fasta_file).exists():
        print(f"Error: File '{args.fasta_file}' not found.")
        sys.exit(1)
    
    # Parse FASTA file
    try:
        exons = parse_fasta(args.fasta_file)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    if not exons:
        print("Error: No exons found in FASTA file.")
        sys.exit(1)
    
    print(f"Loaded {len(exons)} exons from {args.fasta_file}")
    
    # List exons if requested
    if args.list_exons:
        print("\nExons found:")
        for i, (header, seq) in enumerate(exons, 1):
            exon_id = extract_exon_number(header)
            print(f"  {i}: {exon_id} ({len(seq)} bp)")
        print()
    
    # Set output prefix
    if args.output:
        output_prefix = args.output
    else:
        output_prefix = Path(args.fasta_file).stem
    
    # Generate sequences
    generate_sequences(exons, args.exon_number, output_prefix, args.full)

if __name__ == "__main__":
    main()
