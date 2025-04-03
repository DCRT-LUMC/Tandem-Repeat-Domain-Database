#!/usr/bin/env python3

import sys
import os
import math
from collections import defaultdict

def extract_sequence_from_pdb(pdb_file):
    """Extract amino acid sequence from SEQRES records in a PDB file."""
    if not os.path.exists(pdb_file):
        print(f"Error: File {pdb_file} does not exist.")
        return None
    
    sequence = ""
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("SEQRES"):
                    parts = line.split()
                    for aa in parts[4:]:
                        if aa in three_to_one:
                            sequence += three_to_one[aa]
    except Exception as e:
        print(f"Error reading file {pdb_file}: {e}")
        return None
        
    if not sequence:
        # Try to extract from ATOM records if SEQRES is not available
        residue_dict = {}
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith("ATOM"):
                        resnum = int(line[22:26].strip())
                        resname = line[17:20].strip()
                        if resnum not in residue_dict and resname in three_to_one:
                            residue_dict[resnum] = three_to_one[resname]
            
            if residue_dict:
                # Sort by residue number and join
                sequence = ''.join([residue_dict[k] for k in sorted(residue_dict.keys())])
        except Exception as e:
            print(f"Error extracting sequence from ATOM records in {pdb_file}: {e}")
    
    return sequence

def find_sequence_differences(seq1, seq2):
    """Find positions where two sequences differ."""
    differences = []
    
    min_length = min(len(seq1), len(seq2))
    
    # Compare aligned regions
    for i in range(min_length):
        if seq1[i] != seq2[i]:
            differences.append({
                'position': i + 1,  # 1-based position
                'seq1_aa': seq1[i],
                'seq2_aa': seq2[i]
            })
    
    # Handle different lengths
    if len(seq1) > len(seq2):
        for i in range(min_length, len(seq1)):
            differences.append({
                'position': i + 1,
                'seq1_aa': seq1[i],
                'seq2_aa': '-'  # Gap in seq2
            })
    elif len(seq2) > len(seq1):
        for i in range(min_length, len(seq2)):
            differences.append({
                'position': i + 1,
                'seq1_aa': '-',  # Gap in seq1
                'seq2_aa': seq2[i]
            })
    
    return differences

def compare_sequences(seq1, seq2):
    """Compare two amino acid sequences and return detailed statistics."""
    if not seq1 or not seq2:
        return "Cannot compare sequences - one or both sequences are empty."
    
    length1 = len(seq1)
    length2 = len(seq2)
    
    # Calculate percentage identity
    min_length = min(length1, length2)
    matches = sum(a == b for a, b in zip(seq1[:min_length], seq2[:min_length]))
    identity = (matches / min_length) * 100
    
    # Find differences
    differences = find_sequence_differences(seq1, seq2)
    
    result = f"Sequence 1: {length1} amino acids\n"
    result += f"Sequence 2: {length2} amino acids\n"
    result += f"Sequence lengths differ by {abs(length1 - length2)} amino acids\n"
    result += f"Percentage identity over aligned region: {identity:.2f}%\n"
    result += f"Number of differences: {len(differences)}\n\n"
    
    if differences:
        result += "Differences (position: seq1 -> seq2):\n"
        for diff in differences:
            result += f"Position {diff['position']}: {diff['seq1_aa']} -> {diff['seq2_aa']}\n"
    else:
        result += "The sequences are identical.\n"
    
    # Find and report repeat patterns (basic implementation)
    result += "\nRepeat analysis (basic):\n"
    for seq_idx, seq in enumerate([seq1, seq2], 1):
        for repeat_length in range(2, 21):  # Check for repeats of length 2-20
            for i in range(len(seq) - 2*repeat_length):
                pattern = seq[i:i+repeat_length]
                if pattern == seq[i+repeat_length:i+2*repeat_length]:
                    result += f"Sequence {seq_idx}: Potential repeat of '{pattern}' at position {i+1}\n"
                    break
    
    return result

def visualize_sequence_alignment(seq1, seq2, window_size=80):
    """Create a visual representation of the sequence alignment with differences highlighted."""
    min_length = min(len(seq1), len(seq2))
    
    result = ""
    
    # Process the sequence in windows for better readability
    for start in range(0, min_length, window_size):
        end = min(start + window_size, min_length)
        
        # Create match line
        match_line = ""
        for i in range(start, end):
            if seq1[i] == seq2[i]:
                match_line += "|"
            else:
                match_line += " "
        
        result += f"Pos: {start+1}\n"
        result += f"Seq1: {seq1[start:end]}\n"
        result += f"      {match_line}\n"
        result += f"Seq2: {seq2[start:end]}\n\n"
    
    # Handle tails if sequences have different lengths
    if len(seq1) > min_length:
        result += f"Seq1 extra tail: {seq1[min_length:]}\n"
    if len(seq2) > min_length:
        result += f"Seq2 extra tail: {seq2[min_length:]}\n"
        
    return result

def main():
    if len(sys.argv) < 3:
        print("Usage: python compare_residues.py <pdb_file1> <pdb_file2>")
        sys.exit(1)
    
    pdb_file1 = sys.argv[1]
    pdb_file2 = sys.argv[2]
    
    if not os.path.exists(pdb_file1):
        print(f"Error: PDB file '{pdb_file1}' not found")
        sys.exit(1)
    
    if not os.path.exists(pdb_file2):
        print(f"Error: PDB file '{pdb_file2}' not found")
        sys.exit(1)
    
    print(f"Comparing sequences from {os.path.basename(pdb_file1)} and {os.path.basename(pdb_file2)}...")
    
    seq1 = extract_sequence_from_pdb(pdb_file1)
    seq2 = extract_sequence_from_pdb(pdb_file2)
    
    if seq1:
        print(f"\nSequence 1 ({os.path.basename(pdb_file1)}):")
        print(seq1[:50] + "..." if len(seq1) > 50 else seq1)  # Show just the beginning
    else:
        print(f"\nCould not extract sequence from {pdb_file1}")
    
    if seq2:
        print(f"\nSequence 2 ({os.path.basename(pdb_file2)}):")
        print(seq2[:50] + "..." if len(seq2) > 50 else seq2)  # Show just the beginning
    else:
        print(f"\nCould not extract sequence from {pdb_file2}")
    
    if seq1 and seq2:
        print("\nComparison results:")
        print(compare_sequences(seq1, seq2))
        
        print("\nSequence alignment visualization:")
        print(visualize_sequence_alignment(seq1, seq2))

if __name__ == "__main__":
    main()