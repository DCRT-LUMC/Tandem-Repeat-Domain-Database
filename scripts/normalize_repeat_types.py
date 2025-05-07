#!/usr/bin/env python3
import json
import re
import sys
import os
import glob

def normalize_repeat_type(repeat_type):
    if not repeat_type:
        return "Unknown"
    repeat_type = repeat_type.strip()
    if repeat_type.startswith("WD"):
        return "WD"
    if repeat_type.startswith("Alpha"):
        return "Alpha"
    if repeat_type.startswith("PUR"):
        return "PUR"
    if repeat_type.startswith("RCC1"):
        return "RCC1"
    if repeat_type.startswith("SVP"):
        return "SVP"
    if repeat_type.startswith("Spectrin"):
        return "Spectrin"
    # If starts with a number
    if repeat_type and repeat_type[0].isdigit():
        return "Unknown"
    # If the entire string is a roman numeral I-X (case-insensitive)
    if re.match(r"^(I|II|III|IV|V|VI|VII|VIII|IX|X)$", repeat_type, re.IGNORECASE):
        return "Unknown"
    return repeat_type

def normalize_json_file(input_file, output_file=None):
    """
    Normalize repeat types in a JSON file.
    
    Args:
        input_file (str): Path to the input JSON file
        output_file (str, optional): Path to the output JSON file. If None, overwrites the input file.
    """
    # If output_file is not specified, use the input_file
    if output_file is None:
        output_file = input_file
    
    try:
        # Read the JSON file
        with open(input_file, 'r') as f:
            data = json.load(f)
        
        # Ensure the data is a list (array)
        if not isinstance(data, list):
            print(f"Error: Expected a JSON array, but got {type(data).__name__}")
            return False
        
        # Normalize each repeat type
        for item in data:
            if "repeatType" in item:
                item["repeatType"] = normalize_repeat_type(item["repeatType"])
        
        # Write the modified JSON back to the file
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"Successfully normalized repeat types in {input_file}")
        if output_file != input_file:
            print(f"Saved to {output_file}")
        
        return True
    
    except Exception as e:
        print(f"Error: {e}")
        return False

def normalize_directory(input_dir, pattern="*_annotated_repeats.json"):
    """
    Normalize repeat types in all JSON files in a directory that match a pattern.
    
    Args:
        input_dir (str): Path to the directory containing JSON files
        pattern (str, optional): Glob pattern to match files
    """
    # Get all files matching the pattern
    files = glob.glob(os.path.join(input_dir, pattern))
    
    success_count = 0
    for file in files:
        print(f"Processing {file}...")
        if normalize_json_file(file):
            success_count += 1
    
    print(f"Normalized {success_count} out of {len(files)} files")

if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python normalize_repeat_types.py <input_file_or_directory> [output_file]")
        sys.exit(1)
    
    input_path = sys.argv[1]
    
    # If input is a directory, process all JSON files in it
    if os.path.isdir(input_path):
        normalize_directory(input_path)
    # If input is a file, process just that file
    elif os.path.isfile(input_path):
        output_file = sys.argv[2] if len(sys.argv) > 2 else None
        normalize_json_file(input_path, output_file)
    else:
        print(f"Error: {input_path} is not a valid file or directory")
        sys.exit(1)
