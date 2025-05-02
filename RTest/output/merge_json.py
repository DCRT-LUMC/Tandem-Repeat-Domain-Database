import json
import os
import glob
import re

def merge_json_files(directory_path, output_filename):
    # Get all JSON files in the directory
    json_files = glob.glob(os.path.join(directory_path, '*_annotated_repeats.json'))
    
    # Sort files based on the numerical range they contain
    def get_start_number(file_path):
        filename = os.path.basename(file_path)
        # Extract numbers from the filename using regex
        match = re.search(r'(\d+)-\d+_', filename)
        if match:
            return int(match.group(1))
        return 0  # Default value if no match is found
    
    json_files.sort(key=get_start_number)  # Sort files by numerical order
    
    # Initialize an empty list or dict to store the merged data
    merged_data = None
    
    # Process each file
    for file_path in json_files:
        with open(file_path, 'r') as file:
            try:
                data = json.load(file)
                print(f"Processing {os.path.basename(file_path)}...")
                
                # Initialize merged_data based on the first file
                if merged_data is None:
                    merged_data = [] if isinstance(data, list) else {}
                
                # Merge data appropriately based on type
                if isinstance(data, list) and isinstance(merged_data, list):
                    merged_data.extend(data)
                elif isinstance(data, dict) and isinstance(merged_data, dict):
                    merged_data.update(data)
                else:
                    print(f"Warning: Data structure mismatch in {file_path}")
            except json.JSONDecodeError:
                print(f"Error: Could not parse JSON from {file_path}")
    
    # Write the merged data to a new file
    if merged_data is not None:
        output_path = os.path.join(directory_path, output_filename)
        with open(output_path, 'w') as output_file:
            json.dump(merged_data, output_file, indent=2)
        print(f"Merged JSON saved to {output_path}")
    else:
        print("No valid data found to merge")

if __name__ == "__main__":
    directory = "/home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/RTest/output/canonical_v2"
    merge_json_files(directory, "v2_all_cannonical_repeats_annotated.json")