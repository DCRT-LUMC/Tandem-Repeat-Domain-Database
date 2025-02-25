import csv
import json
import re

# This converts the hg38_repeats.txt into a usable .json

input_file = 'data/hg38_repeats.txt'
output_file = 'data/hg38_repeats.json'

json_data = []

# Read the TSV file and adapt to its structure
with open(input_file, 'r') as txt_file:
    # Read the first line to determine the file structure
    first_line = txt_file.readline().strip()
    txt_file.seek(0)  # Reset file pointer to the beginning
    
    # Check if file has headers
    if first_line.startswith('#'):
        # Skip the header line if it starts with #
        next(txt_file)
        headers = first_line[1:].split('\t')
    else:
        # Try to determine if the first line contains headers or data
        first_line_parts = first_line.split('\t')
        try:
            # If first column can be converted to int, it's likely data not headers
            int(first_line_parts[0])
            # Use default column names for repeat data
            headers = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 
                      'period', 'copyNum', 'consensusSize', 'perMatch', 
                      'perIndel', 'alignment']
        except ValueError:
            # First line is probably headers
            headers = first_line_parts
            txt_file.seek(0)  # Reset file pointer to read headers as first row
    
    # Create a CSV reader with the detected or provided headers
    reader = csv.DictReader(txt_file, fieldnames=headers, delimiter='\t')
    
    # Process each row according to detected format
    for row in reader:
        # Convert row to dictionary
        entry = {header: row[header] for header in headers}
        
        # Try to convert numeric fields to appropriate types
        for field in ['chromStart', 'chromEnd', 'score', 'period', 'copyNum',
                     'consensusSize', 'perMatch', 'perIndel']:
            if field in entry and entry[field]:
                try:
                    if '.' in entry[field]:
                        entry[field] = float(entry[field])
                    else:
                        entry[field] = int(entry[field])
                except (ValueError, KeyError):
                    pass  # Keep as string if conversion fails
        
        json_data.append(entry)

# Write the JSON data to the output file
with open(output_file, 'w') as json_file:
    json.dump(json_data, json_file, indent=4)

print(f"Data has been successfully converted to {output_file}")