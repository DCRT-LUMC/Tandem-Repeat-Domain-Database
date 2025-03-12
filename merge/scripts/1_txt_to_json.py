import csv
import json
import re

# This converts the hg38_repeats.txt into a usable .json

input_file = 'merge/data/hg38_repeats.txt'
output_file = 'merge/data/hg38_repeats.json'

json_data = []

# Fields that should be converted from comma-separated strings to arrays
array_fields = ['reserved', 'blockSizes', 'chromStarts', 'aliases']

# Fields that should be converted to integers
int_fields = ['chromStart', 'chromEnd', 'blockCount']

# Fields to remove from output
fields_to_remove = [
    'score', 'name', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames', 'type', 
    'annotationType', 'longName', 'syns', 'subCellLoc', 'pmids', 'thickStart', 'thickEnd'
]

# Field mappings for renaming
field_mappings = {
    'comments': 'repeatType',
    'geneName2': 'aliases'
}

# Read the TSV file
with open(input_file, 'r') as txt_file:
    # Skip header line if it exists
    first_line = txt_file.readline().strip()
    if (first_line.startswith('chrom')):
        headers = first_line.split('\t')
    else:
        # If no header, use the first line as data and define headers
        txt_file.seek(0)  # Reset file pointer to the beginning
        headers = [
            'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
            'thickStart', 'thickEnd', 'reserved', 'blockCount', 'blockSizes',
            'chromStarts', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames',
            'type', 'geneName', 'geneName2', 'geneType', 'status', 'annotationType',
            'position', 'longName', 'syns', 'subCellLoc', 'comments', 'uniProtId', 'pmids'
        ]
    
    reader = csv.DictReader(txt_file, fieldnames=headers, delimiter='\t')
    
    # Process each row
    for row in reader:
        # Create entry dictionary
        entry = {}
        
        for key, value in row.items():
            # Skip fields that should be removed
            if key in fields_to_remove:
                continue
                
            # Rename fields if they're in the mapping
            output_key = field_mappings.get(key, key)
            
            # Convert specified fields to arrays
            if key in array_fields and value:
                entry[output_key] = value.split(',') if ',' in value else [value.strip()]
            # Convert numeric fields to integers
            elif key in int_fields and value:
                try:
                    entry[output_key] = int(value)
                except ValueError:
                    entry[output_key] = value  # Keep as string if conversion fails
            else:
                entry[output_key] = value
        
        json_data.append(entry)

# Write the JSON data to the output file
with open(output_file, 'w') as json_file:
    json.dump(json_data, json_file, indent=2)

print(f"Data has been successfully converted to {output_file}")
print(f"Total entries: {len(json_data)}")