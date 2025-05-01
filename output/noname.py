import json
import os

# Use absolute path or correct relative path
file_path = '/home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/output/DEF_gname_hg38_repeats.json'

# Load the JSON data
with open(file_path, 'r') as f:
    data = json.load(f)

# Initialize counters
missing_gene_name_count = 0
missing_aliases_count = 0
unreviewed_count = 0
missing_gene_and_unreviewed = 0
total_entries = len(data)

# Track unique proteins without gene names
proteins_without_gene_names = set()
unreviewed_proteins = set()
unreviewed_without_gene_names = set()

# Check each entry
for entry in data:
    # Check for missing or empty geneName
    has_no_gene_name = 'geneName' not in entry or entry.get('geneName', '') == ""
    is_unreviewed = 'status' in entry and 'Unreviewed' in entry['status']
    
    if has_no_gene_name:
        missing_gene_name_count += 1
        
        # Add the UniProtId to our set if it exists
        if 'uniProtId' in entry and entry['uniProtId']:
            proteins_without_gene_names.add(entry['uniProtId'])
    
    if is_unreviewed:
        unreviewed_count += 1
        if 'uniProtId' in entry and entry['uniProtId']:
            unreviewed_proteins.add(entry['uniProtId'])
            
    if has_no_gene_name and is_unreviewed:
        missing_gene_and_unreviewed += 1
        if 'uniProtId' in entry and entry['uniProtId']:
            unreviewed_without_gene_names.add(entry['uniProtId'])
    
    # Check for missing or empty aliases - handle different formats
    if ('aliases' not in entry or 
        entry.get('aliases', '') == "" or 
        entry.get('aliases', []) == [] or
        (isinstance(entry.get('aliases', []), list) and len(entry.get('aliases', [])) == 0)):
        missing_aliases_count += 1

# Print results
print(f"Total entries: {total_entries}")
print(f"Entries without geneName: {missing_gene_name_count} ({missing_gene_name_count/total_entries*100:.2f}%)")
print(f"Unique proteins without geneName: {len(proteins_without_gene_names)}")
print(f"Entries with unreviewed status: {unreviewed_count} ({unreviewed_count/total_entries*100:.2f}%)")
print(f"Entries both without geneName and unreviewed: {missing_gene_and_unreviewed} ({missing_gene_and_unreviewed/total_entries*100:.2f}%)")
print(f"Unique unreviewed proteins without gene names: {len(unreviewed_without_gene_names)}")
print(f"Entries without aliases: {missing_aliases_count} ({missing_aliases_count/total_entries*100:.2f}%)")

# Save proteins without gene names to a text file
output_file = '/home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/output/proteins_without_gene_names.txt'
with open(output_file, 'w') as f:
    for protein in sorted(proteins_without_gene_names):
        f.write(f"{protein}\n")

# Save unreviewed proteins list
unreviewed_file = '/home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/output/unreviewed_proteins.txt'
with open(unreviewed_file, 'w') as f:
    for protein in sorted(unreviewed_proteins):
        f.write(f"{protein}\n")

# Save unreviewed proteins without gene names
unreviewed_no_gene_file = '/home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/output/unreviewed_proteins_without_gene_names.txt'
with open(unreviewed_no_gene_file, 'w') as f:
    for protein in sorted(unreviewed_without_gene_names):
        f.write(f"{protein}\n")

print(f"\nAll proteins without gene names have been saved to: {output_file}")
print(f"All unreviewed proteins have been saved to: {unreviewed_file}")
print(f"Unreviewed proteins without gene names have been saved to: {unreviewed_no_gene_file}")

# Calculate percentage of proteins without gene names that are unreviewed
if len(proteins_without_gene_names) > 0:
    percentage = len(unreviewed_without_gene_names) / len(proteins_without_gene_names) * 100
    print(f"\n{percentage:.2f}% of proteins without gene names are unreviewed (TrEMBL entries)")