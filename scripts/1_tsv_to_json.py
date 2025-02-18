import csv
import json
import re

# This converts the .tsv export from UniProt into a usable .json

input_file = 'data/uniprot_gene_alias.tsv'
output_file = 'data/uniprot_gene_alias.json'

json_data = []

# Read the TSV file
with open(input_file, 'r') as tsv_file:
    reader = csv.DictReader(tsv_file, delimiter='\t')
    for row in reader:
        # Split the gene names by spaces or slashes and store them as a list of aliases
        gene_names = re.split(r'[ /]+', row['Gene Names'])
        entry = {
            'From': row['From'],
            'Entry': row['Entry'],
            'Reviewed': row['Reviewed'],
            'Entry Name': row['Entry Name'],
            'Protein names': row['Protein names'],
            'Gene Names': gene_names,
            'Length': row['Length']
        }
        json_data.append(entry)

# Write the JSON data to the output file
with open(output_file, 'w') as json_file:
    json.dump(json_data, json_file, indent=4)

print(f"Data has been successfully converted to {output_file}")