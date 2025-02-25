import requests
import json
from tqdm import tqdm
import time

''' Reads the gene names from the uniprot_gene_alias.json,
    looks up their Ensembl IDs and transcripts,
    and then writes a new JSON file with the results'''

def get_ensembl_id(gene_name, species='homo_sapiens'):
    url = f"https://rest.ensembl.org/xrefs/symbol/{species}/{gene_name}?object_type=gene"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        if data:
            return data[0]['id']
    return None

def get_transcripts(ensembl_id):
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?expand=1"
    headers = {"Content-Type": "application/json"}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        data = response.json()
        transcripts = data.get("Transcript", [])
        return [transcript["id"] for transcript in transcripts]
    else:
        return []

# Read the source JSON file and extract unique gene names
with open('data/uniprot_gene_alias.json', 'r') as f:
    entries = json.load(f)

gene_set = set()
for entry in entries:
    if isinstance(entry, dict) and "Gene Names" in entry:
        for gene in entry["Gene Names"]:
            if gene.strip():
                gene_set.add(gene.strip())

gene_names = sorted(list(gene_set))
total_genes = len(gene_names)

print(f"Found {total_genes} unique gene names to process")

# Create a list to hold the new database records
ensembl_db = []

# Use tqdm to create a progress bar
for gene in tqdm(gene_names, desc="Processing genes", unit="gene"):
    # Small delay to prevent hitting API rate limits
    time.sleep(0.05)
    
    ensembl_id = get_ensembl_id(gene)
    if ensembl_id:
        transcripts = get_transcripts(ensembl_id)
        record = {
            "gene": gene,
            "ensembl_id": ensembl_id,
            "transcripts": transcripts
        }
        ensembl_db.append(record)
        
        # Update the progress bar description with the current gene
        tqdm.write(f"✓ Gene: {gene} (Ensembl ID: {ensembl_id}), Transcripts: {len(transcripts)}")
    else:
        tqdm.write(f"✗ Gene: {gene} not found in Ensembl")

# Write the new JSON database file
with open('data/ensembl_transcripts.json', 'w') as outfile:
    json.dump(ensembl_db, outfile, indent=4)

print(f"\nProcessed {total_genes} genes. Found Ensembl IDs for {len(ensembl_db)} genes.")
print("New Ensembl database created as data/ensembl_transcripts.json")