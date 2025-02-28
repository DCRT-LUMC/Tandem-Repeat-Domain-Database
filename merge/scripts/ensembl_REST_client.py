#!/usr/bin/env python

import sys
import json
import time
import requests
from tqdm import tqdm

class EnsemblRestClient:
    """Client for the Ensembl REST API with rate limiting and error handling.
    Based on the example from the Ensembl REST API documentation.
    custom input or output paths run with python ensembl_REST_client.py data/custom_input.json data/custom_output.json"""
    
    def __init__(self, server='https://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        """Execute a REST request with rate limiting and error handling."""
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        url = self.server + endpoint
        data = None

        # Check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        
        try:
            response = requests.get(url, headers=hdrs, params=params)
            response.raise_for_status()
            if response.content:
                data = response.json()
            self.req_count += 1

        except requests.exceptions.HTTPError as e:
            # Check if we are being rate limited by the server
            if e.response.status_code == 429:
                if 'Retry-After' in e.response.headers:
                    retry = e.response.headers['Retry-After']
                    time.sleep(float(retry))
                    return self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(f'Request failed for {endpoint}: Status code: {e.response.status_code} Reason: {e.response.reason}\n')
           
        return data

    def get_ensembl_id(self, gene_name, species='homo_sapiens'):
        """Get Ensembl ID for a gene symbol."""
        data = self.perform_rest_action(
            endpoint=f'/xrefs/symbol/{species}/{gene_name}',
            params={'object_type': 'gene'}
        )
        if data and len(data) > 0:
            return data[0]['id']
        return None

    def get_transcripts(self, ensembl_id):
        """Get all transcript IDs for an Ensembl gene ID."""
        data = self.perform_rest_action(
            endpoint=f'/lookup/id/{ensembl_id}',
            params={'expand': '1'}
        )
        if data and 'Transcript' in data:
            return [transcript["id"] for transcript in data["Transcript"]]
        return []


def load_gene_names(filename):
    """Load and extract unique gene names from the source JSON file."""
    with open(filename, 'r') as f:
        entries = json.load(f)
    
    gene_set = set()
    for entry in entries:
        if isinstance(entry, dict) and "Gene Names" in entry:
            for gene in entry["Gene Names"]:
                if gene.strip():
                    gene_set.add(gene.strip())
    
    return sorted(list(gene_set))


def process_genes(gene_names, client):
    """Process each gene to retrieve Ensembl IDs and transcripts."""
    ensembl_db = []
    
    # Use tqdm to create a progress bar
    for gene in tqdm(gene_names, desc="Processing genes", unit="gene"):
        ensembl_id = client.get_ensembl_id(gene)
        if ensembl_id:
            transcripts = client.get_transcripts(ensembl_id)
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
    
    return ensembl_db


def write_output(ensembl_db, filename):
    """Write results to a JSON file."""
    with open(filename, 'w') as outfile:
        json.dump(ensembl_db, outfile, indent=4)


def run(input_file='data/uniprot_gene_alias.json', output_file='data/ensembl_transcripts.json'):
    """Main function to run the Ensembl data extraction process."""
    # Initialize the REST client
    client = EnsemblRestClient(reqs_per_sec=15)
    
    # Load and process the gene names
    gene_names = load_gene_names(input_file)
    total_genes = len(gene_names)
    print(f"Found {total_genes} unique gene names to process")
    
    # Process genes
    ensembl_db = process_genes(gene_names, client)
    
    # Save results
    write_output(ensembl_db, output_file)
    
    print(f"\nProcessed {total_genes} genes. Found Ensembl IDs for {len(ensembl_db)} genes.")
    print(f"New Ensembl database created as {output_file}")


if __name__ == '__main__':
    if len(sys.argv) > 2:
        input_file, output_file = sys.argv[1:3]
        run(input_file, output_file)
    else:
        run()  # Use default file paths