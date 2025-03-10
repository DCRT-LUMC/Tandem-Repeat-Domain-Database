#!/usr/bin/env python3

import json
import requests
import argparse
import concurrent.futures
from tqdm import tqdm
import mysql.connector

def get_gene_info(uniprot_id, chrom, start, end):
    """Get gene information from UniProt and/or UCSC"""
    gene_info = {"geneName": "", "geneName2": []}
    
    # Try UniProt first
    if uniprot_id:
        try:
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
            response = requests.get(url, headers={"Accept": "application/json"})
            if response.status_code == 200:
                data = response.json()
                gene_names = []
                gene_synonyms = []
                
                # Extract gene names
                for gene in data.get("genes", []):
                    if "geneName" in gene and "value" in gene["geneName"]:
                        gene_names.append(gene["geneName"]["value"])
                    for syn in gene.get("synonyms", []):
                        if "value" in syn:
                            gene_synonyms.append(syn["value"])
                
                if gene_names:
                    gene_info["geneName"] = ",".join(gene_names)
                    gene_info["geneName2"] = gene_synonyms  # Store as array
                    return gene_info
        except Exception:
            pass
    
    # Fall back to UCSC coordinates
    if chrom and start and end:
        try:
            conn = mysql.connector.connect(
                host="genome-mysql.soe.ucsc.edu", 
                user="genomep", 
                password="genome", 
                database="hg38"
            )
            cursor = conn.cursor(dictionary=True)
            cursor.execute(
                """SELECT DISTINCT r.name2 as gene_name, kgx.geneSymbol, kgx.alias 
                FROM refGene r 
                LEFT JOIN knownGene kg ON r.name = kg.name 
                LEFT JOIN kgXref kgx ON kg.name = kgx.kgID 
                WHERE r.chrom = %s AND r.txStart <= %s AND r.txEnd >= %s""",
                (chrom, end, start)
            )
            results = cursor.fetchall()
            cursor.close()
            conn.close()
            
            gene_names = set(r["gene_name"] for r in results if r["gene_name"])
            gene_aliases = set()
            for r in results:
                if r.get("alias"):
                    gene_aliases.update(a.strip() for a in r["alias"].split(" "))
            
            if gene_names:
                gene_info["geneName"] = ",".join(gene_names)
                gene_info["geneName2"] = list(gene_aliases)  # Store as array
        except Exception:
            pass
    
    return gene_info

def annotate_repeats(input_file, output_file, threads=8):
    """Process all repeats in JSON file to add gene info"""
    with open(input_file, 'r') as f:
        repeats = json.load(f)
    
    # Process in parallel with progress bar
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for repeat in repeats:
            futures.append(
                executor.submit(
                    get_gene_info,
                    repeat.get("uniProtId"),
                    repeat.get("chrom"),
                    repeat.get("chromStart"),
                    repeat.get("chromEnd")
                )
            )
        
        # Update repeats with results
        for i, future in enumerate(tqdm(concurrent.futures.as_completed(futures), 
                                        total=len(futures), 
                                        desc="Annotating genes")):
            repeats[i].update(future.result())
    
    # Write updated data
    with open(output_file, 'w') as f:
        json.dump(repeats, f, indent=2)
    
    print(f"Processed {len(repeats)} records. Output saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add gene information to repeat JSON file")
    parser.add_argument("input", help="Input JSON file with repeat records")
    parser.add_argument("-o", "--output", help="Output JSON file (default: input_annotated.json)")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of parallel threads")
    
    args = parser.parse_args()
    output = args.output or args.input.replace(".json", "_annotated.json")
    
    annotate_repeats(args.input, output, args.threads)