import json
from collections import defaultdict
import os

def check_multiple_proteins(json_file_path):
    # Read the JSON data
    with open(json_file_path, 'r') as file:
        repeats_data = json.load(file)

    # Create a dictionary to map gene names to UniProt IDs
    gene_to_uniprot = defaultdict(set)

    # Populate the dictionary
    for repeat in repeats_data:
        # Skip entries that don't have both geneName and uniProtId fields
        if "geneName" in repeat and "uniProtId" in repeat:
            # Skip empty gene names
            if repeat["geneName"]:
                gene_to_uniprot[repeat["geneName"]].add(repeat["uniProtId"])

    # Find genes with multiple UniProt IDs
    multi_protein_genes = {gene: uniprot_ids for gene, uniprot_ids in gene_to_uniprot.items() 
                        if len(uniprot_ids) > 1}

    # Print the results
    print(f"Found {len(multi_protein_genes)} genes associated with multiple proteins:")
    for gene, uniprot_ids in multi_protein_genes.items():
        print(f"Gene '{gene}' is associated with {len(uniprot_ids)} proteins: {', '.join(uniprot_ids)}")

    return multi_protein_genes

if __name__ == "__main__":
    # Path to the data file
    file_path = r"c:\Users\ojfab\Documents\GitHub\Tandem-Repeat-Domain-Database\merge_test_small\data\gname_hg38_repeats.json"
    
    if not os.path.exists(file_path):
        print(f"Error: File not found at {file_path}")
    else:
        multi_protein_genes = check_multiple_proteins(file_path)
        
        # Save results to a file
        output_dir = r"c:\Users\ojfab\Documents\GitHub\Tandem-Repeat-Domain-Database\merge_test_small\data"
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, "multi_protein_genes.txt")
        
        with open(output_file, 'w') as file:
            file.write(f"Found {len(multi_protein_genes)} genes associated with multiple proteins:\n")
            for gene, uniprot_ids in multi_protein_genes.items():
                file.write(f"Gene '{gene}' is associated with {len(uniprot_ids)} proteins: {', '.join(uniprot_ids)}\n")
        
        print(f"Results saved to {output_file}")
