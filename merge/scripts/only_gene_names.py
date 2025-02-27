import json

#export all gene names and aliases

def extract_gene_names(json_file, output_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    gene_names = set()
    for item in data:
        if isinstance(item, dict) and "Gene Names" in item:
            for name in item["Gene Names"]:
                if name.strip():
                    gene_names.add(name.strip())
    with open(output_file, 'w') as out:
        for name in sorted(gene_names):
            out.write(name + "\n")

if __name__ == "__main__":
    extract_gene_names('data/uniprot_gene_alias.json', 'data/all_gene_names.txt')
    print("Gene names have been written to all_gene_names.txt")