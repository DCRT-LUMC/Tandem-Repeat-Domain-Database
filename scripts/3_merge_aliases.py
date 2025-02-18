import json

# --- Step 1: Build alias groups from the gene alias JSON ---
with open('data/uniprot_gene_alias.json', 'r') as f:
    alias_entries = json.load(f)

parent = {}

def find(x):
    if parent[x] != x:
        parent[x] = find(parent[x])
    return parent[x]

def union(x, y):
    root_x = find(x)
    root_y = find(y)
    if root_x != root_y:
        parent[root_y] = root_x

# For each alias list, initialize and then union all aliases.
for entry in alias_entries:
    if isinstance(entry, dict) and "Gene Names" in entry:
        aliases = [a.strip() for a in entry["Gene Names"] if a.strip()]
        if not aliases:
            continue
        for alias in aliases:
            if alias not in parent:
                parent[alias] = alias
        first = aliases[0]
        for alias in aliases[1:]:
            union(first, alias)

# Build a mapping from the representative alias to the complete set of aliases.
alias_groups = {}
for alias in parent:
    rep = find(alias)
    alias_groups.setdefault(rep, set()).add(alias)

# --- Step 2: Load transcripts file and merge records based on alias groups ---
with open('data/ensembl_transcripts_test.json', 'r') as f:
    transcript_records = json.load(f)

merged = {}  # Key: representative alias; Value: merged record

for record in transcript_records:
    gene = record.get("gene")
    ensembl_id = record.get("ensembl_id")
    transcripts = set(record.get("transcripts", []))
    
    # Determine the group key: if the gene is part of an alias group, use its representative.
    if gene in parent:
        group_key = find(gene)
    else:
        group_key = gene

    if group_key in merged:
        # Merge transcripts and add any new aliases.
        merged[group_key]["transcripts"].update(transcripts)
        merged[group_key]["aliases"].update(alias_groups.get(group_key, {gene}))
    else:
        merged[group_key] = {
            "gene": min(alias_groups.get(group_key, {gene})),  # Define a canonical name
            "ensembl_id": ensembl_id,
            "aliases": set(alias_groups.get(group_key, {gene})),
            "transcripts": transcripts
        }

# Convert sets to sorted lists for final JSON output.
merged_list = []
for rec in merged.values():
    merged_list.append({
        "gene": rec["gene"],
        "ensembl_id": rec["ensembl_id"],
        "aliases": sorted(list(rec["aliases"])),
        "transcripts": sorted(list(rec["transcripts"]))
    })

# --- Step 3: Write merged records to a new JSON file ---
with open('data/ensembl_transcripts_merged.json', 'w') as outfile:
    json.dump(merged_list, outfile, indent=4)

print("Merged transcripts saved to data/ensembl_transcripts_merged.json")