import json
from collections import defaultdict

# --- Step 1: Build alias groups and collect metadata from the gene alias JSON ---
with open('data/uniprot_gene_alias.json', 'r') as f:
    alias_entries = json.load(f)

# Initialize tracking metrics
total_alias_entries = len(alias_entries)
total_aliases_found = 0
aliases_without_gene_names = 0

parent = {}
# Store metadata for each gene
gene_metadata = {}

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
            aliases_without_gene_names += 1
            continue
        
        total_aliases_found += len(aliases)
        
        # Store metadata for each gene name
        metadata = {
            "uniprot_id": entry.get("Entry", ""),
            "protein_name": entry.get("Protein names", ""),
            "entry_name": entry.get("Entry Name", ""),
            "reviewed": entry.get("Reviewed", "")
        }
        
        for alias in aliases:
            if alias not in parent:
                parent[alias] = alias
            
            # Store metadata for this gene
            gene_metadata[alias] = metadata
        
        first = aliases[0]
        for alias in aliases[1:]:
            union(first, alias)

# Build a mapping from the representative alias to the complete set of aliases.
alias_groups = {}
for alias in parent:
    rep = find(alias)
    alias_groups.setdefault(rep, set()).add(alias)

# Calculate statistics about alias groups
total_unique_genes = len(alias_groups)
genes_with_aliases = sum(1 for aliases in alias_groups.values() if len(aliases) > 1)
aliases_per_gene = [len(aliases) for aliases in alias_groups.values()]
max_aliases = max(aliases_per_gene) if aliases_per_gene else 0
avg_aliases = sum(aliases_per_gene) / len(aliases_per_gene) if aliases_per_gene else 0

# --- Step 2: Load transcripts file and merge records based on alias groups ---
with open('data/ensembl_transcripts.json', 'r') as f:
    transcript_records = json.load(f)

total_transcript_records = len(transcript_records)
merged = {}  # Key: representative alias; Value: merged record
merged_count = 0  # Counter for when we merge transcripts

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
        original_transcript_count = len(merged[group_key]["transcripts"])
        merged[group_key]["transcripts"].update(transcripts)
        new_transcript_count = len(merged[group_key]["transcripts"])
        
        if new_transcript_count > original_transcript_count:
            merged_count += 1
            
        merged[group_key]["aliases"].update(alias_groups.get(group_key, {gene}))
    else:
        # Get the canonical name for this group
        canonical_name = min(alias_groups.get(group_key, {gene}))
        
        # Get metadata for canonical name if available
        metadata = gene_metadata.get(canonical_name, {})
        
        merged[group_key] = {
            "gene": canonical_name,
            "ensembl_id": ensembl_id,
            "uniprot_id": metadata.get("uniprot_id", ""),
            "protein_name": metadata.get("protein_name", ""),
            "entry_name": metadata.get("entry_name", ""),
            "reviewed": metadata.get("reviewed", ""),
            "aliases": set(alias_groups.get(group_key, {gene})),
            "transcripts": transcripts
        }

# --- Step 3: Add any alias groups from uniprot_gene_alias.json that did not have transcript entries ---
print("Adding genes with no transcript records...")
added_count = 0

# Add entries for gene groups that have no transcript records
for rep, aliases in alias_groups.items():
    if rep not in merged:
        canonical_name = min(aliases)
        metadata = gene_metadata.get(canonical_name, {})
        
        merged[rep] = {
            "gene": canonical_name,
            "ensembl_id": None,
            "uniprot_id": metadata.get("uniprot_id", ""),
            "protein_name": metadata.get("protein_name", ""),
            "entry_name": metadata.get("entry_name", ""),
            "reviewed": metadata.get("reviewed", ""),
            "aliases": aliases,
            "transcripts": set()
        }
        added_count += 1

print(f"Added {added_count} genes that had aliases but no transcript records")

# Convert sets to sorted lists for final JSON output.
merged_list = []
for rec in merged.values():
    merged_list.append({
        "gene": rec["gene"],
        "ensembl_id": rec["ensembl_id"],
        "uniprot_id": rec["uniprot_id"],
        "protein_name": rec["protein_name"],
        "entry_name": rec["entry_name"],
        "reviewed": rec["reviewed"],
        "aliases": sorted(list(rec["aliases"])),
        "transcripts": sorted(list(rec["transcripts"]))
    })

# Calculate transcript statistics
total_transcripts_before = sum(len(record.get("transcripts", [])) for record in transcript_records)
total_transcripts_after = sum(len(record["transcripts"]) for record in merged_list)
genes_with_transcripts = sum(1 for record in merged_list if record["transcripts"])
genes_without_transcripts = len(merged_list) - genes_with_transcripts

# --- Step 4: Write merged records to a new JSON file ---
with open('data/ensembl_transcripts_merged.json', 'w') as outfile:
    json.dump(merged_list, outfile, indent=4)

# --- Step 5: Print a summary report ---
print("\n" + "="*50)
print("MERGE REPORT")
print("="*50)
print(f"Input Records:")
print(f"- Total alias entries: {total_alias_entries}")
print(f"- Total unique aliases found: {total_aliases_found}")
print(f"- Entries without gene names: {aliases_without_gene_names}")
print(f"- Total transcript records: {total_transcript_records}")
print("\nAlias Groups:")
print(f"- Unique gene groups: {total_unique_genes}")
print(f"- Genes with multiple aliases: {genes_with_aliases} ({genes_with_aliases/total_unique_genes*100:.1f}%)")
print(f"- Maximum aliases per gene: {max_aliases}")
print(f"- Average aliases per gene: {avg_aliases:.2f}")
print("\nMerge Results:")
print(f"- Records after merging: {len(merged_list)} (includes {added_count} genes without transcript data)")
print(f"- Transcript merges performed: {merged_count}")
print(f"- Genes with transcripts: {genes_with_transcripts} ({genes_with_transcripts/len(merged_list)*100:.1f}%)")
print(f"- Genes without transcripts: {genes_without_transcripts} ({genes_without_transcripts/len(merged_list)*100:.1f}%)")
print(f"- Total transcripts before: {total_transcripts_before}")
print(f"- Total transcripts after: {total_transcripts_after}")
print(f"- Duplicate transcripts removed: {total_transcripts_before-total_transcripts_after}")
print("="*50)

print("\nMerged transcripts saved to data/ensembl_transcripts_merged.json")