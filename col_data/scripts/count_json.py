#!/usr/bin/env python3
# filepath: /home/dogdorgesh/Documents/Github/Tandem-Repeat-Domain-Database/merge_test_small/data/count_json.py
import json
import sys
from collections import defaultdict

def analyze_json_file(file_path):
    """Analyze JSON array file for various statistics"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
            if isinstance(data, list):
                total_entries = len(data)
                
                # Count entries without geneName
                entries_without_genename = sum(1 for item in data if "geneName" not in item)
                
                # Count entries with blockCount = 1
                entries_with_blockcount_1 = sum(1 for item in data 
                                               if "blockCount" in item and item["blockCount"] == 1)
                
                # Count entries with missing blockCount field
                entries_without_blockcount = sum(1 for item in data if "blockCount" not in item)
                
                # Count by gene and repeat type to find those with ≥5 repeats of same type
                # BUT ONLY COUNT REPEATS WITH BLOCKCOUNT = 1
                gene_repeat_counts = defaultdict(lambda: defaultdict(int))
                
                # First pass: count repeats by gene and type (only blockCount=1)
                for item in data:
                    if ("geneName" in item and "repeatType" in item and 
                        "blockCount" in item and item["blockCount"] == 1):
                        gene_name = item["geneName"]
                        repeat_type = item["repeatType"]
                        gene_repeat_counts[gene_name][repeat_type] += 1
                
                # Find genes with ≥5 repeats of the same type (with blockCount=1)
                genes_with_5plus_repeats = set()
                for gene, type_counts in gene_repeat_counts.items():
                    for repeat_type, count in type_counts.items():
                        if count >= 5:
                            genes_with_5plus_repeats.add((gene, repeat_type))
                
                # Second pass: count entries satisfying all criteria
                entries_satisfying_all_criteria = 0
                entries_satisfying_all_with_canonical = 0
                entries_with_blockcount_1_and_inframe = 0
                
                for item in data:
                    if "blockCount" in item and item["blockCount"] == 1:
                        # Check if this entry has any in-frame exons
                        has_inframe_exon = False
                        has_inframe_canonical_exon = False
                        
                        # Ensembl exon info is stored in a nested structure
                        if "ensembl_exon_info" in item and isinstance(item["ensembl_exon_info"], dict):
                            exon_info = item["ensembl_exon_info"]
                             
                            # Check transcripts for in-frame exons
                            if "transcripts" in exon_info and isinstance(exon_info["transcripts"], list):
                                for transcript in exon_info["transcripts"]:
                                    is_canonical = False
                                    if "is_canonical" in transcript and transcript["is_canonical"]:
                                        is_canonical = True
                                    
                                    if "containing_exons" in transcript and isinstance(transcript["containing_exons"], list):
                                        for exon in transcript["containing_exons"]:
                                            if "frame_status" in exon and exon["frame_status"] == "in_frame":
                                                has_inframe_exon = True
                                                if is_canonical:
                                                    has_inframe_canonical_exon = True
                                                    break
                                        if has_inframe_canonical_exon:
                                            break
                                    if has_inframe_exon and not has_inframe_canonical_exon:
                                        continue  # Keep checking other transcripts for canonical ones
                        
                        if has_inframe_exon:
                            entries_with_blockcount_1_and_inframe += 1
                            
                            # Check if this is part of a gene with ≥5 repeats of the same type
                            if ("geneName" in item and "repeatType" in item and 
                                (item["geneName"], item["repeatType"]) in genes_with_5plus_repeats):
                                entries_satisfying_all_criteria += 1
                                
                                # Also check if it has a canonical transcript with in-frame exons
                                if has_inframe_canonical_exon:
                                    entries_satisfying_all_with_canonical += 1
                
                # Get the number of unique genes with repeats satisfying all criteria
                genes_satisfying_all = set()
                genes_satisfying_all_with_canonical = set()
                
                for item in data:
                    if ("blockCount" in item and item["blockCount"] == 1 and
                        "geneName" in item and "repeatType" in item and
                        (item["geneName"], item["repeatType"]) in genes_with_5plus_repeats):
                        
                        # Check for in-frame exons
                        has_inframe_exon = False
                        has_inframe_canonical_exon = False
                        
                        if "ensembl_exon_info" in item and isinstance(item["ensembl_exon_info"], dict):
                            exon_info = item["ensembl_exon_info"]
                            if "transcripts" in exon_info and isinstance(exon_info["transcripts"], list):
                                for transcript in exon_info["transcripts"]:
                                    is_canonical = False
                                    if "is_canonical" in transcript and transcript["is_canonical"]:
                                        is_canonical = True
                                        
                                    if "containing_exons" in transcript and isinstance(transcript["containing_exons"], list):
                                        for exon in transcript["containing_exons"]:
                                            if "frame_status" in exon and exon["frame_status"] == "in_frame":
                                                has_inframe_exon = True
                                                if is_canonical:
                                                    has_inframe_canonical_exon = True
                                                    break
                                        if has_inframe_canonical_exon:
                                            break
                        
                        if has_inframe_exon:
                            genes_satisfying_all.add(item["geneName"])
                            
                        if has_inframe_canonical_exon:
                            genes_satisfying_all_with_canonical.add(item["geneName"])
                
                # Add the actual gene names to the results
                results = {
                    "total_entries": total_entries,
                    "entries_without_genename": entries_without_genename,
                    "entries_with_genename": total_entries - entries_without_genename,
                    "entries_with_blockcount_1": entries_with_blockcount_1,
                    "entries_without_blockcount": entries_without_blockcount,
                    "entries_with_blockcount_1_and_inframe": entries_with_blockcount_1_and_inframe,
                    "entries_satisfying_all_criteria": entries_satisfying_all_criteria,
                    "entries_satisfying_all_with_canonical": entries_satisfying_all_with_canonical,
                    "genes_with_5plus_repeats_count": len(genes_with_5plus_repeats),
                    "genes_satisfying_all_criteria_count": len(genes_satisfying_all),
                    "genes_satisfying_all_with_canonical_count": len(genes_satisfying_all_with_canonical),
                    "genes_satisfying_all": sorted(list(genes_satisfying_all)),  # Add the actual gene names
                    "genes_satisfying_all_with_canonical": sorted(list(genes_satisfying_all_with_canonical))  # Add canonical gene names
                }
                
                return results
            else:
                return "Error: JSON root is not a list"
    except json.JSONDecodeError:
        return "Error: Invalid JSON format"
    except Exception as e:
        return f"Error: {str(e)}"

if __name__ == "__main__":
    file_path = sys.argv[1] if len(sys.argv) > 1 else "output/1000_test_exons_hg38_repeats.json"
    results = analyze_json_file(file_path)
    
    if isinstance(results, dict):
        print(f"Analysis of {file_path}:")
        print(f"  Total entries: {results['total_entries']}")
        print(f"  Entries without geneName: {results['entries_without_genename']} ({results['entries_without_genename']/results['total_entries']*100:.2f}%)")
        print(f"  Entries with geneName: {results['entries_with_genename']} ({results['entries_with_genename']/results['total_entries']*100:.2f}%)")
        print(f"  Entries with blockCount = 1: {results['entries_with_blockcount_1']} ({results['entries_with_blockcount_1']/results['total_entries']*100:.2f}%)")
        print(f"  Entries without blockCount field: {results['entries_without_blockcount']} ({results['entries_without_blockcount']/results['total_entries']*100:.2f}%)")
        print(f"  Entries with blockCount = 1 AND in-frame exons: {results['entries_with_blockcount_1_and_inframe']} ({results['entries_with_blockcount_1_and_inframe']/results['total_entries']*100:.2f}%)")
        print(f"  Number of gene-repeatType combinations with ≥5 repeats of blockCount=1: {results['genes_with_5plus_repeats_count']}")
        print(f"  Entries satisfying all criteria: {results['entries_satisfying_all_criteria']} ({results['entries_satisfying_all_criteria']/results['total_entries']*100:.2f}%)")
        print(f"  Entries satisfying all criteria AND in canonical transcripts: {results['entries_satisfying_all_with_canonical']} ({results['entries_satisfying_all_with_canonical']/results['total_entries']*100:.2f}%)")
        print(f"  Number of unique genes satisfying all criteria: {results['genes_satisfying_all_criteria_count']}")
        print(f"  Number of unique genes satisfying all criteria AND in canonical transcripts: {results['genes_satisfying_all_with_canonical_count']}")
        
        # Print the actual gene names
        print(f"\nGenes satisfying all criteria ({results['genes_satisfying_all_criteria_count']}):")
        if results['genes_satisfying_all_criteria_count'] > 0:
            # Print all genes or limit to 20 with indication if there are more
            gene_list = results["genes_satisfying_all"]
            if len(gene_list) <= 20:
                print(", ".join(gene_list))
            else:
                print(", ".join(gene_list[:20]) + f", ... and {len(gene_list) - 20} more genes")
                print(f"(Use '--save-genes filename.txt' option to save the complete list)")
        
        # Print the canonical genes
        print(f"\nGenes satisfying all criteria AND in canonical transcripts ({results['genes_satisfying_all_with_canonical_count']}):")
        if results['genes_satisfying_all_with_canonical_count'] > 0:
            gene_list = results["genes_satisfying_all_with_canonical"]
            if len(gene_list) <= 20:
                print(", ".join(gene_list))
            else:
                print(", ".join(gene_list[:20]) + f", ... and {len(gene_list) - 20} more genes")
                print(f"(Use '--save-genes filename.txt' option to save the complete list)")
    else:
        print(results)