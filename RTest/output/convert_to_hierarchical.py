#!/usr/bin/env python3
import json
import argparse
import os
import sys
from collections import defaultdict
import traceback

def convert_to_hierarchical(input_file, output_file):
    """
    Convert the flat repeat JSON format to a hierarchical gene > protein > repeat structure
    """
    try:
        print(f"Converting {input_file} to hierarchical format...")
        
        # Check if input file exists
        if not os.path.exists(input_file):
            print(f"Error: Input file {input_file} does not exist.")
            return False
        
        # Load input data
        try:
            with open(input_file, 'r') as f:
                repeats_flat = json.load(f)
        except json.JSONDecodeError:
            print(f"Error: Could not parse {input_file} as JSON.")
            return False
        except Exception as e:
            print(f"Error reading input file: {e}")
            return False
            
        print(f"Loaded {len(repeats_flat)} repeats from input file.")
        
        # Initialize hierarchical structure
        hierarchical_data = {"genes": {}}
        
        # Process each repeat
        for repeat in repeats_flat:
            gene_name = repeat.get("geneName", "")
            if not gene_name:
                continue
                
            # Get/create gene entry
            if gene_name not in hierarchical_data["genes"]:
                hierarchical_data["genes"][gene_name] = {
                    "aliases": repeat.get("aliases", ""),
                    "geneType": repeat.get("geneType", ""),
                    "status": repeat.get("status", ""),
                    "proteins": {},
                    "transcripts": {}
                }
            
            # Extract protein
            protein_id = repeat.get("uniProtId", "")
            if not protein_id:
                continue
                
            # Get/create protein entry
            gene = hierarchical_data["genes"][gene_name]
            if protein_id not in gene["proteins"]:
                gene["proteins"][protein_id] = {"repeats": []}
            
            # Extract repeat info (leaving out transcript/exon info)
            repeat_entry = {
                "chrom": repeat.get("chrom", ""),
                "chromStart": repeat.get("chromStart", 0),
                "chromEnd": repeat.get("chromEnd", 0),
                "strand": repeat.get("strand", ""),
                "reserved": repeat.get("reserved", []),
                "blockCount": repeat.get("blockCount", 1),
                "blockSizes": repeat.get("blockSizes", ""),
                "chromStarts": repeat.get("chromStarts", ""),
                "position": repeat.get("position", "").replace(f" on protein {protein_id}", ""),
                "repeatType": repeat.get("repeatType", ""),
                "repeatLength": repeat.get("repeatLength", 0),
                "protein_start": repeat.get("protein_start", 0),
                "protein_end": repeat.get("protein_end", 0)
            }
            
            # Add repeat to protein
            gene["proteins"][protein_id]["repeats"].append(repeat_entry)
            
            # Process transcript information
            if "ensembl_exon_info" in repeat:
                exon_info = repeat.get("ensembl_exon_info", {})
                transcripts = exon_info.get("transcripts", [])
                
                for transcript in transcripts:
                    transcript_id = transcript.get("transcript_id", "")
                    if not transcript_id:
                        continue
                    
                    # Add transcript info if not already present
                    if transcript_id not in gene["transcripts"]:
                        gene["transcripts"][transcript_id] = {
                            "versioned_transcript_id": transcript.get("versioned_transcript_id", ""),
                            "is_canonical": transcript.get("is_canonical", False),
                            "biotype": transcript.get("biotype", ""),
                            "strand": transcript.get("strand", ""),
                            "exon_count": transcript.get("exon_count", 0),
                            "exons": []
                        }
                    
                    # Add exons for this transcript if not already present
                    containing_exons = transcript.get("containing_exons", [])
                    for exon in containing_exons:
                        # Check if exon is already in the list by ID
                        exon_id = exon.get("exon_id", "")
                        if exon_id:
                            # Check if this exon is already in the transcript
                            exon_exists = any(e.get("exon_id") == exon_id for e in gene["transcripts"][transcript_id]["exons"])
                            
                            if not exon_exists:
                                # Add exon without overlap information
                                exon_entry = {
                                    "exon_number": exon.get("exon_number", 0),
                                    "exon_id": exon_id,
                                    "exon_start": exon.get("exon_start", 0),
                                    "exon_end": exon.get("exon_end", 0),
                                    "coding_status": exon.get("coding_status", ""),
                                    "utr_status": exon.get("utr_status", ""),
                                    "coding_percentage": exon.get("coding_percentage", 0),
                                    "phase": exon.get("phase", 0),
                                    "end_phase": exon.get("end_phase", 0),
                                    "frame_status": exon.get("frame_status", "")
                                }
                                gene["transcripts"][transcript_id]["exons"].append(exon_entry)
        
        print(f"Processing complete. Found {len(hierarchical_data['genes'])} genes.")
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_file)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Write output file
        try:
            with open(output_file, 'w') as f:
                json.dump(hierarchical_data, f, indent=2)
            print(f"Conversion complete! Hierarchical format saved to {output_file}")
            return True
        except Exception as e:
            print(f"Error writing output file: {e}")
            return False
            
    except Exception as e:
        print(f"Unexpected error during conversion: {e}")
        traceback.print_exc()
        return False

def main():
    parser = argparse.ArgumentParser(
        description="Convert flat repeat JSON files to hierarchical format"
    )
    parser.add_argument("-i", "--input", required=True, help="Input JSON file or directory")
    parser.add_argument("-o", "--output", help="Output JSON file or directory")
    
    args = parser.parse_args()
    
    print(f"Running conversion with input: {args.input}, output: {args.output}")
    
    # Check if input is file or directory
    if os.path.isfile(args.input):
        # Single file conversion
        output_file = args.output or args.input.replace(".json", "_hierarchical.json")
        success = convert_to_hierarchical(args.input, output_file)
        return 0 if success else 1
    elif os.path.isdir(args.input):
        # Directory conversion
        output_dir = args.output or args.input
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Process each JSON file in directory
        for filename in os.listdir(args.input):
            if filename.endswith(".json"):
                input_path = os.path.join(args.input, filename)
                output_path = os.path.join(output_dir, filename.replace(".json", "_hierarchical.json"))
                convert_to_hierarchical(input_path, output_path)
    else:
        print(f"Error: Input path {args.input} does not exist.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
