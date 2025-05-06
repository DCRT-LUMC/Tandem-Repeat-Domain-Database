import json
from collections import defaultdict
import sys
import re

def get_transcript_obj(repeat):
    transcripts = repeat.get("ensembl_exon_info", {}).get("transcripts", [])
    for t in transcripts:
        if t.get("is_canonical"):
            return t
    if transcripts:
        return transcripts[0]
    return None

def get_transcript_id(transcript):
    if not transcript:
        return "unknown_transcript"
    return transcript.get("transcript_id") or transcript.get("transcript_name") or "unknown_transcript"

def get_gene_name(repeat):
    return repeat.get("geneName") or "unknown_gene"

def get_gene_metadata(repeat):
    return {
        "aliases": repeat.get("aliases"),
        "geneType": repeat.get("geneType"),
    }

def get_protein_id(repeat):
    return repeat.get("uniProtId") or "unknown_protein"

def get_protein_metadata(repeat):
    return {
        "uniProtId": repeat.get("uniProtId"),
        # Remove repeatType from protein-level metadata
        "status": repeat.get("status"),
    }

def get_chrom(repeat):
    return repeat.get("chrom") or "unknown_chrom"

def get_repeat_specific_info(repeat):
    keys_to_remove = {
        "geneName", "uniProtId", "repeatType", "status", "aliases", "geneType",
        "chrom", "ensembl_exon_info"
    }
    return {k: v for k, v in repeat.items() if k not in keys_to_remove}

def exon_key(exon):
    # Use exon_id as unique key
    return exon.get("exon_id")

def normalize_repeat_type(repeat_type):
    if not repeat_type:
        return "Unknown"
    repeat_type = repeat_type.strip()
    if repeat_type.startswith("WD"):
        return "WD"
    if repeat_type.startswith("Alpha"):
        return "Alpha"
    if repeat_type.startswith("PUR"):
        return "PUR"
    if repeat_type.startswith("RCC1"):
        return "RCC1"
    if repeat_type.startswith("SVP"):
        return "SVP"
    if repeat_type.startswith("Spectrin"):
        return "Spectrin"
    # If starts with a number
    if repeat_type and repeat_type[0].isdigit():
        return "Unnamed"
    # If the entire string is a roman numeral I-X (case-insensitive)
    if re.match(r"^(I|II|III|IV|V|VI|VII|VIII|IX|X)$", repeat_type, re.IGNORECASE):
        return "Unnamed"
    return repeat_type

def main(input_json, output_json):
    with open(input_json) as f:
        flat = json.load(f)

    hierarchical = defaultdict(
        lambda: defaultdict(
            lambda: {
                "gene_metadata": None,
                "transcripts": defaultdict(
                    lambda: {
                        "transcript": None,
                        "proteins": defaultdict(lambda: {
                            "protein_metadata": None,
                            "repeat_types": defaultdict(lambda: {"repeats": [], "exons": []})
                        })
                    }
                )
            }
        )
    )

    for repeat in flat:
        chrom = get_chrom(repeat)
        gene = get_gene_name(repeat)
        transcript_obj = get_transcript_obj(repeat)
        transcript_id = get_transcript_id(transcript_obj)
        protein_id = get_protein_id(repeat)

        gene_slot = hierarchical[chrom][gene]
        if gene_slot["gene_metadata"] is None:
            gene_slot["gene_metadata"] = get_gene_metadata(repeat)

        transcript_slot = gene_slot["transcripts"][transcript_id]
        if transcript_slot["transcript"] is None and transcript_obj:
            transcript_meta = dict(transcript_obj)
            transcript_meta.pop("containing_exons", None)
            transcript_slot["transcript"] = transcript_meta

        protein_slot = transcript_slot["proteins"][protein_id]
        if protein_slot["protein_metadata"] is None:
            protein_slot["protein_metadata"] = get_protein_metadata(repeat)

        # --- Group by normalized repeat type ---
        repeat_type = normalize_repeat_type(repeat.get("repeatType"))
        repeat_type_slot = protein_slot["repeat_types"][repeat_type]

        # Build unique exons for this repeat type under this protein
        if "exon_id_to_obj" not in repeat_type_slot:
            repeat_type_slot["exon_id_to_obj"] = {}
        exon_id_to_obj = repeat_type_slot["exon_id_to_obj"]

        repeat_info = get_repeat_specific_info(repeat)
        repeat_exon_ids = []
        if repeat.get("ensembl_exon_info"):
            transcripts = repeat["ensembl_exon_info"].get("transcripts", [])
            for t in transcripts:
                if get_transcript_id(t) == transcript_id and t.get("containing_exons"):
                    for exon in t["containing_exons"]:
                        exon_id = exon_key(exon)
                        if exon_id and exon_id not in exon_id_to_obj:
                            exon_id_to_obj[exon_id] = exon
                        if exon_id:
                            repeat_exon_ids.append(exon_id)
                    break
        repeat_info["containing_exons"] = repeat_exon_ids
        # Store the *normalized* repeatType for reference
        repeat_info["repeatType"] = repeat_type
        repeat_type_slot["repeats"].append(repeat_info)

    # After all repeats, finalize exons under each repeat type and cleanup
    for chrom in hierarchical:
        for gene in hierarchical[chrom]:
            for transcript_id, transcript_slot in hierarchical[chrom][gene]["transcripts"].items():
                for protein_id, protein_slot in transcript_slot["proteins"].items():
                    for repeat_type, repeat_type_slot in protein_slot["repeat_types"].items():
                        if "exon_id_to_obj" in repeat_type_slot:
                            repeat_type_slot["exons"] = list(repeat_type_slot["exon_id_to_obj"].values())
                            del repeat_type_slot["exon_id_to_obj"]
                        if "containing_exons" in repeat_type_slot:
                            del repeat_type_slot["containing_exons"]
                    if "exons" in protein_slot:
                        del protein_slot["exons"]
                    if "repeats" in protein_slot:
                        del protein_slot["repeats"]

    def recursive_defaultdict_to_dict(d):
        if isinstance(d, defaultdict):
            d = {k: recursive_defaultdict_to_dict(v) for k, v in d.items()}
        return d

    result = recursive_defaultdict_to_dict(hierarchical)

    with open(output_json, "w") as f:
        json.dump(result, f, indent=2)

if __name__ == "__main__":
    # Usage: python convert_to_hierarchical.py input.json output.json
    if len(sys.argv) != 3:
        print("Usage: python convert_to_hierarchical.py input.json output.json")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
