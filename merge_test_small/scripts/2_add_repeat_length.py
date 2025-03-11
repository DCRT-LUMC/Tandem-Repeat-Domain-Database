import json
import os

def add_repeat_length_and_filter(input_file, output_file, min_length=60):
    """
    Adds a repeatLength field to each repeat entry and filters out repeats 
    shorter than the minimum length.
    
    Args:
        input_file (str): Path to the input JSON file
        output_file (str): Path to save the modified JSON file
        min_length (int): Minimum repeat length to include in output (default: 60)
    """
    # Load the JSON data
    with open(input_file, 'r') as f:
        repeats = json.load(f)
    
    print(f"Processing {len(repeats)} repeat entries...")
    
    # Create a new list for filtered repeats
    filtered_repeats = []
    excluded_count = 0
    
    # Loop through each repeat entry
    for repeat in repeats:
        # Calculate repeat length
        chrom_start = repeat.get("chromStart", 0)
        chrom_end = repeat.get("chromEnd", 0)
        repeat_length = chrom_end - chrom_start
        
        # Skip repeats shorter than min_length
        if repeat_length < min_length:
            excluded_count += 1
            continue
        
        # Create a new dictionary with the repeatLength field inserted after "comments"
        new_repeat = {}
        for key, value in repeat.items():
            new_repeat[key] = value
            # Insert repeatLength right after comments
            if key == "comments":
                new_repeat["repeatLength"] = repeat_length
        
        # If "comments" wasn't found, add it at the end
        if "repeatLength" not in new_repeat:
            new_repeat["repeatLength"] = repeat_length
        
        # Add the modified repeat to our filtered list
        filtered_repeats.append(new_repeat)
    
    # Save the modified and filtered JSON
    with open(output_file, 'w') as f:
        json.dump(filtered_repeats, f, indent=2)
    
    print(f"Added repeatLength field to entries and excluded {excluded_count} repeats shorter than {min_length} bp.")
    print(f"Saved {len(filtered_repeats)} repeats to {output_file}")

if __name__ == "__main__":
    # Set the input and output file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    input_file = os.path.join(project_root, "data", "gname_hg38_repeats_100.json")
    output_file = os.path.join(project_root, "data", "length_filtered_gname_hg38_repeats.json")
    
    # Process the file, filtering out repeats shorter than 60 bp
    add_repeat_length_and_filter(input_file, output_file, min_length=60)