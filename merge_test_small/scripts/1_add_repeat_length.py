import json
import os
import logging
from tqdm import tqdm

# Set up logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    filename='repeat_processing.log')

def calculate_repeat_length(repeat, repeat_index=None):
    """Calculate the actual repeat length by summing all block sizes"""
    total_length = 0
    
    # Check for missing blockSizes - log and return 0
    if "blockSizes" not in repeat or not repeat["blockSizes"]:
        log_msg = f"Missing blockSizes in repeat: {repeat.get('uniProtId', 'unknown')}"
        if repeat_index is not None:
            log_msg += f" at index {repeat_index}"
        logging.warning(log_msg)
        return 0
    
    # Handle blockSizes as either a list of strings or integers
    block_sizes = repeat["blockSizes"]
    for size in block_sizes:
        # Convert to int if it's a string
        try:
            total_length += int(size)
        except (ValueError, TypeError):
            log_msg = f"Invalid blockSize value '{size}' in repeat: {repeat.get('uniProtId', 'unknown')}"
            if repeat_index is not None:
                log_msg += f" at index {repeat_index}"
            logging.warning(log_msg)
    
    return total_length

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
    missing_blocksize_count = 0
    
    # Loop through each repeat entry with progress bar
    for i, repeat in enumerate(tqdm(repeats)):
        # Calculate repeat length as sum of block sizes
        repeat_length = calculate_repeat_length(repeat, i)
        
        # Count entries with missing blockSizes
        if repeat_length == 0:
            missing_blocksize_count += 1
        
        # Skip repeats shorter than min_length
        if repeat_length < min_length:
            excluded_count += 1
            continue
        
        # Add the repeatLength field to the existing entry
        repeat["repeatLength"] = repeat_length
        
        # Add the repeat to our filtered list
        filtered_repeats.append(repeat)
    
    # Save the modified and filtered JSON
    with open(output_file, 'w') as f:
        json.dump(filtered_repeats, f, indent=2)
    
    print(f"Added repeatLength field to entries.")
    print(f"Found {missing_blocksize_count} repeats with missing blockSizes (see log)")
    print(f"Excluded {excluded_count} repeats shorter than {min_length} bp.")
    print(f"Saved {len(filtered_repeats)} repeats to {output_file}")

if __name__ == "__main__":
    # Set the input and output file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    input_file = os.path.join(project_root, "data", "hg38_repeats.json")
    output_file = os.path.join(project_root, "data", "DEF_length_filtered_hg38_repeats.json")
    
    # Process the file, filtering out repeats shorter than 60 bp
    add_repeat_length_and_filter(input_file, output_file, min_length=60)