def standardize_repeat_types(input_file, output_file):
    """
    Process a JSON file containing tandem repeat data and standardize repeat types.
    
    Args:
        input_file (str): Path to the input JSON file
        output_file (str): Path to save the modified JSON file
    """
    import json
    import re
    import logging
    from tqdm import tqdm
    
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        filename='repeat_standardization.log'
    )
    
    # Load the JSON data
    with open(input_file, 'r') as f:
        repeats = json.load(f)
    
    print(f"Processing {len(repeats)} repeat entries...")
    
    # Keep track of statistics
    missing_type_count = 0
    standardized_count = 0
    repeat_type_mapping = {}
    
    # Process each repeat entry
    for repeat in tqdm(repeats):
        original_type = repeat.get("repeatType", "")
        
        # Handle missing repeat types
        if not original_type:
            missing_type_count += 1
            repeat["repeatType"] = "Unknown"
            logging.warning(f"Missing repeatType for {repeat.get('uniProtId', 'unknown')}")
            continue
        
        # Extract base type without numbers
        standardized_type = re.sub(r'\s+\d+$', '', original_type)
        
        # Special case handling
        if "TNFR-Cys" in standardized_type:
            standardized_type = "TNFR-Cys"
        elif re.search(r'WD\d*', standardized_type):
            standardized_type = "WD"
        elif re.search(r'ANK\d*', standardized_type):
            standardized_type = "ANK"
        elif re.search(r'LRR\d*', standardized_type):
            standardized_type = "LRR"
        # Add any other special cases here
            
        # Update statistics
        if original_type != standardized_type:
            standardized_count += 1
            if original_type not in repeat_type_mapping:
                repeat_type_mapping[original_type] = standardized_type
        
        # Update the repeat type
        repeat["repeatType"] = standardized_type
    
    # Save the modified JSON
    with open(output_file, 'w') as f:
        json.dump(repeats, f, indent=2)
    
    # Print summary
    print(f"Standardized {standardized_count} repeat types")
    print(f"Found {missing_type_count} entries without repeat types")
    print(f"Saved standardized data to {output_file}")
    
    # Print the mapping for review
    print("\nRepeat Type Mapping:")
    for original, standardized in sorted(repeat_type_mapping.items()):
        print(f"  {original} -> {standardized}")

if __name__ == "__main__":
    import os
    
    # Get the directory where the script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    input_file = os.path.join(project_root, "data", "1000_test_exons_hg38_repeats.json")
    output_file = os.path.join(project_root, "data", "1000_standardized_repeats.json")
    
    standardize_repeat_types(input_file, output_file)