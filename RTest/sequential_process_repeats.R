#!/usr/bin/env Rscript

# sequential_process_repeats.R - Simple sequential processing of repeat ranges

# Load the simple process repeats script
source("/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/simple_process_repeats.R")

# Set base paths
input_file <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/DEF_gname_hg38_repeats.json"
output_dir <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/output/canonical_v2"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE)
}



# Process range 501-1000
cat("\n=== Processing range 1001-1500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "501-1000_annotated_repeats.json"),
  range = "501-1000"
)

# Process range 1001-1500
cat("\n=== Processing range 1001-1500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "1001-1500_annotated_repeats.json"),
  range = "1001-1500"
)

# Process range 1501-2000
cat("\n=== Processing range 1501-2000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "1501-2000_annotated_repeats.json"),
  range = "1501-2000"
)

# Process range 2001-2500
cat("\n=== Processing range 2001-2500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "2001-2500_annotated_repeats.json"),
  range = "2001-2500"
)

# Process range 2501-3000
cat("\n=== Processing range 2501-3000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "2501-3000_annotated_repeats.json"),
  range = "2501-3000"
)

# Process range 3001-3500
cat("\n=== Processing range 3001-3500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "3001-3500_annotated_repeats.json"),
  range = "3001-3500"
)

# Process range 3501-4000
cat("\n=== Processing range 3501-4000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "3501-4000_annotated_repeats.json"),
  range = "3501-4000"
)

# Process range 12001-12500
cat("\n=== Processing range 12001-12500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "12001-12500_annotated_repeats.json"),
  range = "12001-12500"
)

# Process range 12501-13000
cat("\n=== Processing range 12501-13000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "12501-13000_annotated_repeats.json"),
  range = "12501-13000"
)

# Process range 13001-13500
cat("\n=== Processing range 13001-13500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "13001-13500_annotated_repeats.json"),
  range = "13001-13500"
)

# Process range 13501-14000
cat("\n=== Processing range 13501-14000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "13501-14000_annotated_repeats.json"),
  range = "13501-14000"
)

# Process range 14001-14500
cat("\n=== Processing range 14001-14500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "14001-14500_annotated_repeats.json"),
  range = "14001-14500"
)

# Process range 14501-15000
cat("\n=== Processing range 14501-15000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "14501-15000_annotated_repeats.json"),
  range = "14501-15000"
)

# Process range 15001-15206
cat("\n=== Processing range 15001-15206 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "15001-15206_annotated_repeats.json"),
  range = "15001-15206"
)

cat("\nAll ranges processed successfully!\n")
