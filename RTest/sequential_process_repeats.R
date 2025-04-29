#!/usr/bin/env Rscript

# sequential_process_repeats.R - Simple sequential processing of repeat ranges

# Load the simple process repeats script
source("/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/simple_process_repeats.R")

# Set base paths
input_file <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/DEF_gname_hg38_repeats.json"
output_dir <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/output"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  cat("Creating output directory:", output_dir, "\n")
  dir.create(output_dir, recursive = TRUE)
}

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

# Process range 4001-4500
cat("\n=== Processing range 4001-4500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "4001-4500_annotated_repeats.json"),
  range = "4001-4500"
)

# Process range 4501-5000
cat("\n=== Processing range 4501-5000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "4501-5000_annotated_repeats.json"),
  range = "4501-5000"
)

# Process range 5001-5500
cat("\n=== Processing range 5001-5500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "5001-5500_annotated_repeats.json"),
  range = "5001-5500"
)

# Process range 5501-6000
cat("\n=== Processing range 5501-6000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "5501-6000_annotated_repeats.json"),
  range = "5501-6000"
)

# Process range 6001-6500
cat("\n=== Processing range 6001-6500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "6001-6500_annotated_repeats.json"),
  range = "6001-6500"
)

# Process range 6501-7000
cat("\n=== Processing range 6501-7000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "6501-7000_annotated_repeats.json"),
  range = "6501-7000"
)

# Process range 7001-7500
cat("\n=== Processing range 7001-7500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "7001-7500_annotated_repeats.json"),
  range = "7001-7500"
)

# Process range 7501-8000
cat("\n=== Processing range 7501-8000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "7501-8000_annotated_repeats.json"),
  range = "7501-8000"
)

# Process range 8001-8500
cat("\n=== Processing range 8001-8500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "8001-8500_annotated_repeats.json"),
  range = "8001-8500"
)

# Process range 8501-9000
cat("\n=== Processing range 8501-9000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "8501-9000_annotated_repeats.json"),
  range = "8501-9000"
)

# Process range 9001-9500
cat("\n=== Processing range 9001-9500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "9001-9500_annotated_repeats.json"),
  range = "9001-9500"
)

# Process range 9501-10000
cat("\n=== Processing range 9501-10000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "9501-10000_annotated_repeats.json"),
  range = "9501-10000"
)

# Process range 10001-10500
cat("\n=== Processing range 10001-10500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "10001-10500_annotated_repeats.json"),
  range = "10001-10500"
)

# Process range 10501-11000
cat("\n=== Processing range 10501-11000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "10501-11000_annotated_repeats.json"),
  range = "10501-11000"
)

# Process range 11001-11500
cat("\n=== Processing range 11001-11500 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "11001-11500_annotated_repeats.json"),
  range = "11001-11500"
)

# Process range 11501-12000
cat("\n=== Processing range 11501-12000 ===\n")
simple_process_repeat_data(
  input_file = input_file,
  output_file = file.path(output_dir, "11501-12000_annotated_repeats.json"),
  range = "11501-12000"
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
