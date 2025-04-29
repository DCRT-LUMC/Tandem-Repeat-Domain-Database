#!/usr/bin/env Rscript

# batch_process_repeats.R - Process repeats in batches of 500

# First, load the necessary libraries
library(jsonlite)
library(ensembldb)
library(AnnotationHub)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(BiocGenerics)
library(IRanges)
library(parallel)
library(BiocParallel)

# Define the function separately to avoid command line argument execution
source_without_executing <- function(path) {
  env <- new.env()
  sys.source(path, env)
  return(env)
}

# Load the script functions without executing the command line portion
script_env <- source_without_executing("/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/simple_process_repeats.R")
simple_process_repeat_data <- script_env$simple_process_repeat_data
read_repeat_json <- script_env$read_repeat_json

# Input file path
input_file <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/DEF_gname_hg38_repeats.json"

# Base output directory
output_dir <- "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/output"

# Batch size
batch_size <- 500

# Get total number of repeats in the file
cat("Determining the total number of repeats...\n")
repeats <- read_repeat_json(input_file)
valid_repeats <- repeats[sapply(repeats, function(r) {
  !is.null(r$chrom) && !is.null(r$chromStart) && !is.null(r$chromEnd)
})]
total_repeats <- length(valid_repeats)
cat("Found", total_repeats, "valid repeats.\n")

# Calculate number of batches
num_batches <- ceiling(total_repeats / batch_size)
cat("Processing will be done in", num_batches, "batches of", batch_size, "repeats each.\n")
cat("Starting from batch 3 (index 1001) since earlier batches are already processed.\n")

# Process each batch starting from batch 3
for (batch in 3:num_batches) {
  # Calculate start and end indices for this batch
  start_idx <- (batch - 1) * batch_size + 1
  end_idx <- min(batch * batch_size, total_repeats)
  
  # Generate range string
  range_str <- paste0(start_idx, "-", end_idx)
  
  # Generate output filename
  output_file <- file.path(output_dir, paste0(range_str, "_annotated_repeats.json"))
  
  cat("\n==========================================================\n")
  cat("Processing batch", batch, "of", num_batches, ":", range_str, "\n")
  cat("Output will be saved to:", output_file, "\n")
  cat("==========================================================\n")
  
  # Process this batch
  simple_process_repeat_data(
    input_file = input_file,
    output_file = output_file,
    range = range_str
  )
  
  # Force garbage collection to free memory between batches
  gc()
  
  cat("Batch", batch, "completed.\n")
}

cat("\nAll batches processed successfully!\n")
