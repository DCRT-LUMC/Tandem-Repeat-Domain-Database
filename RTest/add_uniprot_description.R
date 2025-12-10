#!/usr/bin/env Rscript

# add_uniprot_description.R
# Adds DESCRIPTION and UNIPROTID from EnsDb to repeat data

# Load required libraries
library(ensembldb)
library(AnnotationHub)
library(jsonlite)
library(GenomicRanges)
library(dplyr)
library(purrr)
library(BiocGenerics)
library(IRanges)
library(parallel)
library(BiocParallel)

# Load the Ensembl database
cat("Loading EnsDb.Hsapiens.v113 from AnnotationHub...\n")
hub <- AnnotationHub()
# query_result <- query(hub, c("ensdb", "homo sapiens", "113"))
ensdb <- hub[["AH119325"]]  # EnsDb.Hsapiens.v113
cat("Ensembl database loaded successfully.\n")

# Read repeat JSON file
read_repeat_json <- function(json_path) {
  repeats <- fromJSON(json_path, simplifyDataFrame = FALSE)
  return(repeats)
}

# Get overlapping transcripts for a genomic region
get_overlapping_transcripts <- function(chrom, start, end, strand = NULL, ensdb) {
  # Remove 'chr' prefix if present and convert to Ensembl format
  chrom <- gsub("chr", "", chrom)
  
  # Create genomic range for the region
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
  
  # Add strand information if available
  if (!is.null(strand) && !is.na(strand) && strand %in% c("+", "-")) {
    strand(gr) <- strand
  }
  
  # Get all transcripts overlapping the region
  # We only need gene_id and tx_id really, but keeping others for filtering
  tx <- transcripts(ensdb, filter = GRangesFilter(gr), 
                   columns = c("tx_id", "gene_id", "tx_biotype", "tx_is_canonical"))
                              
  # Filter to keep only canonical transcripts
  if (length(tx) > 0) {
    # Handle potential NA in tx_is_canonical
    tx <- tx[!is.na(tx$tx_is_canonical) & tx$tx_is_canonical == TRUE]
    tx <- tx[!is.na(tx$tx_biotype) & tx$tx_biotype == "protein_coding"]
  }

  # Return NULL if no transcripts found
  if (length(tx) == 0) {
    return(NULL)
  }
  
  return(tx)
}

# Process a single repeat
process_repeat <- function(repeat_data, ensdb) {
  chrom <- repeat_data$chrom
  start <- as.integer(repeat_data$chromStart)
  end <- as.integer(repeat_data$chromEnd)
  strand <- repeat_data$strand
  
  # Get overlapping transcripts
  transcripts_result <- get_overlapping_transcripts(chrom, start, end, strand, ensdb)
  
  description <- NA
  uniProtIdCanonical <- NA
  
  if (!is.null(transcripts_result) && length(transcripts_result) > 0) {
    # Use the first canonical transcript found
    transcript <- transcripts_result[1]
    gene_id <- transcript$gene_id
    
    # Fetch DESCRIPTION and UNIPROTID
    tryCatch({
      # Check available columns to avoid errors
      avail_cols <- columns(ensdb)
      cols_to_fetch <- c()
      if ("DESCRIPTION" %in% avail_cols) cols_to_fetch <- c(cols_to_fetch, "DESCRIPTION")
      if ("UNIPROTID" %in% avail_cols) cols_to_fetch <- c(cols_to_fetch, "UNIPROTID")
      
      if (length(cols_to_fetch) > 0) {
        info <- select(ensdb, keys = gene_id, keytype = "GENEID", columns = cols_to_fetch)
        
        if (nrow(info) > 0) {
          if ("DESCRIPTION" %in% colnames(info)) {
             desc <- info$DESCRIPTION[1]
             if (!is.na(desc)) {
               description <- desc
             }
          }
          
          if ("UNIPROTID" %in% colnames(info)) {
            uids <- unique(info$UNIPROTID[!is.na(info$UNIPROTID)])
            if (length(uids) > 0) {
              uniProtIdCanonical <- paste(uids, collapse = ";")
            }
          }
        }
      }
    }, error = function(e) {
      # Silent fail
    })
  }
  
  # Reconstruct list to insert fields after uniProtId
  new_repeat_data <- list()
  found_target <- FALSE
  
  keys <- names(repeat_data)
  for (k in keys) {
    new_repeat_data[[k]] <- repeat_data[[k]]
    if (k == "uniProtId") {
      new_repeat_data[["description"]] <- description
      new_repeat_data[["uniProtIdCanonical"]] <- uniProtIdCanonical
      found_target <- TRUE
    }
  }
  
  if (!found_target) {
    new_repeat_data[["description"]] <- description
    new_repeat_data[["uniProtIdCanonical"]] <- uniProtIdCanonical
  }
  
  return(new_repeat_data)
}

# Main processing function
simple_process_repeat_data <- function(input_file, output_file, limit = NULL, range = NULL) {
  start_time <- Sys.time()
  
  # Load repeat data
  cat("Loading repeat data from", input_file, "...\n")
  repeats <- read_repeat_json(input_file)
  
  # Filter repeats with valid coordinates
  valid_repeats <- repeats[sapply(repeats, function(r) {
    !is.null(r$chrom) && !is.null(r$chromStart) && !is.null(r$chromEnd)
  })]
  
  # Apply range or limit if specified
  if (!is.null(range)) {
    range_parts <- as.integer(unlist(strsplit(range, "-")))
    if (length(range_parts) == 2) {
      start_idx <- range_parts[1]
      end_idx <- range_parts[2]
      if (start_idx >= 1 && end_idx <= length(valid_repeats) && start_idx <= end_idx) {
        valid_repeats <- valid_repeats[start_idx:end_idx]
        cat(sprintf("Processing repeats %d to %d out of %d repeats...\n", start_idx, end_idx, length(repeats)))
      }
    }
  } else if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    valid_repeats <- valid_repeats[1:min(limit, length(valid_repeats))]
    cat("Processing first", length(valid_repeats), "out of", length(repeats), "repeats...\n")
  }
  
  # Parallel processing setup
  num_cores <- 12
  cat("Using", num_cores, "cores for processing with BiocParallel...\n")
  
  processed_repeats <- NULL
  tryCatch({
    bp_param <- MulticoreParam(workers = num_cores, progressbar = TRUE, stop.on.error = FALSE, tasks = min(100, length(valid_repeats)))
    register(bp_param)
    
    cat("Processing repeats in parallel...\n")
    processed_repeats <- bplapply(valid_repeats, function(repeat_data) {
      result <- process_repeat(repeat_data, ensdb)
      return(result)
    }, BPPARAM = bp_param)
    
    cat("Parallel processing completed successfully.\n")
  }, error = function(e) {
    cat("Parallel processing failed:", conditionMessage(e), "\nFalling back to sequential...\n")
    processed_repeats <<- lapply(valid_repeats, function(repeat_data) {
      process_repeat(repeat_data, ensdb)
    })
  })
  
  cat("Writing results to", output_file, "...\n")
  write_json(processed_repeats, output_file, auto_unbox = TRUE, pretty = TRUE)
  
  end_time <- Sys.time()
  cat("Processing completed in", format(end_time - start_time), "\n")
}

# Command line args
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript add_uniprot_description.R <input_json> <output_json> [limit|range]")
  }
  
  input_file <- args[1]
  output_file <- args[2]
  
  if (length(args) >= 3) {
    if (grepl("-", args[3])) {
      simple_process_repeat_data(input_file, output_file, range = args[3])
    } else {
      simple_process_repeat_data(input_file, output_file, limit = as.numeric(args[3]))
    }
  } else {
    simple_process_repeat_data(input_file, output_file)
  }
}
