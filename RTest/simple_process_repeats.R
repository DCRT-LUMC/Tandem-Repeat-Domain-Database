#!/usr/bin/env Rscript

# simple_process_repeats.R - A simplified version for your environment

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
query_result <- query(hub, c("ensdb", "homo sapiens", "113"))
ensdb <- hub[["AH119325"]]  # EnsDb.Hsapiens.v113
cat("Ensembl database loaded successfully.\n")

# Clean repeat type
clean_repeat_type <- function(repeat_type) {
  if (is.null(repeat_type) || is.na(repeat_type)) {
    return(repeat_type)
  }
  
  cleaned_type <- repeat_type
  
  # Check for semicolon
  if (grepl(";", cleaned_type)) {
    cleaned_type <- trimws(strsplit(cleaned_type, ";")[[1]][1])
  }
  
  # Check for space followed by a number
  parts <- strsplit(cleaned_type, " ")[[1]]
  if (length(parts) > 1 && any(grepl("^[0-9]+(\\.[0-9]+)?$", parts[-1]))) {
    cleaned_type <- parts[1]
  }
  
  # Check if result is just a number
  if (grepl("^[0-9]+(\\.[0-9]+)?$", cleaned_type)) {
    return("Unknown")
  }
  
  return(cleaned_type)
}

# Read repeat JSON file
read_repeat_json <- function(json_path) {
  repeats <- fromJSON(json_path, simplifyDataFrame = FALSE)
  return(repeats)
}

# Flatten complex objects for JSON serialization
flatten_object <- function(x) {
  if (is(x, "DataFrame") || is.data.frame(x)) {
    as.data.frame(x) %>% map(~if(is.list(.x)) paste(.x, collapse=", ") else .x)
  } else if (is.list(x)) {
    map(x, flatten_object)
  } else {
    x
  }
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
  tx <- transcripts(ensdb, filter = GRangesFilter(gr), 
                   columns = c("tx_id", "gene_id", "gene_name", "tx_biotype", 
                              "tx_name", "tx_cds_seq_start", "tx_cds_seq_end",
                              "tx_id_version", "tx_support_level", "tx_is_canonical"))
                              
  # Filter to keep only canonical transcripts
  if (length(tx) > 0) {
    tx <- tx[tx$tx_is_canonical == TRUE]
    tx <- tx[tx$tx_biotype == "protein_coding"]
  }

  # Return NULL if no transcripts found
  if (length(tx) == 0) {
    return(NULL)
  }
  
  return(tx)
}

# Get exons for a transcript and add protein coordinates
get_exons_for_transcript <- function(tx_id, ensdb) {
  # Get exons by transcript ID with all available metadata
  exons_list <- tryCatch({
      exonsBy(ensdb, by = "tx", filter = TxIdFilter(tx_id))
  }, error = function(e) {
      cat(sprintf("Warning: exonsBy failed for tx_id %s: %s\n", tx_id, conditionMessage(e)))
      NULL
  })
  
  if (is.null(exons_list) || length(exons_list) == 0 || !(tx_id %in% names(exons_list))) {
    return(NULL) # Return NULL if no exons found for this transcript
  }
  
  exons_with_meta <- exons_list[[tx_id]]
  
  # Ensure exon_rank is present, if not, add it based on order
  if (!("exon_rank" %in% colnames(mcols(exons_with_meta)))) {
      warning(paste("exon_rank missing for", tx_id, "- adding based on order."))
      exons_with_meta$exon_rank <- 1:length(exons_with_meta)
  }
  
  # Order by rank just in case
  exons_with_meta <- exons_with_meta[order(exons_with_meta$exon_rank)]
  
  # Get additional metadata if needed (like exon_id, phase - though phase calculation is separate)
  # Note: exonsBy usually returns exon_id, check if other columns are needed
  # exon_meta <- exons(ensdb, filter = TxIdFilter(tx_id), 
  #                    columns = c("exon_id", "phase", "end_phase")) # Example
  # mcols(exons_with_meta) <- cbind(mcols(exons_with_meta), mcols(exon_meta)) # If needed

  # Get protein coordinates for each exon, FILTERING by tx_id
  protein_coords_list <- lapply(1:length(exons_with_meta), function(i) {
    exon <- exons_with_meta[i]
    
    # Translate genomic coordinates to protein
    prot_list <- tryCatch({
      genomeToProtein(exon, ensdb)
    }, error = function(e) {
      # Optional: Log error for specific exon
      # cat(sprintf("Warning: genomeToProtein failed for exon rank %d (tx: %s): %s\n", mcols(exon)$exon_rank, tx_id, conditionMessage(e)))
      NULL
    })
    
    start_prot <- NA_integer_
    end_prot <- NA_integer_
    
    # Process the result (a GRangesList)
    if (!is.null(prot_list) && length(prot_list) > 0 && length(prot_list[[1]]) > 0) {
      
      # Get the GRanges object containing ALL mappings for this genomic exon
      all_prot_gr <- prot_list[[1]] 
      
      # *** FILTERING STEP ***
      # Keep only the mappings that match our target transcript ID (tx_id)
      if ("tx_id" %in% colnames(mcols(all_prot_gr))) {
          filtered_prot_gr <- all_prot_gr[mcols(all_prot_gr)$tx_id == tx_id]
      } else {
          warning(paste("tx_id column missing in genomeToProtein result for exon rank", mcols(exon)$exon_rank, "tx:", tx_id))
          filtered_prot_gr <- GRanges() # Use empty GRanges
      }
      # *** END FILTERING STEP ***

      # Now, get the IRanges from the *filtered* GRanges object
      prot_ranges <- ranges(filtered_prot_gr) 
      
      # Calculate min/max ONLY on the filtered ranges
      if (length(prot_ranges) > 0) {
        # If filtering resulted in multiple ranges *for the target transcript* 
        # (less common but possible), min/max still applies *to those*.
        start_prot <- min(start(prot_ranges)) 
        end_prot <- max(end(prot_ranges))
      } 
      # If prot_ranges is empty after filtering, start/end_prot remain NA
      
    } 
    # If prot_list was NULL or empty initially, start/end_prot remain NA

    return(data.frame(protein_start = start_prot, protein_end = end_prot))
  })
  
  # Combine the list of data frames into single columns
  coords_df <- dplyr::bind_rows(protein_coords_list)
  
  # Add protein coordinates to exons_with_meta
  mcols(exons_with_meta)$protein_start <- coords_df$protein_start
  mcols(exons_with_meta)$protein_end <- coords_df$protein_end
  
  return(exons_with_meta)
}

# Get coding status of an exon
get_coding_status <- function(exon, tx_cds_start, tx_cds_end) {
  exon_start <- start(exon)
  exon_end <- end(exon)
  
  # If no CDS information is available
  if (is.na(tx_cds_start) || is.na(tx_cds_end)) {
    return(list(status = "non_coding", utr_status = "non_coding_transcript", percentage = 0))
  }
  
  # Completely non-coding exon
  if (exon_end <= tx_cds_start || exon_start >= tx_cds_end) {
    # Determine if it's 5' or 3' UTR based on strand
    exon_strand <- as.character(strand(exon))
    
    if (exon_end <= tx_cds_start) {
      utr_type <- if (exon_strand == "+") "5'_UTR" else "3'_UTR"
    } else {
      utr_type <- if (exon_strand == "+") "3'_UTR" else "5'_UTR"
    }
    
    return(list(status = "non_coding", utr_status = utr_type, percentage = 0))
  }
  
  # Partially coding exon
  if (exon_start < tx_cds_start || exon_end > tx_cds_end) {
    utr_regions <- c()
    exon_strand <- as.character(strand(exon))
    
    if (exon_start < tx_cds_start) {
      utr_regions <- c(utr_regions, if (exon_strand == "+") "5'_UTR" else "3'_UTR")
    }
    
    if (exon_end > tx_cds_end) {
      utr_regions <- c(utr_regions, if (exon_strand == "+") "3'_UTR" else "5'_UTR")
    }
    
    # Calculate coding percentage
    exon_length <- exon_end - exon_start
    coding_start <- max(exon_start, tx_cds_start)
    coding_end <- min(exon_end, tx_cds_end)
    coding_length <- coding_end - coding_start
    coding_percentage <- (coding_length / exon_length) * 100
    
    utr_status <- paste(utr_regions, collapse = "+")
    
    return(list(status = "partial_coding", utr_status = utr_status, percentage = round(coding_percentage, 2)))
  }
  
  # Fully coding exon
  return(list(status = "fully_coding", utr_status = "none", percentage = 100))
}

# Determine frame status based on phase information
get_frame_status <- function(phase, end_phase, coding_status) {
  if (coding_status == "non_coding") {
    return("non_coding")
  }
  
  if (is.na(phase) || is.na(end_phase) || phase == -1 || end_phase == -1) {
    return("unknown")
  }
  
  if (phase == end_phase) {
    return("in_frame")
  }
  
  return("out_of_frame")
}

# Extract protein coordinates from position field
extract_protein_position <- function(position_text) {
  if (is.null(position_text) || is.na(position_text) || position_text == "") {
    return(list(protein_start = NA_integer_, protein_end = NA_integer_))
  }
  
  # Match pattern like "amino acids 438-483 on protein Q6TDP4"
  pattern <- "amino acids ([0-9]+)-([0-9]+) on protein"
  match <- regexec(pattern, position_text)
  if (match[[1]][1] != -1) {
    matches <- regmatches(position_text, match)[[1]]
    protein_start <- as.integer(matches[2])
    protein_end <- as.integer(matches[3])
    return(list(protein_start = protein_start, protein_end = protein_end))
  }
  
  return(list(protein_start = NA_integer_, protein_end = NA_integer_))
}

# Calculate phase from CDS positions - more accurate than database retrieval
calculate_phase <- function(exon, cds_regions) {
  # Check if there are any CDS regions
  if (length(cds_regions) == 0) {
    return(list(phase = -1, end_phase = -1))
  }
  
  # Find overlapping CDS regions
  overlaps <- findOverlaps(exon, cds_regions)
  
  if (length(overlaps) == 0) {
    # Non-coding exon
    return(list(phase = -1, end_phase = -1))
  }
  
  # Get the cumulative CDS length before this exon
  all_cds_before <- cds_regions[1:subjectHits(overlaps)[1]-1]
  cds_pos <- sum(width(all_cds_before))
  
  # Calculate phase
  phase <- cds_pos %% 3
  
  # Calculate end phase
  overlap <- pintersect(exon, cds_regions[subjectHits(overlaps)])
  overlap_length <- sum(width(overlap))
  end_phase <- (cds_pos + overlap_length) %% 3
  
  return(list(phase = phase, end_phase = end_phase))
}

# Process a single repeat
process_repeat <- function(repeat_data, ensdb) {
  # Get CDS regions by transcript - add this line!
  cds_by_tx <- cdsBy(ensdb, by="tx")
  
  # Extract basic information
  chrom <- repeat_data$chrom
  start <- as.integer(repeat_data$chromStart)
  end <- as.integer(repeat_data$chromEnd)
  strand <- repeat_data$strand
  
  # Extract protein coordinates from position field
  if (!is.null(repeat_data$position)) {
    protein_pos <- extract_protein_position(repeat_data$position)
    repeat_data$protein_start <- protein_pos$protein_start
    repeat_data$protein_end <- protein_pos$protein_end
  }
  
  # Clean the repeat type
  if (!is.null(repeat_data$repeatType)) {
    repeat_data$repeatType <- clean_repeat_type(repeat_data$repeatType)
  }
  
  # Get overlapping transcripts
  transcripts_result <- get_overlapping_transcripts(chrom, start, end, strand, ensdb)
  
  # Prepare result structure
  ensembl_exon_info <- list(
    transcripts_count = 0,
    has_canonical_transcript = FALSE,
    location_summary = "unknown",
    transcripts = list()
  )
  
  # If no transcripts overlap, return empty result
  if (is.null(transcripts_result) || length(transcripts_result) == 0) {
    repeat_data$ensembl_exon_info <- ensembl_exon_info
    return(repeat_data)
  }
  
  # Process each transcript
  transcript_info <- list()
  locations <- c()
  canonical_found <- FALSE
  
  for (i in 1:length(transcripts_result)) {
    transcript <- transcripts_result[i]
    tx_id <- transcript$tx_id
    
    # Skip if strand doesn't match (when strand is specified)
    if (!is.null(strand) && !is.na(strand) && strand %in% c("+", "-")) {
      tx_strand <- as.character(strand(transcript))
      if (strand != tx_strand) {
        next
      }
    }
    
    # Get exons for this transcript
    exons <- get_exons_for_transcript(tx_id, ensdb)
    if (is.null(exons) || length(exons) == 0) {
      next
    }
    
    # Basic transcript info
    exon_count <- length(exons)
    tx_strand <- as.character(strand(transcript))
    tx_biotype <- transcript$tx_biotype
    
    # Check if this is a canonical transcript
    is_canonical <- FALSE
    if (!is.null(transcript$tx_is_canonical) && !is.na(transcript$tx_is_canonical)) {
      is_canonical <- as.logical(transcript$tx_is_canonical)
      if (is_canonical) {
        canonical_found <- TRUE
      }
    }
    
    # Use the directly fetched versioned transcript ID
    versioned_tx_id <- transcript$tx_id_version # Use tx_id_version directly
    
    # Determine location (exonic, intronic, outside)
    tx_start <- start(transcript)
    tx_end <- end(transcript)
    location <- "unknown"
    
    # Check if outside transcript bounds
    if (end <= tx_start || start >= tx_end) {
      location <- "outside"
    } else {
      # Check if overlapping any exon
      is_exonic <- FALSE
      for (j in 1:length(exons)) {
        exon_start <- start(exons[j])
        exon_end <- end(exons[j])
        
        if (max(start, exon_start) < min(end, exon_end)) {
          is_exonic <- TRUE
          break
        }
      }
      
      if (is_exonic) {
        location <- "exonic"
      } else {
        location <- "intronic"
      }
    }
    
    # Skip intronic transcripts - only continue processing for exonic or outside transcripts
    if (location != "intronic") {
      locations <- c(locations, location)
      
      # Get CDS information
      tx_cds_start <- transcript$tx_cds_seq_start
      tx_cds_end <- transcript$tx_cds_seq_end
      
      # Get CDS regions for this transcript
      cds_regions <- NULL
      if (tx_id %in% names(cds_by_tx)) {
        cds_regions <- cds_by_tx[[tx_id]]
      } else {
        cds_regions <- GRanges()  # Empty GRanges if no CDS found
      }
      
      # Process exons only for exonic transcripts
      containing_exons <- list()
      if (location == "exonic") {
        for (j in 1:length(exons)) {
          exon <- exons[j]
          exon_start <- start(exon)
          exon_end <- end(exon)
          
          # Check if repeat overlaps this exon
          if (max(start, exon_start) < min(end, exon_end)) {
            exon_number <- j
            
            # Calculate overlap statistics
            overlap_start <- max(start, exon_start)
            overlap_end <- min(end, exon_end)
            overlap_length <- overlap_end - overlap_start
            exon_length <- exon_end - exon_start
            overlap_percentage <- (overlap_length / exon_length) * 100
            
            # Determine position in transcript
            if (exon_count == 1) {
              position <- "single_exon"
            } else if (j == 1) {
              position <- "first_exon"
            } else if (j == exon_count) {
              position <- "last_exon"
            } else {
              position <- paste0("middle_exon_", exon_number)
            }
            
            # Get coding status for this exon
            coding_info <- get_coding_status(exon, tx_cds_start, tx_cds_end)
            
            # Calculate phase information from CDS positions instead of database retrieval
            phase_info <- calculate_phase(exon, cds_regions)
            phase <- phase_info$phase
            end_phase <- phase_info$end_phase
            
            # Get frame status
            frame_status <- get_frame_status(phase, end_phase, coding_info$status)
            
            # Create exon info object
            exon_info <- list(
              exon_number = exon_number,
              exon_id = exon$exon_id,
              exon_start = exon_start,
              exon_end = exon_end,
              protein_start = exon$protein_start,
              protein_end = exon$protein_end,
              overlap_bp = overlap_length,
              position = position,
              overlap_percentage = round(overlap_percentage, 2),
              coding_status = coding_info$status,
              utr_status = coding_info$utr_status,
              coding_percentage = coding_info$percentage,
              phase = phase,
              end_phase = end_phase,
              frame_status = frame_status
            )
            
            containing_exons <- c(containing_exons, list(exon_info))
          }
        }
      }
      
      # Create transcript info object
      tx_info <- list(
        transcript_id = tx_id,
        versioned_transcript_id = versioned_tx_id,
        transcript_name = transcript$tx_name,
        is_canonical = is_canonical,
        biotype = tx_biotype,
        location = location,
        exon_count = exon_count,
        strand = tx_strand,
        containing_exons = containing_exons
      )
      
      transcript_info <- c(transcript_info, list(tx_info))
    }
  }
  
  # Summarize location (prioritize exonic > outside/unknown)
  location_summary <- "unknown"
  if ("exonic" %in% locations) {
    location_summary <- "exonic"
  } else if ("outside" %in% locations) {
    location_summary <- "intergenic"
  }
  
  # Update ensembl_exon_info
  ensembl_exon_info$transcripts_count <- length(transcript_info)
  ensembl_exon_info$has_canonical_transcript <- canonical_found
  ensembl_exon_info$location_summary <- location_summary
  ensembl_exon_info$transcripts <- transcript_info
  
  # Add to repeat data
  repeat_data$ensembl_exon_info <- ensembl_exon_info
  
  return(repeat_data)
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
    # Parse range in format "start-end"
    range_parts <- as.integer(unlist(strsplit(range, "-")))
    if (length(range_parts) == 2) {
      start_idx <- range_parts[1]
      end_idx <- range_parts[2]
      
      # Validate range
      if (start_idx >= 1 && end_idx <= length(valid_repeats) && start_idx <= end_idx) {
        valid_repeats <- valid_repeats[start_idx:end_idx]
        cat(sprintf("Processing repeats %d to %d out of %d repeats...\n", 
                   start_idx, end_idx, length(repeats)))
      } else {
        cat("Invalid range specified. Using all repeats.\n")
      }
    } else {
      cat("Invalid range format. Use 'start-end' (e.g., '1-1000'). Using all repeats.\n")
    }
  } else if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    valid_repeats <- valid_repeats[1:min(limit, length(valid_repeats))]
    cat("Processing first", length(valid_repeats), "out of", length(repeats), "repeats...\n")
  } else {
    cat("Processing", length(valid_repeats), "out of", length(repeats), "repeats with valid coordinates...\n")
  }
  
  # Try parallel processing with BiocParallel
  num_cores <- 4
  #num_cores <- min(detectCores() - 1, 16)  # More conservative core usage
  cat("Using", num_cores, "cores for processing with BiocParallel...\n")  # Fix typo: num_ores → num_cores
  
  # First try parallel processing
  processed_repeats <- NULL
  tryCatch({
    # Configure BiocParallel to use multiple cores
    bp_param <- MulticoreParam(
      workers = num_cores,
      progressbar = TRUE,
      stop.on.error = FALSE,
      tasks = min(100, length(valid_repeats))  # Process in batches
    )
    
    # Register the parameters
    register(bp_param)
    
    # Process repeats in parallel with BiocParallel
    cat("Processing repeats in parallel...\n")
    processed_repeats <- bplapply(valid_repeats, function(repeat_data) {
      # Each worker will have its own environment with necessary functions
      result <- process_repeat(repeat_data, ensdb)
      return(result)
    }, BPPARAM = bp_param)
    
    cat("Parallel processing completed successfully.\n")
  }, error = function(e) {
    cat("Parallel processing failed with error:", conditionMessage(e), "\n")
    cat("Falling back to sequential processing...\n")
    
    # Process sequentially
    processed_repeats <<- lapply(valid_repeats, function(repeat_data) {
      cat(".")  # Simple progress indicator
      result <- process_repeat(repeat_data, ensdb)
      return(result)
    })
    cat("\nSequential processing completed.\n")
  })
  
  # Write results to JSON file
  cat("Writing results to", output_file, "...\n")
  write_json(processed_repeats, output_file, auto_unbox = TRUE, pretty = TRUE)
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("Processing completed in", format(time_taken), "\n")
  
  # Return invisible to prevent automatic printing
  return(invisible(processed_repeats))
}

# Add command line argument handling
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript simple_process_repeats.R <input_json> <output_json> [limit|range]
          Where range is in format 'start-end' (e.g., '1-1000')")
  }
  
  input_file <- args[1]
  output_file <- args[2]
  
  # Check if third argument is a range or a limit
  if (length(args) >= 3) {
    if (grepl("-", args[3])) {
      # It's a range parameter
      range <- args[3]
      simple_process_repeat_data(input_file, output_file, range = range)
    } else {
      # It's a limit parameter
      limit <- as.numeric(args[3])
      simple_process_repeat_data(input_file, output_file, limit = limit)
    }
  } else {
    simple_process_repeat_data(input_file, output_file)
  }
}

# Example usage when running directly in R
# To process repeats in the range 100-200:
# simple_process_repeat_data(
#   input_file = "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/DEF_gname_hg38_repeats.json",
#   output_file = "/home/dogdorgesh/Documents/GitHub/Tandem-Repeat-Domain-Database/RTest/test_annotated_repeats.json",
#   range = "1-10"  # IMPORTANT: range must be a string with quotes
# )

# Incorrect usage (will cause error):
# simple_process_repeat_data(
#   input_file = "path/to/input.json",
#   output_file = "path/to/output.json",
#   range = 1-10  # ERROR: This evaluates to -9, not the string "1-10"
# )
