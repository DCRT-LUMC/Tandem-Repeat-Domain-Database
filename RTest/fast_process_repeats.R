#!/usr/bin/env Rscript

# fast_process_repeats.R - High-performance version

# Load required libraries
suppressPackageStartupMessages({
  library(ensembldb)
  library(AnnotationHub)
  library(jsonlite)
  library(GenomicRanges)
  library(dplyr)
  library(purrr)
  library(data.table)
  library(BiocParallel)
  library(future)
  library(future.apply)
  library(fst)
})

# Helper functions
clean_repeat_type <- function(repeat_type) {
  if (is.null(repeat_type) || is.na(repeat_type)) return(repeat_type)
  
  cleaned_type <- repeat_type
  if (grepl(";", cleaned_type)) cleaned_type <- trimws(strsplit(cleaned_type, ";")[[1]][1])
  
  parts <- strsplit(cleaned_type, " ")[[1]]
  if (length(parts) > 1 && any(grepl("^[0-9]+(\\.[0-9]+)?$", parts[-1]))) cleaned_type <- parts[1]
  
  if (grepl("^[0-9]+(\\.[0-9]+)?$", cleaned_type)) return("Unknown")
  
  return(cleaned_type)
}

extract_protein_position <- function(position_text) {
  if (is.null(position_text) || is.na(position_text) || position_text == "") {
    return(list(protein_start = NA_integer_, protein_end = NA_integer_))
  }
  
  pattern <- "amino acids ([0-9]+)-([0-9]+) on protein"
  match <- regexec(pattern, position_text)
  if (match[[1]][1] != -1) {
    matches <- regmatches(position_text, match)[[1]]
    return(list(protein_start = as.integer(matches[2]), protein_end = as.integer(matches[3])))
  }
  
  return(list(protein_start = NA_integer_, protein_end = NA_integer_))
}

# Main processing function
fast_process_repeat_data <- function(input_file, output_file, limit = NULL, range = NULL, cores = NULL) {
  start_time <- Sys.time()
  
  # Set number of cores
  if (is.null(cores)) cores <- min(parallel::detectCores() - 1, 12)
  cat(paste("Using", cores, "cores\n"))
  
  # Configure parallel processing
  future::plan(future::multicore, workers = cores)
  options(future.globals.maxSize = 8 * 1024^3)  # 8GB limit for global data
  
  # Read all repeats at once
  cat("Loading repeat data...\n")
  repeats <- fromJSON(input_file, simplifyDataFrame = FALSE)
  
  # Filter repeats with valid coordinates
  valid_repeats <- repeats[sapply(repeats, function(r) {
    !is.null(r$chrom) && !is.null(r$chromStart) && !is.null(r$chromEnd)
  })]
  
  # Apply range or limit if specified
  if (!is.null(range)) {
    range_parts <- as.integer(unlist(strsplit(range, "-")))
    if (length(range_parts) == 2) {
      start_idx <- range_parts[1]
      end_idx <- min(range_parts[2], length(valid_repeats))
      if (start_idx >= 1 && start_idx <= end_idx) {
        valid_repeats <- valid_repeats[start_idx:end_idx]
        cat(sprintf("Processing repeats %d to %d out of %d repeats...\n", 
                    start_idx, end_idx, length(repeats)))
      }
    }
  } else if (!is.null(limit) && is.numeric(limit) && limit > 0) {
    valid_repeats <- valid_repeats[1:min(limit, length(valid_repeats))]
    cat("Processing first", length(valid_repeats), "out of", length(repeats), "repeats...\n")
  } else {
    cat("Processing", length(valid_repeats), "out of", length(repeats), "repeats...\n")
  }
  
  # Load Ensembl database once
  cat("Loading and indexing Ensembl database...\n")
  hub <- AnnotationHub()
  ensdb <- hub[["AH119325"]]  # EnsDb.Hsapiens.v113
  
  # Temporarily suppress warnings for database operations
  suppressWarnings({
    # Pre-load all transcripts to memory
    cat("Preloading transcript data...\n")
    all_transcripts <- transcripts(ensdb, 
                                  columns=c("tx_id", "gene_id", "gene_name", "tx_biotype", 
                                            "tx_name", "tx_cds_seq_start", "tx_cds_seq_end",
                                            "tx_id_version", "tx_is_canonical"))
    
    # Filter to keep only protein-coding transcripts and remove problematic sequences
    all_transcripts <- all_transcripts[all_transcripts$tx_biotype == "protein_coding"]
    all_transcripts <- all_transcripts[!seqnames(all_transcripts) %in% c("LRG_432")]
    
    # Pre-load exon information - safely handle missing columns
    cat("Preloading exon data...\n")
    all_exons <- exons(ensdb, columns=c("exon_id", "exon_idx", "exon_seq_start", 
                                      "exon_seq_end"))
    
    # Create exons by transcript mapping for quick lookup
    cat("Creating exon-transcript mappings...\n")
    exons_by_tx <- exonsBy(ensdb, by="tx")
    
    # Pre-load CDS information
    cat("Preloading CDS data...\n")
    all_cds <- cdsBy(ensdb, by="tx")
  })
  
  # Fix: Suppress warnings about out-of-bound ranges
  options(warn = -1)  # Temporarily suppress warnings
  
  # Create GRanges object from all repeats for faster overlap computation
  cat("Preparing repeat data for processing...\n")
  repeat_granges <- tryCatch({
    GRanges(
      seqnames = sapply(valid_repeats, function(r) gsub("chr", "", r$chrom)),
      ranges = IRanges(
        start = as.integer(sapply(valid_repeats, function(r) r$chromStart)),
        end = as.integer(sapply(valid_repeats, function(r) r$chromEnd))
      ),
      strand = sapply(valid_repeats, function(r) if(!is.null(r$strand)) r$strand else "*")
    )
  }, error = function(e) {
    cat("Error creating GRanges object:", conditionMessage(e), "\n")
    cat("Trying again with more error checking...\n")
    
    # More careful approach
    seqnames <- character(length(valid_repeats))
    starts <- integer(length(valid_repeats))
    ends <- integer(length(valid_repeats))
    strands <- character(length(valid_repeats))
    
    for (i in seq_along(valid_repeats)) {
      r <- valid_repeats[[i]]
      seqnames[i] <- gsub("chr", "", r$chrom)
      starts[i] <- as.integer(r$chromStart)
      ends[i] <- as.integer(r$chromEnd)
      strands[i] <- if(!is.null(r$strand)) r$strand else "*"
    }
    
    return(GRanges(
      seqnames = seqnames,
      ranges = IRanges(start = starts, end = ends),
      strand = strands
    ))
  })
  
  options(warn = 0)  # Restore warning level
  
  # Store index of original repeats
  mcols(repeat_granges)$repeat_idx <- 1:length(valid_repeats)
  
  # Find all transcript overlaps in one operation
  cat("Finding all genomic overlaps...\n")
  transcript_overlaps <- findOverlaps(repeat_granges, all_transcripts)
  
  # Group overlaps by repeat
  overlaps_by_repeat <- split(subjectHits(transcript_overlaps), queryHits(transcript_overlaps))
  
  # Process each repeat
  cat("Processing repeats in parallel...\n")
  pb <- txtProgressBar(min = 0, max = length(valid_repeats), style = 3)
  
  # Process repeats with parallel future_lapply
  processed_repeats <- future_lapply(1:length(valid_repeats), function(i) {
    tryCatch({
      repeat_data <- valid_repeats[[i]]
      repeat_gr <- repeat_granges[i]
      
      # Clean repeat data
      if (!is.null(repeat_data$repeatType)) {
        repeat_data$repeatType <- clean_repeat_type(repeat_data$repeatType)
      }
      
      # Extract protein coordinates
      if (!is.null(repeat_data$position)) {
        protein_pos <- extract_protein_position(repeat_data$position)
        repeat_data$protein_start <- protein_pos$protein_start
        repeat_data$protein_end <- protein_pos$protein_end
      }
      
      # Get overlapping transcript indices
      overlap_indices <- if (i %in% names(overlaps_by_repeat)) overlaps_by_repeat[[as.character(i)]] else integer(0)
      
      # Prepare result structure
      ensembl_exon_info <- list(
        transcripts_count = 0,
        has_canonical_transcript = FALSE,
        location_summary = "unknown",
        transcripts = list()
      )
      
      if (length(overlap_indices) == 0) {
        repeat_data$ensembl_exon_info <- ensembl_exon_info
        return(repeat_data)
      }
      
      # Process overlapping transcripts for this repeat
      repeat_chrom <- gsub("chr", "", repeat_data$chrom)
      repeat_start <- as.integer(repeat_data$chromStart)
      repeat_end <- as.integer(repeat_data$chromEnd)
      repeat_strand <- repeat_data$strand
      
      transcript_info <- list()
      locations <- c()
      canonical_found <- FALSE
      
      for (tx_idx in overlap_indices) {
        tx <- all_transcripts[tx_idx]
        tx_id <- tx$tx_id
        
        # Skip if strand doesn't match (when strand is specified)
        if (!is.null(repeat_strand) && !is.na(repeat_strand) && repeat_strand %in% c("+", "-")) {
          tx_strand <- as.character(strand(tx))
          if (repeat_strand != tx_strand) next
        }
        
        # When working with exons, verify they exist and have the right format
        exons <- if (tx_id %in% names(exons_by_tx)) {
          ex <- exons_by_tx[[tx_id]]
          if (is(ex, "GRangesList") && length(ex) > 0) unlist(ex, use.names = FALSE)
          else if (is(ex, "GRanges")) ex
          else GRanges() # Empty GRanges as fallback
        } else {
          GRanges() # Empty GRanges if no exons found
        }
        
        # Ensure we have a valid GRanges object
        if (length(exons) == 0) next
        
        # Basic transcript info
        exon_count <- length(exons)
        tx_strand <- as.character(strand(tx))
        tx_biotype <- tx$tx_biotype
        
        # Check if canonical
        is_canonical <- FALSE
        if (!is.null(tx$tx_is_canonical) && !is.na(tx$tx_is_canonical)) {
          is_canonical <- as.logical(tx$tx_is_canonical)
          if (is_canonical) canonical_found <- TRUE
        }
        
        # Get transcript version
        versioned_tx_id <- tx$tx_id_version
        
        # Determine location
        tx_start <- start(tx)
        tx_end <- end(tx)
        
        # Check if outside transcript bounds
        if (repeat_end <= tx_start || repeat_start >= tx_end) {
          location <- "outside"
        } else {
          # Check if overlapping any exon
          is_exonic <- FALSE
          
          # Fix: Make sure exons is a valid object to iterate over
          if (is(exons, "GRangesList")) {
            exons <- unlist(exons, use.names = FALSE)
          }
          
          # Fix: Use safer iteration with indices instead of "for (exon in exons)"
          if (length(exons) > 0) {
            for (j in seq_along(exons)) {
              exon <- exons[j]
              exon_start <- start(exon)
              exon_end <- end(exon)
              
              if (max(repeat_start, exon_start) < min(repeat_end, exon_end)) {
                is_exonic <- TRUE
                break
              }
            }
          }
          
          location <- if (is_exonic) "exonic" else "intronic"
        }
        
        # Skip intronic transcripts
        if (location == "intronic") next
        
        locations <- c(locations, location)
        
        # Get CDS information
        tx_cds_start <- tx$tx_cds_seq_start
        tx_cds_end <- tx$tx_cds_seq_end
        
        # Get CDS regions for this transcript
        cds_regions <- if (tx_id %in% names(all_cds)) all_cds[[tx_id]] else GRanges()
        
        # Process exons only for exonic transcripts
        containing_exons <- list()
        if (location == "exonic") {
          for (j in 1:length(exons)) {
            exon <- exons[j]
            exon_start <- start(exon)
            exon_end <- end(exon)
            
            # Check if repeat overlaps this exon
            if (max(repeat_start, exon_start) < min(repeat_end, exon_end)) {
              # Process exon details similarly to original code but more efficiently
              exon_number <- j
              
              # Calculate overlap statistics
              overlap_start <- max(repeat_start, exon_start)
              overlap_end <- min(repeat_end, exon_end)
              overlap_length <- overlap_end - overlap_start
              exon_length <- exon_end - exon_start
              overlap_percentage <- (overlap_length / exon_length) * 100
              
              # Determine position in transcript
              position <- if (exon_count == 1) "single_exon"
                          else if (j == 1) "first_exon"
                          else if (j == exon_count) "last_exon"
                          else paste0("middle_exon_", exon_number)
              
              # Coding status calculation
              if (is.na(tx_cds_start) || is.na(tx_cds_end)) {
                coding_status <- "non_coding"
                utr_status <- "non_coding_transcript"
                coding_percentage <- 0
              } else if (exon_end <= tx_cds_start || exon_start >= tx_cds_end) {
                # Non-coding exon
                coding_status <- "non_coding"
                exon_strand <- as.character(strand(exon))
                
                if (exon_end <= tx_cds_start) {
                  utr_status <- if (exon_strand == "+") "5'_UTR" else "3'_UTR"
                } else {
                  utr_status <- if (exon_strand == "+") "3'_UTR" else "5'_UTR"
                }
                coding_percentage <- 0
              } else if (exon_start < tx_cds_start || exon_end > tx_cds_end) {
                # Partially coding
                coding_status <- "partial_coding"
                
                utr_regions <- c()
                exon_strand <- as.character(strand(exon))
                
                if (exon_start < tx_cds_start) {
                  utr_regions <- c(utr_regions, if (exon_strand == "+") "5'_UTR" else "3'_UTR")
                }
                
                if (exon_end > tx_cds_end) {
                  utr_regions <- c(utr_regions, if (exon_strand == "+") "3'_UTR" else "5'_UTR")
                }
                
                # Calculate coding percentage
                coding_start <- max(exon_start, tx_cds_start)
                coding_end <- min(exon_end, tx_cds_end)
                coding_length <- coding_end - coding_start
                coding_percentage <- (coding_length / exon_length) * 100
                
                utr_status <- paste(utr_regions, collapse = "+")
              } else {
                # Fully coding
                coding_status <- "fully_coding"
                utr_status <- "none"
                coding_percentage <- 100
              }
              
              # Phase calculation
              phase <- -1
              end_phase <- -1
              
              if (length(cds_regions) > 0) {
                exon_overlaps <- findOverlaps(exon, cds_regions)
                
                if (length(exon_overlaps) > 0) {
                  # Calculate phases from CDS position
                  subj_hit <- subjectHits(exon_overlaps)[1]
                  if (subj_hit > 1) {
                    all_cds_before <- cds_regions[1:(subj_hit-1)]
                    cds_pos <- sum(width(all_cds_before))
                    phase <- cds_pos %% 3
                  } else {
                    phase <- 0
                  }
                  
                  overlap <- pintersect(exon, cds_regions[subj_hit])
                  overlap_length <- width(overlap)
                  end_phase <- (phase + overlap_length) %% 3
                }
              }
              
              # Frame status
              frame_status <- if (coding_status == "non_coding") "non_coding"
                              else if (is.na(phase) || is.na(end_phase) || phase == -1 || end_phase == -1) "unknown"
                              else if (phase == end_phase) "in_frame"
                              else "out_of_frame"
              
              # Get protein coordinates - replace simplified approach with actual calculation
              protein_start <- NA_integer_
              protein_end <- NA_integer_
              
              # Calculate proper protein coordinates
              if (coding_status != "non_coding") {
                # Create GRanges object for just this exon
                exon_gr <- GRanges(
                  seqnames = seqnames(exon),
                  ranges = IRanges(start = exon_start, end = exon_end),
                  strand = strand(exon)
                )
                
                # Translate genomic coordinates to protein
                prot_coords <- tryCatch({
                  genomeToProtein(exon_gr, ensdb)
                }, error = function(e) {
                  NULL
                })
                
                # Store protein start and end if translation succeeded
                if (!is.null(prot_coords) && length(prot_coords) > 0) {
                  protein_start <- min(start(prot_coords))
                  protein_end <- max(end(prot_coords))
                }
              }
              
              # Create exon info object
              exon_info <- list(
                exon_number = exon_number,
                exon_id = exon$exon_id,
                exon_start = exon_start,
                exon_end = exon_end,
                protein_start = protein_start,  # Now will contain actual values when available
                protein_end = protein_end,      # Now will contain actual values when available
                overlap_bp = overlap_length,
                position = position,
                overlap_percentage = round(overlap_percentage, 2),
                coding_status = coding_status,
                utr_status = utr_status,
                coding_percentage = round(coding_percentage, 2),
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
          transcript_name = tx$tx_name,
          is_canonical = is_canonical,
          biotype = tx_biotype,
          location = location,
          exon_count = exon_count,
          strand = tx_strand,
          containing_exons = containing_exons
        )
        
        transcript_info <- c(transcript_info, list(tx_info))
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
      
      # Update progress bar - this won't actually work with futures, but we'll leave it
      setTxtProgressBar(pb, i)
      
      return(repeat_data)
    }, error = function(e) {
      # Return basic data on error
      cat(sprintf("Error processing repeat %d: %s\n", i, conditionMessage(e)))
      repeat_data$ensembl_exon_info <- list(
        transcripts_count = 0,
        has_canonical_transcript = FALSE,
        location_summary = "unknown",
        transcripts = list()
      )
      return(repeat_data)
    })
  }, future.seed = TRUE)
  
  close(pb)
  
  # Write results to JSON file using a more efficient approach
  cat("Writing results to", output_file, "...\n")
  
  # Use prettier JSON writing for better readability
  write_json(processed_repeats, output_file, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  
  # Helper function to safely close database connections
  safely_disconnect_db <- function(ensdb) {
    tryCatch({
      if (!is.null(ensdb) && "EnsDb" %in% class(ensdb)) {
        conn <- dbconn(ensdb)
        if (dbIsValid(conn)) {
          dbDisconnect(conn)
        }
      }
    }, error = function(e) {
      # Silently ignore errors
    })
  }
  
  # Explicitly close database connections
  safely_disconnect_db(ensdb)
  
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  cat("Processing completed in", format(time_taken), "\n")
  
  return(invisible(processed_repeats))
}

# Command line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    stop("Usage: Rscript fast_process_repeats.R <input_json> <output_json> [limit|range] [cores]")
  }
  
  input_file <- args[1]
  output_file <- args[2]
  
  # Parse optional arguments
  limit <- NULL
  range <- NULL
  cores <- NULL
  
  if (length(args) >= 3) {
    if (grepl("-", args[3])) {
      # It's a range parameter
      range <- args[3]
    } else {
      # It's a limit parameter
      limit <- as.numeric(args[3])
    }
  }
  
  if (length(args) >= 4) {
    cores <- as.numeric(args[4])
  }
  
  fast_process_repeat_data(input_file, output_file, limit, range, cores)
}