#!/usr/bin/env Rscript

# check_uniprot_description.R
# Adds canonical UniProt ID and UniProt description to each entry in a repeat JSON

library(ensembldb)
library(AnnotationHub)
library(jsonlite)
library(dplyr)

# Define output file path
output_file <- "C:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\RTest\\check_uniprot_results2.txt"

# Load the Ensembl database
cat("Loading EnsDb.Hsapiens.v113 from AnnotationHub...\n")
hub <- AnnotationHub()
ensdb <- hub[["AH119325"]]
cat("Ensembl database loaded successfully.\n")

# Path to the single repeat file
json_path <- "C:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\docs\\all_annotated_repeats.json"

if (!file.exists(json_path)) {
  stop(paste("File not found:", json_path))
}

# Read the JSON
cat("Reading", json_path, "...\n")
repeats <- fromJSON(json_path, simplifyDataFrame = FALSE)

if (length(repeats) == 0) {
  stop("JSON file is empty.")
}

# Where to write the updated JSON
updated_json_path <- "C:\\Users\\Okke\\Documents\\GitHub\\Tandem-Repeat-Domain-Database\\RTest\\all_repeats_canonical_uniprot_description.json"

# Helper: normalize UniProt IDs
normalize_uniprot <- function(x) {
  if (is.null(x) || is.na(x)) return(NA_character_)
  x <- trimws(as.character(x))
  x <- strsplit(x, "[;, ]+")[[1]][1]
  x <- trimws(x)
  if (identical(x, "")) return(NA_character_)
  x
}

uniprot_base <- function(x) {
  x <- normalize_uniprot(x)
  if (is.na(x)) return(NA_character_)
  x <- sub("-.*$", "", x)
  x <- sub("\\..*$", "", x)
  x
}

# Helper: insert fields right after uniProtId (and avoid duplicates)
insert_after_uniprot <- function(entry, canonical_uniprot_id, gene_description) {
  insert_fields <- list(
    canonicalUniProtId = canonical_uniprot_id,
    uniProtDescription = gene_description
  )

  # Drop existing keys to avoid duplicates
  entry_names <- names(entry)
  entry <- entry[setdiff(entry_names, names(insert_fields))]

  entry_names <- names(entry)
  insert_after <- match("uniProtId", entry_names)

  if (is.na(insert_after)) {
    c(entry, insert_fields)
  } else {
    c(entry[seq_len(insert_after)], insert_fields, entry[(insert_after + 1):length(entry)])
  }
}

# Process every entry
for (i in seq_along(repeats)) {
  entry <- repeats[[i]]

  gene_symbol <- entry$geneName
  json_uniprot <- normalize_uniprot(entry$uniProtId)
  json_base <- uniprot_base(json_uniprot)

  # Default behavior: canonical is the UniProt ID from JSON (do NOT replace with a different protein)
  canonical_uniprot_id <- json_uniprot
  gene_description <- NA_character_

  # If there is no UniProt ID in JSON, fall back to NA (or later you could use gene_symbol logic if desired)
  if (is.null(canonical_uniprot_id) || is.na(canonical_uniprot_id) || canonical_uniprot_id == "") {
    repeats[[i]] <- insert_after_uniprot(entry, NA_character_, NA_character_)
    next
  }

  # Use EnsDb only to enrich description and (optionally) upgrade to an isoform for the SAME base accession
  avail_cols <- columns(ensdb)

  # 1) Description: prefer matching by UniProtID base using local filtering (UNIPROTID isn't reliably queryable by base)
  if (!is.na(json_base) && "UNIPROTID" %in% avail_cols) {
    cols_small <- intersect(c("UNIPROTID", "DESCRIPTION", "GENENAME", "SYMBOL"), avail_cols)

    # Try: if SYMBOL available, constrain by gene symbol to keep this reasonably small.
    # (But never use SYMBOL to override canonical_uniprot_id.)
    res_small <- NULL
    if (!is.null(gene_symbol) && !is.na(gene_symbol) && gene_symbol != "" && "SYMBOL" %in% avail_cols) {
      tryCatch({
        res_small <- ensembldb::select(
          ensdb,
          keys = gene_symbol,
          keytype = "SYMBOL",
          columns = cols_small
        )
      }, error = function(e) {
        res_small <- NULL
      })
    }

    if (!is.null(res_small) && nrow(res_small) > 0) {
      # Keep only UniProt IDs sharing the same base accession
      res_small$UNIPROT_BASE <- sub("\\..*$", "", sub("-.*$", "", res_small$UNIPROTID))
      hit <- res_small[!is.na(res_small$UNIPROTID) & res_small$UNIPROT_BASE == json_base, , drop = FALSE]

      if (nrow(hit) > 0) {
        # Description
        if ("DESCRIPTION" %in% colnames(hit)) {
          d <- hit$DESCRIPTION[!is.na(hit$DESCRIPTION)][1]
          if (!is.na(d)) gene_description <- trimws(gsub("\\s*\\[.*\\]\\s*$", "", d))
        }

        # Optional: if JSON provided base (e.g., Q9NY46) and EnsDb has isoform(s) for SAME base,
        # prefer an isoform (contains '-') but do not switch to a different base accession.
        if (!grepl("-", canonical_uniprot_id, fixed = TRUE)) {
          iso <- unique(hit$UNIPROTID[!is.na(hit$UNIPROTID) & grepl("-", hit$UNIPROTID, fixed = TRUE)])
          if (length(iso) > 0) canonical_uniprot_id <- iso[[1]]
        }
      }
    }
  }

  repeats[[i]] <- insert_after_uniprot(entry, canonical_uniprot_id, gene_description)
}

cat("\n--- Writing updated JSON ---\n")
cat("Output:", updated_json_path, "\n")
writeLines(toJSON(repeats, auto_unbox = TRUE, pretty = TRUE), updated_json_path)
