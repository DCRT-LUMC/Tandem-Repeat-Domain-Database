#!/usr/bin/env Rscript

# ensdb_lookup_uniprot.R
# Lookup a UniProt ID in EnsDb.Hsapiens.v113 and print all matching rows/columns.

suppressPackageStartupMessages({
  library(ensembldb)
  library(AnnotationHub)
})

# UniProt ID to lookup (edit this)
up_in <- "Q9NY46"  # e.g. "Q9NY46", "Q9NY46-1", "Q9NY46.204"

cat("Loading EnsDb.Hsapiens.v113 from AnnotationHub...\n")
hub <- AnnotationHub()
ensdb <- hub[["AH119325"]]
cat("Loaded.\n\n")

avail_cols <- columns(ensdb)
cat("EnsDb columns (", length(avail_cols), "):\n", sep = "")
cat(paste(sort(avail_cols), collapse = ", "), "\n\n", sep = "")

if (!"UNIPROTID" %in% avail_cols) {
  stop("This EnsDb does not expose a UNIPROTID column.")
}

# Normalize input
up_in <- trimws(up_in)
up_first <- strsplit(up_in, "[;, ]+")[[1]][1]
up_base_iso <- sub("-.*$", "", up_first)
up_base_dot <- sub("\\..*$", "", up_base_iso)
base <- up_base_dot

cat("Input: ", up_in, "\n", sep = "")
cat("Base accession used for prefix search: ", base, "\n\n", sep = "")

# 1) Exact match attempts (EnsDb requires exact keys)
keys_to_try <- unique(c(up_first, up_base_iso, up_base_dot))
cat("Exact keys to try: ", paste(keys_to_try, collapse = ", "), "\n\n", sep = "")

for (k in keys_to_try) {
  cat("=== Exact query UNIPROTID = ", k, " ===\n", sep = "")
  res <- NULL
  tryCatch({
    res <- ensembldb::select(
      ensdb,
      keys = k,
      keytype = "UNIPROTID",
      columns = avail_cols
    )
  }, error = function(e) {
    cat("ERROR: ", conditionMessage(e), "\n\n", sep = "")
  })

  if (is.null(res) || nrow(res) == 0) {
    cat("No rows returned.\n\n")
  } else {
    cat("Rows: ", nrow(res), "\n", sep = "")
    print(res)
    cat("\n")
  }
}

# 2) Prefix match: fetch a small mapping table and filter locally
cat("=== Prefix search for UNIPROTID starting with '", base, "' ===\n", sep = "")
cols_small <- intersect(c("UNIPROTID", "UNIPROTDB", "UNIPROTMAPPINGTYPE", "GENEID", "SYMBOL", "TXID", "PROTEINID"), avail_cols)

# Pull mappings by asking for rows keyed by UNIPROTID is not possible without exact keys,
# so we instead pull a broad mapping table (may be large).
# To constrain, we only request the minimal columns.
all_map <- NULL
tryCatch({
  all_map <- ensembldb::select(
    ensdb,
    keys = keys(ensdb, keytype = "GENEID"),
    keytype = "GENEID",
    columns = cols_small
  )
}, error = function(e) {
  cat("ERROR building mapping table: ", conditionMessage(e), "\n", sep = "")
})

if (is.null(all_map) || nrow(all_map) == 0) {
  cat("Could not build mapping table for prefix search.\n")
} else {
  hit <- all_map[!is.na(all_map$UNIPROTID) & startsWith(all_map$UNIPROTID, base), , drop = FALSE]
  if (nrow(hit) == 0) {
    cat("No UNIPROTID values starting with '", base, "' found in EnsDb (for the requested mapping columns).\n", sep = "")
  } else {
    cat("Matches: ", nrow(hit), "\n", sep = "")
    print(unique(hit))
  }
}
