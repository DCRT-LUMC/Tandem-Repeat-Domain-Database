# Helper to normalize/compare UniProt IDs (ignore isoform/EnsDb suffixes)
strip_uniprot_suffix <- function(x) {
  # EnsDb sometimes stores values like "Q9C0D5.171" or isoforms like "Q9C0D5-1"
  sub("[.-].*$", "", x)
}

# Wherever the canonical UniProt ID is decided, change to:
# - Trust the JSON UniProt ID as canonical.
# - Only use EnsDb UniProt IDs for *validation/diagnostics*, not for substitution.

# json_uniprot_id <- <value read from JSON>
# symbol <- <gene symbol>
# ensdb_uniprot_ids <- <vector of UniProt IDs from EnsDb for this symbol>

# Replace the old fallback block with this logic:
canonical_uniprot_id <- json_uniprot_id

if (!is.null(ensdb_uniprot_ids) && length(ensdb_uniprot_ids) > 0 && !is.na(json_uniprot_id) && nzchar(json_uniprot_id)) {
  json_norm <- strip_uniprot_suffix(json_uniprot_id)
  ensdb_norm <- strip_uniprot_suffix(ensdb_uniprot_ids)

  if (!(json_norm %in% ensdb_norm)) {
    warning(sprintf(
      "The UniProt ID from JSON (%s) was NOT found in EnsDb for symbol %s. Keeping JSON UniProt as canonical.\nAvailable UniProt IDs for %s are:\n%s",
      json_uniprot_id,
      symbol,
      symbol,
      paste0(capture.output(print(ensdb_uniprot_ids)), collapse = "\n")
    ))
  }
}

# IMPORTANT: ensure any downstream code uses `canonical_uniprot_id` (the JSON value)
# and does not overwrite it with EnsDb-derived UniProt IDs.