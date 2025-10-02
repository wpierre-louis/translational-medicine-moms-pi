#!/usr/bin/env Rscript
# HMP2 MOMS-PI 16S → CSVs (vaginal subset), no phyloseq

args <- commandArgs(trailingOnly = TRUE)
outdir <- ifelse(length(args) >= 1, args[1], "data/processed")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  if (!requireNamespace("HMP2Data", quietly = TRUE))
    stop("Install with: BiocManager::install('HMP2Data')")
})

# Load the package-shipped data objects directly (avoids memory spikes)
data("momspi16S_mtx",  package = "HMP2Data", envir = environment())
data("momspi16S_samp", package = "HMP2Data", envir = environment())
data("momspi16S_tax",  package = "HMP2Data", envir = environment())

mtx  <- get("momspi16S_mtx")                                 # taxa x samples
samp <- as.data.frame(get("momspi16S_samp"), stringsAsFactors = FALSE)
tax  <- as.data.frame(get("momspi16S_tax"),  stringsAsFactors = FALSE)

# Requirements for filtering and IDs
stopifnot("file_name" %in% names(samp), "sample_body_site" %in% names(samp))

# NA-safe vaginal filter
keep_vals <- c("vagina", "vaginal", "vaginal swab")
sb <- tolower(trimws(as.character(samp$sample_body_site)))
flag <- !is.na(sb) & (sb %in% keep_vals)

# Candidate sample IDs from metadata
cand_ids <- samp$file_name[flag]

# Ensure matrix orientation is samples in columns; flip if needed
if (!all(cand_ids %in% colnames(mtx))) {
  if (all(cand_ids %in% rownames(mtx))) {
    mtx <- t(mtx)
  } else {
    write.csv(
      as.data.frame(sort(table(replace(sb, is.na(sb), "<NA>")), decreasing = TRUE)),
      file.path(outdir, "audit_body_site_counts.csv"),
      row.names = FALSE
    )
    stop("Candidate sample IDs did not match count matrix dimnames.")
  }
}
# Final intersect
cand_ids <- intersect(cand_ids, colnames(mtx))
if (!length(cand_ids)) stop("Zero vaginal samples after intersecting with count matrix.")

# Subset counts and drop zero-abundance taxa
mtx_sub <- mtx[, cand_ids, drop = FALSE]
keep_taxa <- rowSums(mtx_sub != 0) > 0
mtx_sub <- mtx_sub[keep_taxa, , drop = FALSE]

# Subset taxonomy to remaining taxa and add feature key
tax_sub <- tax[rownames(mtx_sub), , drop = FALSE]
tax_sub$feature_id <- rownames(tax_sub)

# Build vaginal metadata for those samples
samp_vag <- samp[match(cand_ids, samp$file_name), , drop = FALSE]
samp_vag$sample_id <- samp_vag$file_name

# Race and gender mapping (HMP2 columns)
if (!("subject_gender" %in% names(samp_vag))) samp_vag$subject_gender <- NA_character_
if (!("subject_race"   %in% names(samp_vag))) samp_vag$subject_race   <- NA_character_

norm_lower <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "n/a", "null", "none", "unknown")] <- NA
  x
}
map_vec <- function(x, map) { out <- x; if (length(map)) for (k in names(map)) out[x == k] <- map[[k]]; out }

gender_raw <- norm_lower(samp_vag$subject_gender)
race_raw   <- norm_lower(samp_vag$subject_race)

gender_map <- c("female" = "Female", "male" = "Male")
race_map <- c(
  "caucasian"                        = "White / Caucasian",
  "african american"                 = "Black / African American",
  "asian"                            = "Asian",
  "hispanic or latino"               = "Hispanic / Latino",
  "ethnic other"                     = "Other",
  "american indian or alaska native" = "American Indian / Alaska Native"
)

samp_vag$subject_gender <- map_vec(gender_raw, gender_map)
samp_vag$subject_race   <- map_vec(race_raw,   race_map)

# Keep only the requested metadata columns (drop URLs and everything else)
wanted <- c(
  "file_name", "sample_id", "subject_id", "visit_number",
  "subject_gender", "subject_race",
  "sample_body_site", "study_full_name", "project_name"
)
have <- intersect(wanted, names(samp_vag))
missing_cols <- setdiff(wanted, have)
if (length(missing_cols)) message("Missing metadata columns (omitted): ", paste(missing_cols, collapse = ", "))
samp_vag <- samp_vag[, have, drop = FALSE]

# Prepare outputs
otu_out <- as.data.frame(mtx_sub, stringsAsFactors = FALSE)
otu_out$feature_id <- rownames(otu_out)

# Audits
aud <- function(v) as.data.frame(sort(table(v, useNA = "ifany"), decreasing = TRUE))
write.csv(aud(gender_raw), file.path(outdir, "audit_subject_gender_counts.csv"), row.names = FALSE)
write.csv(aud(race_raw),   file.path(outdir, "audit_subject_race_counts.csv"),   row.names = FALSE)
write.csv(data.frame(
  total_samples_all       = nrow(samp),
  vaginal_flag_true       = sum(flag, na.rm = TRUE),
  vaginal_after_intersect = length(cand_ids),
  taxa_after_filter       = nrow(otu_out)
), file.path(outdir, "debug_vaginal_counts.csv"), row.names = FALSE)

# Write CSVs
write.csv(otu_out,  file.path(outdir, "vaginal_otu.csv"),      row.names = FALSE)
write.csv(tax_sub,  file.path(outdir, "vaginal_taxonomy.csv"), row.names = FALSE)
write.csv(samp_vag, file.path(outdir, "vaginal_metadata.csv"), row.names = FALSE)
cat("Exported:\n",
    file.path(outdir, "vaginal_otu.csv"), "\n",
    file.path(outdir, "vaginal_taxonomy.csv"), "\n",
    file.path(outdir, "vaginal_metadata.csv"), "\n", sep = "")