#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(HMP2Data)
  library(phyloseq)
})

# Setup output path
outdir <- "data/processed"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Load HMP2 MomSPI dataset 
ps <- momspi16S()

# Calculate alpha diversity metrics
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson"))

# Clean sample IDs 
alpha$sample_id <- tolower(trimws(rownames(alpha)))

# Load pre-exported vaginal metadata 
meta_vag <- read.csv(file.path(outdir, "vaginal_metadata.csv"), stringsAsFactors = FALSE)

# --- Normalize IDs for joining ---
meta_vag$file_name <- tolower(trimws(meta_vag$file_name))
meta_vag$sample_id <- tolower(trimws(meta_vag$sample_id))

# Debug overlap
ov <- length(intersect(alpha$sample_id, meta_vag$sample_id))
cat("Overlap detected:", ov, "of", nrow(alpha), "samples\n")

# Merge alpha diversity results with demographics 
merged <- merge(
  alpha,
  meta_vag[, c(
    "sample_id","subject_id","visit_number","subject_gender",
    "subject_race","sample_body_site","study_full_name","project_name"
  )],
  by = "sample_id",
  all.x = TRUE
)
merged <- merged[merged$sample_id %in% meta_vag$sample_id, ]

# Export 
out_csv <- file.path(outdir, "tableau_alpha_diversity.csv")
write.csv(merged, out_csv, row.names = FALSE)
cat("Wrote merged alpha diversity with", nrow(merged), "rows\n")