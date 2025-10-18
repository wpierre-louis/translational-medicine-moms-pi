#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(HMP2Data)
  library(phyloseq)
})

# --- Setup output path ---
outdir <- "data/processed"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load dataset ---
ps <- momspi16S()
meta <- as.data.frame(sample_data(ps))

# --- Filter to vaginal samples only ---
keep_vals <- c("vagina","vaginal","vaginal swab")
keep_idx <- tolower(trimws(meta$sample_body_site)) %in% keep_vals
cand_ids <- meta$file_name[keep_idx]
cand_ids <- intersect(cand_ids, sample_names(ps))
stopifnot(length(cand_ids) > 0)
ps_vag <- prune_samples(cand_ids, ps)

# --- Calculate alpha diversity metrics ---
alpha <- estimate_richness(ps_vag, measures = c("Shannon", "Simpson"))
alpha$sample_id <- rownames(alpha)

# --- Merge with sample metadata ---
meta_vag <- meta[match(alpha$sample_id, meta$file_name), 
                 c("file_name","subject_id","visit_number","subject_gender","subject_race",
                   "sample_body_site","study_full_name","project_name"),
                 drop = FALSE]
meta_vag$sample_id <- meta_vag$file_name
alpha_merged <- merge(alpha, meta_vag, by = "sample_id", all.x = TRUE)

# --- Write outputs ---
outfile <- file.path(outdir, "tableau_alpha_diversity.csv")
write.csv(alpha_merged, outfile, row.names = FALSE)
cat(" Wrote:", outfile, "with", nrow(alpha_merged), "rows\n")