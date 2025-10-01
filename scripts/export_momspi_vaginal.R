ul <- Sys.getenv("R_LIBS_USER")
if (!nzchar(ul)) ul <- file.path(Sys.getenv("HOME"), "R",
  paste(R.version$platform, paste(R.version$major, R.version$minor, sep="."), "library", sep="-"))
if (!dir.exists(ul)) dir.create(ul, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(ul, .libPaths()))

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", lib=ul)
  if (!requireNamespace("HMP2Data", quietly=TRUE)) BiocManager::install("HMP2Data", ask=FALSE, lib=ul, INSTALL_opts="--no-build-vignettes")
  if (!requireNamespace("phyloseq",  quietly=TRUE)) BiocManager::install("phyloseq",  ask=FALSE, lib=ul)
  library(HMP2Data); library(phyloseq)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- ifelse(length(args) >= 1, args[1], "data/processed")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ps   <- momspi16S()
meta <- as.data.frame(sample_data(ps))

# Filter by vaginal site
keep_vals <- c("vagina","vaginal","vaginal swab")
keep <- tolower(trimws(meta$sample_body_site))

# IMPORTANT: sample_names(ps) == meta$file_name
stopifnot("file_name" %in% names(meta))
cand_ids <- meta$file_name[ keep %in% keep_vals ]
cand_ids <- intersect(cand_ids, sample_names(ps))

# Debug
write.csv(data.frame(
  total_meta=nrow(meta),
  total_ps_samples=length(sample_names(ps)),
  matches_by_value=sum(keep %in% keep_vals),
  matches_after_intersect=length(cand_ids),
  chosen_id_col="file_name"
), file.path(outdir, "debug_vaginal_counts.csv"), row.names=FALSE)
if (length(cand_ids) == 0) stop("Zero matches after intersect; check debug_vaginal_counts.csv")

# Prune counts/taxa only; metadata we’ll take from 'meta' directly
ps_vag <- prune_samples(cand_ids, ps)

# OTU matrix (taxa as rows)
otu <- as(otu_table(ps_vag), "matrix")
if (!taxa_are_rows(ps_vag)) otu <- t(otu)
otu <- as.data.frame(otu); otu$feature_id <- rownames(otu)

# Taxonomy
tax <- as.data.frame(tax_table(ps_vag)); tax$feature_id <- rownames(tax)

# Metadata: subset original meta by file_name, not from ps_vag
meta_vag <- meta[match(cand_ids, meta$file_name), , drop=FALSE]
meta_vag$sample_id <- meta_vag$file_name  # explicit key for downstream merges

# Write
write.csv(otu,      file.path(outdir, "vaginal_otu.csv"),      row.names=FALSE)
write.csv(tax,      file.path(outdir, "vaginal_taxonomy.csv"), row.names=FALSE)
write.csv(meta_vag, file.path(outdir, "vaginal_metadata.csv"), row.names=FALSE)
cat("Exported vaginal_otu.csv, vaginal_taxonomy.csv, vaginal_metadata.csv\n")
