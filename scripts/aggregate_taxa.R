ul <- Sys.getenv("R_LIBS_USER")
if (!nzchar(ul)) ul <- file.path(Sys.getenv("HOME"), "R",
  paste(R.version$platform, paste(R.version$major, R.version$minor, sep="."), "library", sep="-"))
if (!dir.exists(ul)) dir.create(ul, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(ul, .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: aggregate_taxa.R <otu.csv> <tax.csv> <out.csv>")
otu_path <- args[1]; tax_path <- args[2]; out_path <- args[3]

otu <- read.csv(otu_path, check.names = FALSE)
tax <- read.csv(tax_path, check.names = FALSE)
stopifnot("feature_id" %in% colnames(otu), "feature_id" %in% colnames(tax))

if (!requireNamespace("tidyr", quietly=TRUE)) install.packages("tidyr", repos="https://cloud.r-project.org", lib=ul)
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr", repos="https://cloud.r-project.org", lib=ul)
library(tidyr); library(dplyr)

otu_long <- otu %>% pivot_longer(-feature_id, names_to="sample_id", values_to="abundance")
otu_tax <- left_join(otu_long, tax[, c("feature_id","Genus")], by="feature_id")
otu_tax$Genus[is.na(otu_tax$Genus)] <- "Unassigned"

genus <- otu_tax %>%
  group_by(Genus, sample_id) %>%
  summarize(abundance = sum(abundance), .groups="drop") %>%
  pivot_wider(names_from=sample_id, values_from=abundance, values_fill = 0)

write.csv(genus, out_path, row.names = FALSE)
cat("Exported:", out_path, "\n")
