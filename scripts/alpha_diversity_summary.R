#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# --- Load alpha diversity dataset ---
df <- read_csv("data/processed/tableau_alpha_diversity.csv", show_col_types = FALSE)

# --- Ensure numeric and categorical fields ---
df <- df %>%
  mutate(
    visit_number = as.numeric(visit_number),
    Shannon = as.numeric(Shannon),
    Simpson = as.numeric(Simpson)
  )

# --- Summarize alpha diversity by visit number ---
summary_by_visit <- df %>%
  filter(!is.na(visit_number)) %>%
  group_by(visit_number) %>%
  summarise(
    mean_Shannon = mean(Shannon, na.rm = TRUE),
    sd_Shannon = sd(Shannon, na.rm = TRUE),
    mean_Simpson = mean(Simpson, na.rm = TRUE),
    sd_Simpson = sd(Simpson, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(visit_number)

# --- Write output ---
out_path <- "data/processed/tableau_alpha_by_visit.csv"
write_csv(summary_by_visit, out_path)

cat("Wrote:", out_path, "| Rows:", nrow(summary_by_visit), "\n")
