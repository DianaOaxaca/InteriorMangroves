#!/usr/bin/env Rscript

# ============================================================
# 04_reviewer_controls_depth_interior_genus.R
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

OUTPUT_DIR <- "results/reviewer_controls"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

files <- list(
  depth_profile = "source_data/species_asv/depth_profile/09_ridge_depth_trend_stats.tsv",
  interior_only = "source_data/interior_only/14_interior_location_network_tables/14_pairwise_location_contrasts.tsv",
  genus_level = "source_data/genus_level/18_genus_ridge_network_tables/18_pairwise_group_contrasts.tsv"
)

for (nm in names(files)) {
  if (!file.exists(files[[nm]])) stop("Missing required input: ", files[[nm]])
}

write_tsv(read_tsv(files$depth_profile, show_col_types = FALSE), file.path(OUTPUT_DIR, "depth_profile_trend_stats.tsv"))
write_tsv(read_tsv(files$interior_only, show_col_types = FALSE), file.path(OUTPUT_DIR, "interior_only_pairwise_summary.tsv"))
write_tsv(read_tsv(files$genus_level, show_col_types = FALSE), file.path(OUTPUT_DIR, "genus_level_pairwise_summary.tsv"))

cat("Reviewer-control summaries written to: ", OUTPUT_DIR, "\n", sep = "")
