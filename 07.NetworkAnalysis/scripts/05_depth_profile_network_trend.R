#!/usr/bin/env Rscript

# ============================================================
# 05_depth_profile_network_trend.R
#
# Repository inputs:
#   source_data/species_asv/02_network_metrics_bootstrap.tsv
#   source_data/species_asv/depth_profile/09_ridge_depth_trend_stats.tsv
#   source_data/species_asv/depth_profile/09_ridge_source_data_depth_profile.tsv
#
# Outputs:
#   results/depth_profile/depth_trend_stats_recomputed.tsv
#   results/depth_profile/depth_profile_metric_summary_recomputed.tsv
#   results/depth_profile/depth_profile_source_data_from_bootstrap.tsv
#   results/depth_profile/depth_trend_stats_from_publication_tables.tsv
#   results/depth_profile/depth_profile_source_data_from_publication_tables.tsv
#
# Algorithmic provenance:
#   - Upstream networks were inferred as ridge-regularized partial
#     correlations on CLR-transformed common top-200 ASV/species tables.
#   - Depth-profile testing uses Spearman rank correlation between
#     ordinal sediment depth horizon and bootstrap network metrics.
#   - Kruskal-Wallis tests summarize distributional differences among
#     depth horizons.
#   - Bootstrap resampling follows Efron & Tibshirani (1993).
#   - CLR transformation follows Aitchison (1986).
#   - Ridge partial-correlation logic follows Schäfer & Strimmer (2005).
#   - Network metrics were computed with igraph following Csardi & Nepusz (2006).
#   - Code curation/editing algorithm: AI-assisted refactoring with
#     GPT-5.5 Thinking (OpenAI, 2026).
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

OUTPUT_DIR <- "results/depth_profile"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

BOOT_FILE <- "source_data/species_asv/02_network_metrics_bootstrap.tsv"
PUBLISHED_STATS_FILE <- "source_data/species_asv/depth_profile/09_ridge_depth_trend_stats.tsv"
PUBLISHED_SOURCE_FILE <- "source_data/species_asv/depth_profile/09_ridge_source_data_depth_profile.tsv"

depth_order <- c("Depth_0_15", "Depth_16_30", "Depth_31_45", "Depth_50_75")
depth_midpoint <- c(Depth_0_15 = 7.5, Depth_16_30 = 23, Depth_31_45 = 38, Depth_50_75 = 62.5)

metrics_to_test <- c(
  "largest_component_prop",
  "modularity_louvain",
  "natural_connectivity",
  "transitivity",
  "robustness_auc_random",
  "robustness_auc_targeted"
)

check_file <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

check_file(BOOT_FILE)
check_file(PUBLISHED_STATS_FILE)
check_file(PUBLISHED_SOURCE_FILE)

boot <- readr::read_tsv(BOOT_FILE, show_col_types = FALSE)

required_cols <- c("dataset", "group")
missing_cols <- setdiff(required_cols, names(boot))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

depth_df <- boot %>%
  filter(dataset == "interior_depth", group %in% depth_order) %>%
  mutate(
    depth_group = factor(group, levels = depth_order),
    depth_rank = as.numeric(depth_group),
    depth_mid_cm = unname(depth_midpoint[as.character(group)])
  )

if (nrow(depth_df) == 0) stop("No rows found for dataset == interior_depth and expected depth groups.")

available_metrics <- intersect(metrics_to_test, names(depth_df))
if (length(available_metrics) == 0) stop("None of the expected network metrics were found.")

spearman_one <- function(metric_name) {
  x <- depth_df$depth_rank
  y <- depth_df[[metric_name]]
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 4) {
    return(tibble(
      metric = metric_name, n = sum(ok),
      spearman_rho = NA_real_, spearman_p = NA_real_,
      kruskal_chisq = NA_real_, kruskal_df = NA_real_, kruskal_p = NA_real_,
      trend_direction = "insufficient_data"
    ))
  }

  sp <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  kw <- suppressWarnings(kruskal.test(y[ok] ~ depth_df$depth_group[ok]))
  rho <- unname(sp$estimate)

  tibble(
    metric = metric_name,
    n = sum(ok),
    spearman_rho = rho,
    spearman_p = sp$p.value,
    kruskal_chisq = unname(kw$statistic),
    kruskal_df = unname(kw$parameter),
    kruskal_p = kw$p.value,
    trend_direction = case_when(
      rho < 0 ~ "decreases_with_depth",
      rho > 0 ~ "increases_with_depth",
      TRUE ~ "no_monotonic_trend"
    )
  )
}

trend_stats <- bind_rows(lapply(available_metrics, spearman_one))

depth_summary <- depth_df %>%
  group_by(depth_group, depth_mid_cm) %>%
  summarise(
    n_boot = n(),
    across(
      all_of(available_metrics),
      list(
        median = ~median(.x, na.rm = TRUE),
        q025 = ~quantile(.x, 0.025, na.rm = TRUE),
        q975 = ~quantile(.x, 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

published_stats <- readr::read_tsv(PUBLISHED_STATS_FILE, show_col_types = FALSE)
published_source <- readr::read_tsv(PUBLISHED_SOURCE_FILE, show_col_types = FALSE)

readr::write_tsv(trend_stats, file.path(OUTPUT_DIR, "depth_trend_stats_recomputed.tsv"))
readr::write_tsv(depth_summary, file.path(OUTPUT_DIR, "depth_profile_metric_summary_recomputed.tsv"))
readr::write_tsv(depth_df, file.path(OUTPUT_DIR, "depth_profile_source_data_from_bootstrap.tsv"))
readr::write_tsv(published_stats, file.path(OUTPUT_DIR, "depth_trend_stats_from_publication_tables.tsv"))
readr::write_tsv(published_source, file.path(OUTPUT_DIR, "depth_profile_source_data_from_publication_tables.tsv"))

cat("\nDepth-profile network trend analysis complete.\n")
cat("Input bootstrap metrics: ", BOOT_FILE, "\n", sep = "")
cat("Output directory: ", OUTPUT_DIR, "\n\n", sep = "")
print(trend_stats, n = Inf, width = Inf)
