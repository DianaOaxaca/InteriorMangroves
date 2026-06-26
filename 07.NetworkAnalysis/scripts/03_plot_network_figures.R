#!/usr/bin/env Rscript

# ============================================================
# 03_plot_network_figures_pdf_only.R
#
# Inputs:
#   source_data/genus_level/19_final_figure_tables/19_source_data_main_contrasts_long.tsv
#   source_data/genus_level/19_final_figure_tables/19_source_data_robustness_curves.tsv
#
# Outputs:
#   figures/genus_level/genus_level_main_contrasts.pdf
#   figures/genus_level/genus_level_robustness_curves.pdf
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(scales)
})

INPUT_DIR <- "source_data/genus_level/19_final_figure_tables"
OUTPUT_DIR <- "figures/genus_level"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

group_colors <- c(
  Interior_relict = "#1B6B3A",
  Coastal_natural = "#2C7FB8",
  Coastal_impaired = "#D95F02"
)

group_labels <- c(
  Interior_relict = "Interior relict",
  Coastal_natural = "Coastal natural",
  Coastal_impaired = "Coastal impaired"
)

theme_manuscript <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black", size = base_size),
      axis.title = element_text(size = 9.5),
      axis.text = element_text(size = 8.5, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 8.5),
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(size = 9),
      plot.margin = margin(5, 6, 5, 6)
    )
}

plot_df <- readr::read_tsv(file.path(INPUT_DIR, "19_source_data_main_contrasts_long.tsv"), show_col_types = FALSE)
curves <- readr::read_tsv(file.path(INPUT_DIR, "19_source_data_robustness_curves.tsv"), show_col_types = FALSE)

p_main <- ggplot(plot_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.30, color = "black") +
  geom_point(aes(color = group), position = position_jitter(width = 0.12), size = 0.35, alpha = 0.16, stroke = 0) +
  facet_grid(metric_label ~ dataset_label, scales = "free_y") +
  scale_fill_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  scale_color_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  labs(x = NULL, y = "Bootstrap network metric") +
  theme_manuscript(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "genus_level_main_contrasts.pdf"), p_main, width = 8.2, height = 8.8, bg = "white", useDingbats = FALSE)

p_curves <- ggplot(curves, aes(x = removal_fraction, y = median_largest_component_prop, color = group, fill = group)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 0.75) +
  facet_grid(mode_label ~ dataset_label) +
  scale_color_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  scale_fill_manual(values = group_colors, labels = group_labels, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 0.9, by = 0.15)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
  labs(x = "Fraction of nodes removed", y = "Remaining largest component") +
  theme_manuscript(9)

ggsave(file.path(OUTPUT_DIR, "genus_level_robustness_curves.pdf"), p_curves, width = 8.0, height = 5.8, bg = "white", useDingbats = FALSE)

cat("Finished PDF plotting.\n")
