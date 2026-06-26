#!/usr/bin/env Rscript

# ============================================================
# 01_build_network_matrices.R
#
# Inputs:
#   data/physeq_objeto_Estero.rds
#   data/interior_mangroves/phyloseq.rds
#   data/interior_impaired/physeq.rds
#
# Outputs:
#   source_data/species_asv/network_matrices/*_counts_common_top200.tsv
#   source_data/species_asv/network_matrices/*_relative_common_top200.tsv
#   source_data/species_asv/network_matrices/*_clr_common_top200.tsv
#   source_data/species_asv/network_matrices/*_metadata_clean.tsv
#
# Algorithmic provenance:
#   Phyloseq object handling follows McMurdie & Holmes (2013).
#   Centered log-ratio transformation follows Aitchison (1986).
#   Code curation/editing algorithm: AI-assisted refactoring with
#   GPT-5.5 Thinking (OpenAI, 2026).
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tibble)
})

OUTPUT_DIR <- "source_data/species_asv/network_matrices"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N <- 200L
PSEUDOCOUNT <- 1

PHYSEQ_ESTEROS <- "data/physeq_objeto_Estero.rds"
PHYSEQ_INTERIOR <- "data/interior_mangroves/phyloseq.rds"
PHYSEQ_IMPAIRED <- "data/interior_impaired/physeq.rds"

read_ps <- function(path) {
  if (!file.exists(path)) stop("Missing input phyloseq object: ", path)
  readRDS(path)
}

count_matrix <- function(ps) {
  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  otu
}

metadata_table <- function(ps) {
  as.data.frame(phyloseq::sample_data(ps), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("sample_id")
}

clr_transform <- function(x, pseudocount = 1) {
  logx <- log(x + pseudocount)
  sweep(logx, 1, rowMeans(logx), "-")
}

write_matrix <- function(x, path) {
  as.data.frame(x, check.names = FALSE) %>%
    tibble::rownames_to_column("sample_id") %>%
    readr::write_tsv(path)
}

build_comparison <- function(ps_a, ps_b, group_a, group_b, comparison_name) {
  mat_a <- count_matrix(ps_a)
  mat_b <- count_matrix(ps_b)

  common_taxa <- intersect(colnames(mat_a), colnames(mat_b))
  if (length(common_taxa) < 10) stop("Too few common taxa for ", comparison_name)

  pooled <- rbind(mat_a[, common_taxa, drop = FALSE], mat_b[, common_taxa, drop = FALSE])
  top_taxa <- names(sort(colSums(pooled), decreasing = TRUE))[seq_len(min(TOP_N, ncol(pooled)))]

  counts <- pooled[, top_taxa, drop = FALSE]
  relative <- sweep(counts, 1, rowSums(counts), "/")
  relative[is.na(relative)] <- 0
  clr <- clr_transform(counts, pseudocount = PSEUDOCOUNT)

  meta <- bind_rows(
    metadata_table(ps_a) %>% mutate(network_group = group_a),
    metadata_table(ps_b) %>% mutate(network_group = group_b)
  )

  write_matrix(counts, file.path(OUTPUT_DIR, paste0(comparison_name, "_counts_common_top200.tsv")))
  write_matrix(relative, file.path(OUTPUT_DIR, paste0(comparison_name, "_relative_common_top200.tsv")))
  write_matrix(clr, file.path(OUTPUT_DIR, paste0(comparison_name, "_clr_common_top200.tsv")))
  readr::write_tsv(meta, file.path(OUTPUT_DIR, paste0(comparison_name, "_metadata_clean.tsv")))
}

ps_estero <- read_ps(PHYSEQ_ESTEROS)
ps_interior <- read_ps(PHYSEQ_INTERIOR)
ps_impaired <- read_ps(PHYSEQ_IMPAIRED)

build_comparison(ps_interior, ps_estero, "Interior_relict", "Coastal_natural", "depth_natural")
build_comparison(ps_interior, ps_impaired, "Interior_relict", "Coastal_impaired", "interior_impaired")

cat("Finished building ASV/species-level network matrices.\n")
