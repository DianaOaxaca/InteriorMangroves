#!/usr/bin/env Rscript

# ============================================================
# 02_infer_ridge_networks.R
# Inputs:
#   source_data/species_asv/network_matrices/*_clr_common_top200.tsv
#   source_data/species_asv/network_matrices/*_metadata_clean.tsv
#
# Outputs:
#   results/ridge_networks/species_asv/network_metrics_bootstrap.tsv
#   results/ridge_networks/species_asv/network_metrics_summary.tsv
#   results/ridge_networks/species_asv/edge_detection_frequency.tsv
#
# Algorithmic provenance:
#   Ridge partial correlations follow Schäfer & Strimmer (2005).
#   Bootstrap follows Efron & Tibshirani (1993).
#   Louvain modularity follows Blondel et al. (2008).
#   Network metrics use igraph following Csardi & Nepusz (2006).
#   Code curation/editing algorithm: AI-assisted refactoring with
#   GPT-5.5 Thinking (OpenAI, 2026).
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(igraph)
})

INPUT_DIR <- "source_data/species_asv/network_matrices"
OUTPUT_DIR <- "results/ridge_networks/species_asv"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

COMPARISONS <- c("depth_natural", "interior_impaired")
N_BOOT <- 500L
RIDGE_LAMBDA <- 0.10
EDGE_DENSITY <- 0.02
set.seed(123)

read_matrix <- function(path) {
  x <- readr::read_tsv(path, show_col_types = FALSE)
  ids <- x[[1]]
  mat <- as.matrix(x[, -1, drop = FALSE])
  rownames(mat) <- ids
  storage.mode(mat) <- "numeric"
  mat
}

ridge_partial_cor <- function(x, lambda = 0.10) {
  x <- x[, apply(x, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  x <- scale(x, center = TRUE, scale = TRUE)
  x[is.na(x)] <- 0
  s <- cov(x)
  s[is.na(s)] <- 0
  precision <- solve(s + diag(lambda, ncol(s)))
  d <- diag(precision)
  pc <- -precision / sqrt(outer(d, d))
  diag(pc) <- 0
  pc[is.na(pc) | is.infinite(pc)] <- 0
  pmax(pmin(pc, 0.999), -0.999)
}

make_edges <- function(pc, density = 0.02) {
  taxa <- colnames(pc)
  ij <- which(upper.tri(pc), arr.ind = TRUE)
  edges <- tibble(
    from = taxa[ij[, 1]],
    to = taxa[ij[, 2]],
    weight = pc[ij],
    abs_weight = abs(pc[ij])
  ) %>% arrange(desc(abs_weight), from, to)

  n_edges <- max(1L, floor(density * nrow(edges)))
  dplyr::slice_head(edges, n = n_edges)
}

graph_metrics <- function(edges, nodes) {
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = tibble(name = nodes))
  comps <- igraph::components(g)
  w <- abs(igraph::E(g)$weight)
  w[w <= 0 | is.na(w)] <- 1e-6
  comm <- igraph::cluster_louvain(g, weights = w)
  a <- as.matrix(igraph::as_adjacency_matrix(g, attr = NULL, sparse = FALSE))
  eig <- eigen(a, symmetric = TRUE, only.values = TRUE)$values

  tibble(
    n_nodes = igraph::vcount(g),
    n_edges = igraph::ecount(g),
    average_degree = mean(igraph::degree(g)),
    largest_component_prop = max(comps$csize) / igraph::vcount(g),
    modularity_louvain = igraph::modularity(comm),
    natural_connectivity = log(mean(exp(eig - max(eig)))) + max(eig),
    transitivity = igraph::transitivity(g, type = "global", isolates = "zero"),
    prop_positive = mean(igraph::E(g)$weight > 0)
  )
}

run_comparison <- function(comparison) {
  mat <- read_matrix(file.path(INPUT_DIR, paste0(comparison, "_clr_common_top200.tsv")))
  meta <- readr::read_tsv(file.path(INPUT_DIR, paste0(comparison, "_metadata_clean.tsv")), show_col_types = FALSE)
  id_col <- if ("sample_id" %in% names(meta)) "sample_id" else names(meta)[1]

  meta <- meta %>% filter(.data[[id_col]] %in% rownames(mat))
  mat <- mat[meta[[id_col]], , drop = FALSE]

  groups <- sort(unique(meta$network_group))
  n_sub <- min(table(meta$network_group))

  metrics <- list()
  edges_all <- list()

  for (b in seq_len(N_BOOT)) {
    for (grp in groups) {
      pool <- which(meta$network_group == grp)
      idx <- sample(pool, n_sub, replace = TRUE)
      pc <- ridge_partial_cor(mat[idx, , drop = FALSE], RIDGE_LAMBDA)
      edges <- make_edges(pc, EDGE_DENSITY)

      metrics[[paste(comparison, grp, b, sep = "__")]] <- graph_metrics(edges, colnames(pc)) %>%
        mutate(comparison = comparison, group = grp, bootstrap_id = b)

      edges_all[[paste(comparison, grp, b, sep = "__")]] <- edges %>%
        mutate(comparison = comparison, group = grp, bootstrap_id = b)
    }
  }

  list(metrics = bind_rows(metrics), edges = bind_rows(edges_all))
}

res <- lapply(COMPARISONS, run_comparison)
metrics <- bind_rows(lapply(res, `[[`, "metrics"))
edges <- bind_rows(lapply(res, `[[`, "edges"))

metrics_summary <- metrics %>%
  group_by(comparison, group) %>%
  summarise(
    across(
      c(n_nodes, n_edges, average_degree, largest_component_prop, modularity_louvain, natural_connectivity, transitivity, prop_positive),
      list(median = median, q025 = ~quantile(.x, 0.025), q975 = ~quantile(.x, 0.975)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

edge_freq <- edges %>%
  group_by(comparison, group, from, to) %>%
  summarise(
    detection_frequency = n() / N_BOOT,
    median_weight = median(weight),
    median_abs_weight = median(abs_weight),
    .groups = "drop"
  )

write_tsv(metrics, file.path(OUTPUT_DIR, "network_metrics_bootstrap.tsv"))
write_tsv(metrics_summary, file.path(OUTPUT_DIR, "network_metrics_summary.tsv"))
write_tsv(edge_freq, file.path(OUTPUT_DIR, "edge_detection_frequency.tsv"))

cat("Finished ridge network inference.\n")
