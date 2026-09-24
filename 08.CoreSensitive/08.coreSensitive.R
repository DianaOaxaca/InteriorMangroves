## =============================================================================
## 08.coreSensitive.R
## Sensitivity of the core microbiome to the prevalence threshold
## (Reviewer 2, comment 6 — Communications Biology COMMSBIO-26-2302-T)
##
## Reproduces the operational core definition used in 02.SurfaceComparison and
## 03.4.Coremicrobiome  (feature present in >= X% of samples within each group;
## binary detection; unrarefied 97%-clustered features) and sweeps the threshold
## across 30 / 50 / 70 / 90 %. Writes Supplementary Table S18.
##
## At >=50% it reproduces the values already reported in the manuscript:
##   Surface: Fossil Lagoon 1427|San Pedro River 1419|Terminos 70| Celestun 1230
##            global surface core = 52
##   Depth  : Fossil Lagoon 1481 | San Pedro Martir 1695 | Terminos impaired 123
##
## Base R only (plus qiime2R to read the .qza); no tidyverse dependency.
## =============================================================================

library(qiime2R)

## ---- 0. Inputs (same files read from data/ by the other scripts) ------------
table_qza      <- "data/cluster_table_filter_freq218_emcEstero.qza"
meta_file      <- "data/metadata14.tsv"
out_tsv        <- "08.CoreSensitive/SupplementaryTable_S18_core_sensitivity.tsv"
thresholds     <- c(0.30, 0.50, 0.70, 0.90)
surface_depths <- c("0-15", "0.15", "5")      # surface horizon per dataset

## ---- 1. Load feature table + metadata ---------------------------------------
otu <- read_qza(table_qza)$data               # features x samples (counts)
otu <- as.matrix(otu)

meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(meta) && grepl("^#", meta[[1]][1])) meta <- meta[-1, ]
rownames(meta) <- meta$id

common <- intersect(colnames(otu), meta$id)   # align samples present in both
otu    <- otu[, common, drop = FALSE]
meta   <- meta[common, ]

## Study_zone -> manuscript system labels
zone_to_system <- c("Laguna Cacahuate" = "Fossil Lagoon",
                    "Rio San Pedro"    = "San Pedro River",
                    "Estero_pargo"     = "Terminos Lagoon",
                    "Celestún"         = "Celestun")
meta$Mangrove_system <- zone_to_system[meta$Study_zone]

## ---- 2. Core helpers (identical definition to the main analyses) ------------
core_ids  <- function(ids, thr) {
  prev <- rowSums(otu[, ids, drop = FALSE] > 0) / length(ids)
  names(prev)[prev >= thr]
}
core_size <- function(ids, thr) length(core_ids(ids, thr))

pct <- paste0(">=", thresholds * 100, "%")

## ---- 3. Panel A — surface cross-system (Fringe, flood, surface horizon) -----
surf_sys <- c("Fossil Lagoon", "San Pedro River", "Terminos Lagoon", "Celestun")
surf_ids <- function(sys) meta$id[meta$Mangrove_system == sys &
                                  meta$Ecological_type == "Fringe" &
                                  meta$season          == "flood"  &
                                  meta$depth %in% surface_depths]

A <- sapply(surf_sys, function(sys) sapply(thresholds,
                      function(t) core_size(surf_ids(sys), t)))
rownames(A) <- pct                                   # thresholds x systems

## global surface core (core in ALL four systems)
global_core <- sapply(thresholds, function(t)
  length(Reduce(intersect, lapply(surf_sys,
                          function(sys) core_ids(surf_ids(sys), t)))))

panelA <- data.frame(Prevalence = pct, A, `Global core (all 4)` = global_core,
                     check.names = FALSE, row.names = NULL)

## ---- 4. Panel B — depth profiles (flood; all horizons) ----------------------
depth_ids <- function(zone, etype = NA) {
  keep <- meta$Study_zone == zone & meta$season == "flood"
  if (!is.na(etype)) keep <- keep & meta$Ecological_type == etype
  meta$id[keep]
}
B_groups <- list(
  "Fossil Lagoon (inland)"       = depth_ids("Laguna Cacahuate"),
  "San Pedro Martir R. (inland)" = depth_ids("Rio San Pedro"),
  "Terminos impaired (coastal)"  = depth_ids("Estero_pargo", "Impaired"))
B <- sapply(B_groups, function(ids) sapply(thresholds, function(t) core_size(ids, t)))
rownames(B) <- pct
panelB <- data.frame(Prevalence = pct, B, check.names = FALSE, row.names = NULL)

## ---- 5. Sample sizes (for the legend) ---------------------------------------
n_surf  <- sapply(surf_sys,   function(sys) length(surf_ids(sys)))
n_depth <- sapply(B_groups, length)

## ---- 6. Write ---------------------------------------------------------------
con <- file(out_tsv, "w")
writeLines("Supplementary Table S18. Sensitivity of core microbiome size
           to the prevalence threshold", con)
writeLines("Core = features present in >=X% of samples within each group;
           binary detection; unrarefied 97%-clustered features", con)
writeLines("", con)
writeLines(paste0("Panel A. Surface sediments (Fringe, flood, surface horizon). n: ",
                  paste(surf_sys, n_surf, sep="=", collapse=", ")), con)
close(con)
suppressWarnings(write.table(panelA, out_tsv, sep = "\t", quote = FALSE,
                             row.names = FALSE, append = TRUE))
con <- file(out_tsv, "a")
writeLines("", con)
writeLines(paste0("Panel B. Full depth profiles (flood; all horizons). n: ",
                  paste(names(n_depth), n_depth, sep="=", collapse=", ")), con)
close(con)
suppressWarnings(write.table(panelB, out_tsv, sep = "\t", quote = FALSE,
                             row.names = FALSE, append = TRUE))
message("Wrote ", out_tsv)

## ---- 7. Console output + sanity check ---------------------------------------
cat("\n== Panel A: surface (n =", paste(n_surf, collapse=","), ") ==\n"); print(panelA)
cat("\n== Panel B: depth  (n =", paste(n_depth, collapse=","), ") ==\n"); print(panelB)
cat("\n-- check @50% must match manuscript: FL 1427 | SP 1419 | Ter 70 | Cel 1230 | global 52;",
    "depth 1481 | 1695 | 123 --\n")

