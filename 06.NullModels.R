###############################################################################
#  ANÁLISIS DE PROCESOS DE ENSAMBLAJE — versión iCAMP (RÁPIDA)
#  βNTI (selección) + RCbray (procesos estocásticos)  vía qpen()
#
#  iCAMP::qpen implementa Stegen et al. 2013 con backend optimizado en C,
#  maneja miles de ASVs sin filtrar y sin colgarse. SIN recorte de taxones.
#
#  Para: manuscrito "Long-term environmental stability..." (Communications Biology)
#  Responde a R1 (assembly processes), R3 (selection vs dispersal vs drift), R2.
###############################################################################

## ---------------------------------------------------------------------------
## 0. PAQUETES
## ---------------------------------------------------------------------------
# install.packages("iCAMP")    # CRAN
# install.packages(c("phyloseq","qiime2R","tidyverse"))  # qiime2R desde GitHub
library(iCAMP)
library(phyloseq)
library(qiime2R)
library(tidyverse)

set.seed(123)
nperm  <- 999      # randomizaciones (iCAMP las maneja eficientemente)
ncores <- 4       # qpen paraleliza internamente; 14 está bien en tu equipo

## ---------------------------------------------------------------------------
## 1. IMPORTAR DATOS  (idéntico a tu pipeline InteriorMangroves)
## ---------------------------------------------------------------------------
ps <- qza_to_phyloseq(
  features = "data/ASV_table_filter_freq218_emc.qza",
  tree     = "data/rooted-tree-iqtree.qza",
  taxonomy = "data/taxonomy.qza",
  metadata = "data/metadata2.tsv"
)

# Renombrar sistemas como en 03.1.RelativeAbundances
metadata <- as(sample_data(ps), "data.frame") %>%
  mutate(Location = case_when(
    Location == "El Cacahuate"    ~ "Fossil Lagoon",
    Location == "La Piedad"       ~ "Reforma Waterfalls",
    Location == "Dique Miguelito" ~ "Miguelito Dike",
    TRUE ~ Location))
sample_data(ps) <- sample_data(metadata)

## ---------------------------------------------------------------------------
## 2. PREPARAR MATRICES PARA qpen()
## ---------------------------------------------------------------------------
# qpen necesita:
#   comm : matriz comunidad, MUESTRAS en filas, ASVs en columnas (conteos)
#   pd   : matriz de distancias filogenéticas entre ASVs (cophenetic del árbol)
# SIN filtrar ASVs.

comm <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) comm <- t(comm)        # -> muestras en filas
tree <- phy_tree(ps)

# Emparejar árbol y tabla (mismo set de tips)
comunes <- intersect(colnames(comm), tree$tip.label)
comm <- comm[, comunes, drop = FALSE]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, comunes))

cat("Muestras:", nrow(comm), " ASVs:", ncol(comm), "\n")

# Matriz de distancia filogenética (una sola vez; qpen la reutiliza)
# cophenetic.phylo puede ser pesada en RAM pero se calcula UNA vez, no por permutación
pd <- cophenetic(tree)
pd <- pd[colnames(comm), colnames(comm)]      # mismo orden que comm

## ---------------------------------------------------------------------------
## 3. CORRER qpen POR GRUPO  (cada sitio interior + pooled)
## ---------------------------------------------------------------------------
meta_df <- data.frame(sample_data(ps))
interiores <- c("Fossil Lagoon", "Reforma Waterfalls", "Miguelito Dike")

run_qpen <- function(samples, etiqueta) {
  message(">>> qpen: ", etiqueta, " (", length(samples), " muestras)")
  cm <- comm[rownames(comm) %in% samples, , drop = FALSE]
  cm <- cm[, colSums(cm) > 0, drop = FALSE]   # quita ASVs ausentes en el subset
  pdi <- pd[colnames(cm), colnames(cm)]
  
  res <- qpen(
    comm = cm,
    pd   = pdi,
    rand = nperm,
    nworker = ncores,
    sig.bNTI = 2,        # umbral |βNTI| = 2
    sig.rc   = 0.95      # umbral |RCbray| = 0.95
  )
  saveRDS(res, paste0("qpen_", gsub("[^A-Za-z0-9]","_",etiqueta), ".rds"))
  res$result$grupo <- etiqueta
  res$result
}

resultados <- list()

# por sitio
for (sis in interiores) {
  ids <- rownames(meta_df)[meta_df$Location == sis]
  resultados[[sis]] <- run_qpen(ids, sis)
}

# pooled (todos los interiores juntos) — clave para tu hipótesis espacial
ids_all <- rownames(meta_df)[meta_df$Location %in% interiores]
resultados[["pooled"]] <- run_qpen(ids_all, "All interiors (pooled)")

todos <- bind_rows(resultados)
write_csv(todos, "06.null-models-results/qpen_pairwise_all.csv")

## ---------------------------------------------------------------------------
## 4. RESUMEN: % de cada proceso por grupo
## ---------------------------------------------------------------------------
# qpen ya clasifica cada par en la columna 'process':
#   Heterogeneous.Selection / Homogeneous.Selection /
#   Dispersal.Limitation / Homogenizing.Dispersal / Drift.and.Others
resumen <- todos %>%
  count(grupo, process) %>%
  group_by(grupo) %>%
  mutate(porcentaje = round(100 * n / sum(n), 1)) %>%
  ungroup()
write_csv(resumen, "06.null-models-results/qpen_process_summary.csv")
print(resumen)

## ---------------------------------------------------------------------------
## 5. FIGURAS
## ---------------------------------------------------------------------------
# distribución de βNTI por grupo
p1 <- ggplot(todos, aes(x = grupo, y = bNTI, fill = grupo)) +
  geom_hline(yintercept = c(-2,2), linetype = "dashed", color = "grey40") +
  geom_boxplot(width = 0.6, outlier.alpha = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.6) +
  labs(x = NULL, y = expression(beta*"NTI"),
       title = "Community assembly within inland relict mangroves",
       caption = "βNTI < -2 = homogeneous selection (expected under long-term stability)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
ggsave("06.null-models-results/Fig_bNTI_qpen.pdf", p1, width = 7.5, height = 5)

# stacked de procesos
p2 <- ggplot(resumen, aes(x = grupo, y = porcentaje, fill = process)) +
  geom_col() +
  labs(x = NULL, y = "% of pairwise comparisons", fill = "Assembly process") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave("06.null-models-results/Fig_assembly_stacked_qpen.pdf", p2, width = 8, height = 5)

message("LISTO. Revisa qpen_process_summary.csv y las figuras.")
###############################################################################