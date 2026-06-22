###############################################################################
#  SIGNIFICANCIA DE LAS DIFERENCIAS EN PROCESOS DE ENSAMBLAJE ENTRE SITIOS
#  Test apropiado para datos de qpen: las comparaciones son PAREADAS y NO
#  independientes, por lo que NO se usa LMM (pseudoreplicación). Se usa test
#  exacto de Fisher / chi-cuadrado sobre la tabla de conteos grupo x proceso,
#  que es el estándar en la literatura de ensamblaje (Stegen, Ning/iCAMP).
#
#  Correr DESPUÉS de tener 'todos' (qpen_pairwise_all.csv) con las 999 perms.
###############################################################################

library(tidyverse)

# cargar resultados (ajustá la ruta si hace falta)
todos <- read_csv("06.null-models-results/qpen_pairwise_all.csv")

# Trabajamos SOLO con los tres sitios individuales para comparar entre ellos
# (el 'pooled' es un superset, no entra en la comparación entre sitios)
sitios <- c("Fossil Lagoon", "Reforma Waterfalls", "Miguelito Dike")
dat <- todos %>% filter(grupo %in% sitios)

## ---------------------------------------------------------------------------
## 1. Tabla de contingencia: sitio x proceso (conteo de pares)
## ---------------------------------------------------------------------------
tabla <- dat %>%
  count(grupo, process) %>%
  pivot_wider(names_from = process, values_from = n, values_fill = 0) %>%
  column_to_rownames("grupo") %>%
  as.matrix()

cat("Tabla de contingencia (sitio x proceso):\n")
print(tabla)

## ---------------------------------------------------------------------------
## 2. Test global: ¿difiere la composición de procesos entre los 3 sitios?
## ---------------------------------------------------------------------------
# Fisher exacto es preferible con conteos chicos en algunas celdas.
# Si Fisher tarda o falla por tabla grande, usa chi-cuadrado con simulación.
fisher_global <- tryCatch(
  fisher.test(tabla, simulate.p.value = TRUE, B = 1e5),
  error = function(e) NULL
)
chi_global <- chisq.test(tabla, simulate.p.value = TRUE, B = 1e5)

cat("\n--- Test global (los 3 sitios difieren en composición de procesos?) ---\n")
if (!is.null(fisher_global))
  cat("Fisher exacto (simulado): p =", signif(fisher_global$p.value, 4), "\n")
cat("Chi-cuadrado (simulado):  p =", signif(chi_global$p.value, 4), "\n")

## ---------------------------------------------------------------------------
## 3. Comparaciones pareadas entre sitios (post-hoc) con corrección
## ---------------------------------------------------------------------------
pares <- combn(sitios, 2, simplify = FALSE)
posthoc <- map_dfr(pares, function(pr) {
  sub <- tabla[pr, ]
  sub <- sub[, colSums(sub) > 0, drop = FALSE]   # quitar procesos ausentes
  p <- fisher.test(sub, simulate.p.value = TRUE, B = 1e5)$p.value
  tibble(sitio_1 = pr[1], sitio_2 = pr[2], p_raw = p)
})
posthoc$p_adj <- p.adjust(posthoc$p_raw, method = "BH")  # Benjamini-Hochberg

cat("\n--- Comparaciones pareadas entre sitios (Fisher, BH-corregido) ---\n")
print(posthoc)
write_csv(posthoc, "06.null-models-results/posthoc_process_differences.csv")

## ---------------------------------------------------------------------------
## 4. (Opcional) ¿El βNTI medio de cada sitio difiere de 0?
##    Esto SÍ es interpretable: |media βNTI| y si cruza ±2 indica si la
##    selección domina. Reportar media + IC bootstrap (no t-test, por la
##    no-independencia de pares). Bootstrap a nivel de muestra sería lo ideal;
##    aquí damos media y rango como descripción, no como test inferencial.
## ---------------------------------------------------------------------------
resumen_bnti <- dat %>%
  group_by(grupo) %>%
  summarise(
    n_pares      = n(),
    bNTI_medio   = round(mean(bNTI), 2),
    bNTI_sd      = round(sd(bNTI), 2),
    pct_seleccion = round(100*mean(abs(bNTI) > 2), 1),   # % pares bajo selección
    .groups = "drop"
  )
cat("\n--- Resumen descriptivo de βNTI por sitio ---\n")
print(resumen_bnti)
write_csv(resumen_bnti, "06.null-models-results/bNTI_summary_by_site.csv")

message("\nLISTO. Archivos en 06.null-models-results/")
###############################################################################