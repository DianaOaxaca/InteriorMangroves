# 07.NetworkAnalysis_Reproducible

This folder contains the reproducible network-analysis files used for the manuscript revision and reviewer/editor requests.

All paths used by the public scripts are relative to this folder.

## Real input phyloseq objects

```text
data/
  physeq_objeto_Estero.rds
  interior_mangroves/
    phyloseq.rds
  interior_impaired/
    physeq.rds
```

These files are included here again so that this analysis folder is self-contained.

## Scripts

```text
scripts/01_build_network_matrices_real_inputs.R
scripts/02_infer_ridge_networks_real_inputs.R
scripts/03_plot_network_figures_pdf_only.R
scripts/04_reviewer_controls_depth_interior_genus.R
scripts/05_depth_profile_network_trend.R
```

## What each script does

- `01_build_network_matrices_real_inputs.R`: builds ASV/species-level network matrices from the real phyloseq objects.
- `02_infer_ridge_networks_real_inputs.R`: infers ridge-regularized partial-correlation networks from the included matrices.
- `03_plot_network_figures_pdf_only.R`: regenerates PDF-only genus-level network figures from included source tables.
- `04_reviewer_controls_depth_interior_genus.R`: exports concise reviewer-control summaries.
- `05_depth_profile_network_trend.R`: recomputes the interior-only depth-profile Spearman/Kruskal trend analysis from bootstrap network metrics.

## Source data included

```text
source_data/species_asv/
source_data/species_asv/network_matrices/
source_data/species_asv/depth_profile/
source_data/interior_only/
source_data/genus_level/
source_data/genus_level/network_matrices/
```

## Figures

Only PDF versions are included.

```text
figures/species_asv/depth_profile/
figures/interior_only/
figures/genus_level/
```

## Analyses represented

1. ASV/species-level coastal-vs-interior ridge network comparison.
2. Interior-only depth-profile network trend.
3. Interior-only named-site control.
4. Genus-level coastal-vs-interior sensitivity.

## Recommended execution order

From inside `07.NetworkAnalysis_Reproducible/`:

```bash
Rscript scripts/01_build_network_matrices_real_inputs.R
Rscript scripts/02_infer_ridge_networks_real_inputs.R
Rscript scripts/05_depth_profile_network_trend.R
Rscript scripts/03_plot_network_figures_pdf_only.R
Rscript scripts/04_reviewer_controls_depth_interior_genus.R
```
