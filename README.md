# InteriorMangroves
# Source data — figure-to-file guide

This repository archives the code and processed data for the paper
**"Long-term environmental stability might sustain reproducible depth-structured
sediment microbiomes in Mexican inland relict mangroves"**
(archived on Zenodo: https://doi.org/10.5281/zenodo.22945165).

## How the figures are built

Each main-text figure is assembled from several panels. **Every panel originates in a
subsection script** of the corresponding analysis section (e.g. `AlphaDiversity`,
`BetaDiversity`, `Coremicrobiome`, `DifferentialAbundance`, `MetabolicInference`),
which computes the underlying values and saves the panel as an `.rds` object in `rds/`
(numerical results are also written to `Tables/` and, for assembly and networks, to
`06.null-models-results/` and `07.NetworkAnalysis/`). The `FigN.*.qmd` script then
**reads those `.rds` objects and combines them into a single multi-panel figure with
cowplot**; final cosmetic editing is done in Inkscape. The numerical source data for
each panel therefore live in its subsection outputs (`rds/` and `Tables/`), not in the
combining script.

> **Note on numbering.** Some combining-script filenames carry `FigN` labels from an
> earlier figure ordering. The mapping below uses the **final manuscript figure
> numbers**; scripts are identified by their analysis section, which reflects content.

## Figure 1 — Study design and sampling map
No numerical source data (cartographic figure).

## Figure 2 — Surface cross-system comparison (inland relict vs coastal)
Section `02.SurfaceComparison`; combined in `02.5.Fig2.InteriorCoastalCompare.qmd`.

| Panel | Content | Originating subsection |
|---|---|---|
| 2a | Alpha diversity (q1, q2; q0 in Suppl.) | `02.1.AlphaDiversity.qmd` |
| 2b | NMDS (Bray–Curtis) + PERMANOVA/PERMDISP | `02.2.BetaDiversity.qmd` |
| 2c | Differential abundance (ANCOM-BC2, LFC) | `02.4.DifferentialAbundance.qmd` |
| 2d | Resident core UpSet + relative-abundance heatmap | `02.3.Coremicrobiome.qmd` |

## Figure 3 — Depth-resolved comparison (inland relict vs natural Términos)
Section `04.DepthComparison`; combined in `04.5.Fig4.InteriorCoastalDepthCompare.qmd`.

| Panel | Content | Originating subsection |
|---|---|---|
| 3a, 3b | Alpha diversity across depth | `04.1.AlphaDiversity.qmd` |
| 3c | NMDS (Bray–Curtis) + statistics | `04.2.BetaDiversity.qmd` |
| 3d | Depth-profile core UpSet + heatmap | `04.3.Coremicrobiome.qmd` |
| 3e | Differential abundance (ANCOM-BC2, LFC) | `04.4.DifferentialAbundance.qmd` |

## Figure 4 — Inland relict vs hydrologically impaired coastal zone
Section `05.InteriorImpairedCoastalComparison`; combined in `05.5.Fig5.InteriorImpairedCoastalComparison.qmd`.

| Panel | Content | Originating subsection |
|---|---|---|
| 4a, 4b | Alpha diversity (observed richness) | `05.1.AlphaDiversity.qmd` |
| 4c | NMDS (Bray–Curtis) + statistics | `05.2.BetaDiversity.qmd` |
| 4d | Resident core UpSet + heatmap | `05.3.Coremicrobiome.qmd` |

## Figure 5 — Inland relict only: vertical taxonomic and inferred functional organization
Section `03.InteriorMangroves`; combined in `03.7.Fig3.InteriorMangroves.qmd`.

| Panel | Content | Originating subsection |
|---|---|---|
| 5a | Alpha diversity (q0, Faith's PD) | `03.2.AlphaDiversity.qmd` |
| 5b | Differential abundance across depth | `03.5.DifferentialAbundance.qmd` |
| 5c, 5d | Weighted / unweighted UniFrac NMDS | `03.3.BetaDiversity.qmd` |
| 5e | db-RDA (constrained by environmental variables) | `03.3.BetaDiversity.qmd` |
| 5f | Inferred functional pathways heatmap | `03.6.MetabolicInference.qmd` |
| 5g | Inland resident core UpSet | `03.4.Coremicrobiome.qmd` |

(Relative-abundance profiles used in the Supplementary figures: `03.1.RelativeAbundances.qmd`.)

## Figure 6 — Community assembly and microbial association networks
Section `06.ConvergentEvidence`; combined in `06.3.Fig6.ConvergentEvidence.qmd`.

| Panel | Content | Originating subsection |
|---|---|---|
| 6a | βNTI distributions | `06.1.AssemblyProcesses.qmd` |
| 6b | Assembly-process proportions (βNTI/RCbray) | `06.1.AssemblyProcesses.qmd` |
| 6c | Network metrics across depth | `06.2.Networks.qmd` |

Assembly and network numerical outputs are in `06.null-models-results/` and
`07.NetworkAnalysis/`.

## Supplementary items
- **Supplementary Data 1 and 2** — interactive Krona charts (hierarchical taxonomic
  composition of all systems and of the inland depth profiles).
- **Supplementary Table S17** — core-microbiome prevalence-threshold sensitivity
  analysis (`08.coreSensitive.R`).

## Data sources
- Inland relict mangroves (this study): NCBI BioProject **PRJNA1243457**.
- Términos Lagoon (natural and impaired zones): **10.5281/zenodo.14885371**.
- Celestún Lagoon: NCBI BioProject **PRJNA550111**.

Software versions and the parameters applied at each step are reported in the Methods.

---
Raw sequence reads are available from the accessions above; the processed feature tables, metadata, diversity and core tables, differential-abundance results, and null-model and network outputs used to build the figures are contained in this repository.
---
