# InteriorMangroves
# Source data — figure-to-file guide

This repository archives the code and processed data for the paper
**"Long-term environmental stability might sustain reproducible depth-structured
sediment microbiomes in Mexican inland relict mangroves"**
(archived on Zenodo: https://doi.org/10.5281/zenodo.22945165).

The table below maps each **main-text figure** to the script that generates it and
to the processed data used to produce it. Software versions and the specific
parameters applied at each step are reported in the **Methods** of the paper.

> **Note on numbering.** The analysis scripts were written before the figures were
> finalized, so the `FigN` labels embedded in some script filenames follow an earlier
> ordering. The mapping below uses the **final manuscript figure numbers**; the
> relevant script is identified by its analysis series (folder/series name), which
> reflects its actual content.

| Manuscript figure | What it shows | Figure-generating script | Underlying analysis scripts | Processed data / outputs |
|---|---|---|---|---|
| **Fig. 1** | Study design and sampling map | — (map; no numerical data) | — | — |
| **Fig. 2** | Surface cross-system comparison: alpha diversity, NMDS, differential abundance, resident core | `02.5.Fig2.InteriorCoastalCompare.qmd` | `02.1.AlphaDiversity`, `02.2.BetaDiversity`, `02.3.Coremicrobiome`, `02.4.DifferentialAbundance` | `data/` (metadata, feature table), `rds/`, `Tables/` |
| **Fig. 3** | Depth-resolved comparison, inland relict vs natural Términos Lagoon | `04.5.Fig4.InteriorCoastalDepthCompare.qmd` | `04.1.AlphaDiversity`, `04.2.BetaDiversity`, `04.3.Coremicrobiome`, `04.4.DifferentialAbundance` | `data/`, `rds/`, `Tables/` |
| **Fig. 4** | Inland relict vs hydrologically impaired coastal zone | `05.5.Fig5.InteriorImpairedCoastalComparison.qmd` | `05.1.AlphaDiversity`, `05.2.BetaDiversity`, `05.3.Coremicrobiome`, `05.4.DifferentialAbundance` | `data/`, `rds/`, `Tables/` |
| **Fig. 5** | Inland relict only: vertical taxonomic and inferred functional organization | `03.7.Fig3.InteriorMangroves.qmd` | `03.1.RelativeAbundances`, `03.2.AlphaDiversity`, `03.3.BetaDiversity`, `03.4.Coremicrobiome`, `03.5.DifferentialAbundance`, `03.6.MetabolicInference` | `data/`, `rds/`, `Tables/` |
| **Fig. 6** | Community assembly (βNTI / RCbray) and microbial association networks | `06.3.Fig6.ConvergentEvidence.qmd` | `06.1.AssemblyProcesses`, `06.2.Networks` | `06.null-models-results/`, `07.NetworkAnalysis/`, `rds/` |

## Supplementary items

- **Supplementary Data 1 and 2** — interactive Krona charts (hierarchical taxonomic
  composition of all systems and of the inland depth profiles).
- **Supplementary Table S17** — core-microbiome prevalence-threshold sensitivity
  analysis (`08.coreSensitive.R`).

## Data sources

- Inland relict mangroves (this study): NCBI BioProject **PRJNA1243457**.
- Términos Lagoon (natural and impaired zones): **10.5281/zenodo.14885371**.
- Celestún Lagoon: NCBI BioProject **PRJNA550111**.

Raw sequence reads are available from the accessions above; the processed feature
tables, metadata, diversity and core tables, differential-abundance results, and
null-model and network outputs used to build the figures are contained in this
repository.
