# ATN021B
Analysis of the ATN021B cohort


Supports Metabolomics Investigation of Antiretroviral Therapy and HIV-associated Bone Mineral Density Loss in Male Adolescents with HIV

## Structure

- `data/` — input data
- `R/preprocessing/` — data cleaning, run first
- `R/analysis/` — scripts for manuscript tables
- `R/figures/` — scripts for manuscript figures

## Scripts and corresponding outputs

### Tables

| Output | Script |
|---|---|
| Table 1a | `R/analysis/tableone.R` |
| Table 1b | `R/analysis/tableone_collapsedTDF.R` |
| Supp Table 2, 3 | `R/analysis/bone_linearregression.R` |
| Supp Table 4, 10 | `R/analysis/emmeans.R` |
| Supp Table 5 | `R/analysis/boruta_Group1and2_noARVs.R` |
| Supp Table 6 | `R/analysis/emmeansborutaoverlap_12.R` |
| Supp Table 7 | `R/analysis/emmeans_2vs4ab_noARVs.R` |
| Supp Table 8 | `R/analysis/Boruta_2vs4ab_noadj.R` |
| Supp Table 9 | `R/analysis/emmeansborutaoverlap_24ab.R` |
| Supp Table 11 | `R/analysis/Unstratified_Bone.vs.Metabolites_Correlation.R` |

### Figures

| Output | Script |
|---|---|
| Figure 1a | Made manually, not in repo |
| Figure 1b | `R/figures/bone_summary.R` |
| Figure 1c | `R/figures/bone_summary_collapsedTDF.R` |
| Figure 2 | `R/figures/complexheatmap_emmeans12.R` |
| Figure 3a | `R/figures/pca_plots.R` |
| Figure 3b | `R/figures/volcano_plots.R` |
| Figure 3c | `R/figures/bubble_plots.R` |
| Figure 4a | `R/figures/metabolite_plotting.R` |
| Figure 4b | `R/figures/collapsedTDF_mixOmics.R` |
| Supp Figure 1 | `R/figures/contrasts_heatmap.R` |
| Supp Figure 2a | `R/analysis/emmeansborutaoverlap_12.R` |
| Supp Figure 2b | `R/analysis/emmeansborutaoverlap_24ab.R` |
| Supp Figure 2c | `R/figures/emmeans_spearman_intersect_2vs4ab.R` |

## Preprocessing

`R/preprocessing/remove_ARVs_code.R` removes antiretroviral (ARV)-related features prior to all analysis

## Requirements

R (≥ 4.0). Key packages: `mixOmics`, `Boruta`, `ComplexHeatmap`, `emmeans`, `tableone`.

## Contact

Iris Gu — irisgu@g.ucla.edu
