# Ternary Partitioning of Continental-Scale Water Clarity

Reproducible R and Bookdown pipeline for **"Ternary partitioning of continental scale water clarity"** (Jawitz, Cohen, & Ross), built on the harmonized **AquaMatch** Water Quality Portal datasets (`1970–2025`) from the Environmental Data Initiative (EDI).

## Live Documentation & Technical Report

- **Live Bookdown Site:** <https://matthewross07.github.io/light_partitioning/>
- **Repository:** <https://github.com/matthewross07/light_partitioning>

The compiled Bookdown technical report is located in [`docs/index.html`](docs/index.html) and walks through:

1. **Harmonized Continental Data Sources (`R/00_download.R`)**: Programmatic retrieval of AquaMatch Secchi disk depth (`edi.1856.1`), Chlorophyll *a* (`edi.1756.2`), Dissolved Organic Carbon (`edi.1809.1`), and Total Suspended Solids / Suspended Sediment Concentration (`edi.2048.2`), alongside U.S. EPA / CEC Level I–IV ecoregion boundaries and RiverSR NHDPlusV2 flowlines.
2. **Simultaneous Matchups & Initial WQP Hygiene (`R/prep.R`)**: Co-located same-day matching, exclusion of non-Secchi turbidity-tube approximations (`sdd_approx == 0`) and MDL floor placeholders (`sdd_mdl <= 1`), gross unit-error screening, and exact same-day duplicate removal (`281,439` simultaneous 4-variable records; `490,209` simultaneous 3-LAC records).
3. **Non-Algal Particle (NAP) Derivation & Method Sensitivity (`R/prep.R`)**: Nonlinear phytoplankton dry-biomass subtraction ($\text{NAP} = \text{TSS} - 0.234\,\text{Chl }a^{0.57}$; Gohin et al. 2005) and sensitivity analyses across C:Chl conversion formulations and analytical method tiers.
4. **Two-Stage Outlier Screening & Beer–Lambert Plausibility Diagnostics (`R/prep.R`)**: Marginal `0.1–99.9%` quantile trimming (`n = 279,052`) paired with an *a priori* Beer–Lambert physical plausibility screen (`|log10(kd / kd,lit)| <= 0.80`) using independent published optical cross-sections (Bannister 1974; Kirk 1994; Morris et al. 1995; Gohin et al. 2005), removing the `1.33%` discordant WQP records while preserving 100% of true environmental extremes (`n = 275,349`).
5. **Hierarchical Ecoregion Stratification & Bio-Optical Models 1–10 (`R/ecoregions.R`, `R/models.R`)**: Block-diagonal OLS estimation across 10 spatial/typological models (`R² = 0.37 to 0.50`), Monte Carlo sample-size stability cones (`Figure S1`), and selection of Model 7 (`Level III ecoregions`, `46` groups, `R² = 0.45`, $\Delta\text{BIC} = -37{,}264$).
6. **Continental Ternary Partitioning & Manuscript Figures (`R/ternary.R`, `R/figures.R`, `R/claims.R`)**: Projection of `469,551` observations across `39,041` monitoring stations into ternary attenuation space (`Figures 1–4` and `Figure S1`), dominance-regime statistical tests, and automated 13-point pipeline verification (`R/verify.R`).

## Reproducing the Workflow

```r
source("R/00_setup.R")
fetch_all()                  # Downloads raw AquaMatch EDI datasets & spatial layers if missing
bookdown::render_book()      # Runs full pipeline, writes figs/*.png & *.pdf, renders docs/index.html
```
