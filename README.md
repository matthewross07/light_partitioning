# Ternary Partitioning of Continental-Scale Water Clarity

Open-source R and Bookdown pipeline for **"Ternary partitioning of continental scale water clarity"** (Jawitz, Cohen, & Ross), built on the harmonized **AquaMatch** Water Quality Portal datasets (`1970–2025`) from the Environmental Data Initiative (EDI).

## Live Documentation & Technical Report

- **Live Bookdown Site:** <https://matthewross07.github.io/light_partitioning/>
- **Repository:** <https://github.com/matthewross07/light_partitioning>

The compiled Bookdown technical report is located in [`docs/index.html`](docs/index.html) and covers:

1. **Harmonized Continental Data Sources (`R/00_download.R`)**: Programmatic retrieval of AquaMatch Secchi disk depth (`edi.1856.1`), Chlorophyll *a* (`edi.1756.2`), Dissolved Organic Carbon (`edi.1809.1`), and Total Suspended Solids / Suspended Sediment Concentration (`edi.2048.2`), alongside U.S. EPA / CEC Level I–IV ecoregion boundaries and RiverSR NHDPlusV2 flowlines.
2. **Simultaneous Matchups & Initial WQP Hygiene (`R/prep.R`)**: Co-located same-day matching, exclusion of non-Secchi turbidity-tube approximations (`sdd_approx == 0`) and MDL floor placeholders (`sdd_mdl <= 1`), gross unit-error screening, and same-day duplicate removal (`281,439` simultaneous 4-variable matchups; `486,689` simultaneous 3-constituent matchups).
3. **Non-Algal Particle (NAP) Derivation & Method Sensitivity (`R/prep.R`)**: Nonlinear phytoplankton dry-biomass subtraction ($\text{NAP} = \text{TSS} - 0.234\,\text{Chl }a^{0.57}$; Gohin et al. 2005) and sensitivity analyses across C:Chl formulations and analytical method tiers.
4. **Analytical Reporting Floors & Beer–Lambert Residual Screen (`R/prep.R`)**: Sub-detection analytical floor substitution (`ifelse(x < floor, floor, x)`) paired with an *a priori* Beer–Lambert residual screen (`|z_r| <= 2.5 SD`, a `~5-fold` tolerance on $r = \log_{10}(k_d / k_{d,\text{lit}})$) using independent published optical cross-sections (Bannister 1974; Kirk 1994; Morris et al. 1995; Gohin et al. 2005). Retains `97.51%` (`n = 274,445`) of training records while flagging `2.49%` (`6,994`) discordant outliers—rescuing the valid 3-constituent chemistry of `5,177` (`1.84%`) shallow streambed bottom-hit records (`Secchi` set to `NA`) in the ternary application set (`n = 479,666`).
5. **Hierarchical Ecoregion Stratification & Bio-Optical Models 1–10 (`R/ecoregions.R`, `R/models.R`)**: Block-diagonal OLS estimation across 10 spatial/typological models (`R² = 0.45 to 0.57`), Monte Carlo sample-size stability cones (`Figure S1`), and selection of Model 7 (`Level III ecoregions`, `47` groups, `R² = 0.54`, $\Delta\text{BIC} = -55{,}486$).
6. **Continental Ternary Partitioning & Manuscript Figures (`R/ternary.R`, `R/figures.R`, `R/claims.R`)**: Projection of `479,666` observations across `39,423` monitoring stations into ternary attenuation space (`Figures 1–4` and `Figure S1`), dominance-regime statistical tests, `{sonicray}` 3D GPU optical titration visuals, and automated 13-point pipeline verification (`R/verify.R`).

## Reproducing the Workflow

```r
source("R/00_setup.R")
fetch_all()                  # Downloads raw AquaMatch EDI datasets & spatial layers if missing
bookdown::render_book()      # Runs full pipeline, writes figs/*.png & *.pdf, renders docs/index.html
```
