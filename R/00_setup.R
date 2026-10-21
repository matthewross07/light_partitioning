# -----------------------------------------------------------------------------
# 00_setup.R  -- packages, paths, and EVERY analysis decision as a named value
#
# Every chapter of the book sources this file. If you want to know "what number
# did we use for X?", the answer is in `P` below, with a comment documenting
# its scientific rationale and literature provenance.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(feather)
  library(patchwork)
  library(scales)
  library(knitr)
  library(kableExtra)
})

sf_use_s2(FALSE)   # all spatial ops are on projected (Albers) ecoregion layers
options(dplyr.summarise.inform = FALSE, readr.show_col_types = FALSE)

# ---- Paths -------------------------------------------------------------------
PATH <- list(
  raw      = 'data/raw',
  derived  = 'data/derived',
  figs     = 'figs'
)
for (d in PATH[c('raw', 'derived', 'figs')]) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Analysis parameters (decisions) ----------------------------------------
P <- list(
  # -- Initial screening of AquaMatch harmonized data:
  #    Hard caps applied before quantile trimming to remove gross unit errors:
  max_secchi_m   = 15,      # secchi < 15 m
  max_chl_ugL    = 1000,    # chl_a < 1000 ug/L
  max_tss_mgL    = 1000,    # tss < 1000 mg/L
  max_doc_mgL    = 50,      # doc < 50 mg/L

  # -- AquaMatch Secchi MDL floor & WQP QA/QC handling (Brousil et al. 2025, edi.1856.1):
  #    In AquaMatch harmonize_sdd.R, WQP Secchi entries <= 0.01 m (including
  #    0.0 m placeholders) were mutated to 0.01 m with mdl_flag >= 2. Treating
  #    sdd_mdl >= 2 as missing (NA) prevents the artificial 0.01 m floor from
  #    corrupting the 0.1% quantile trim while keeping valid 3-LAC rows in the
  #    application set. Non-Secchi turbidity-tube approximations (sdd_approx > 0,
  #    Tier 3) are also excluded from Secchi training records, and exact same-day
  #    duplicate WQP uploads at the same station are deduplicated.
  sdd_max_mdl_flag = 1L,
  exclude_sdd_approx = TRUE,
  dedup_same_day     = TRUE,

  # -- A priori bio-optical physical plausibility screen (Rules 1 & 2):
  #    Rule 1 (Standard Deviation Screen): Across simultaneous 4-variable matchups,
  #    compute the log Beer-Lambert residual r = log10(kd / kd_lit) using
  #    independent published bio-optical constants (zero fitted parameters):
  #      k_w,lit   = 0.080 m^-1      (Kirk 1994; Lee et al. 2015)
  #      k_NAP,lit = 0.040 m^2 g^-1  (Kirk 1994; Gohin et al. 2005)
  #      k_DOC,lit = 0.070 m^2 g^-1  (Morris et al. 1995; Griffin et al. 2018)
  #      k_Chl,lit = 0.016 m^2 mg^-1 (Bannister 1974; Kirk 1994)
  #    Across the 281,439 matchups, mean(r) = 0.013 and SD(r) = 0.270 log10 units.
  #    Standardizing z_r = (r - mean(r)) / SD(r), we require |z_r| <= 2.5 SD
  #    (equivalent to |r - mean(r)| <= 0.675 log10 units, a ~4.7-fold [~5x] tolerance),
  #    retaining 97.51% (274,445) of training observations and flagging 2.49% (6,994).
  #    Rule 2 ("Odd-One-Out" Disposition in Application):
  #      - Upper tail (z_r > +2.5 SD, 1.84%): Secchi is the odd one out (shallow
  #        streambed bottom hit in clear water where Chl a, DOC, and TSS agree).
  #        Excluded from regression training, but in ternary application only
  #        Secchi is set to NA so the valid 3-constituent chemistry is retained.
  #      - Lower tail (z_r < -2.5 SD, 0.65%): One of the chemistry constituents
  #        (Chl a, DOC, or TSS) is the odd one out (e.g., 1000x unit error
  #        contradicting a clear Secchi reading). Excluded from both training
  #        and ternary application.
  lit_kw                   = 0.080,
  lit_alpha_nap            = 0.040,
  lit_alpha_doc            = 0.070,
  lit_alpha_chl            = 0.016,
  kd_discordance_sd        = 2.5,
  max_log10_kd_discordance = 0.675, # 2.5 * SD(r) = 2.5 * 0.2700 = 0.675 log10 units (~4.7-fold)

  # -- Chl a -> algal dry biomass (Gohin et al. 2005, Int. J. Remote Sens. 26:1639-1663)
  #    biomass [mg/L] = (234/1000) * chl_a[ug/L]^0.57
  #    NAP ("tss_dead") = TSS - biomass
  chl_biomass_a  = 234,
  chl_biomass_b  = 0.57,
  chl_ratio_sensitivity = c(50, 100, 200, 234),  # explored; 234*Chl^0.57 chosen

  # -- Analytical reporting floors (replacing marginal quantile trimming):
  #    Rather than dropping the 0.1% and 99.9% tails of each variable (which
  #    discards valid oligotrophic clear-water and storm/bloom extremes), values
  #    reported as zero or below standard analytical method detection limits (MDLs)
  #    are floored via ifelse(x < floor, floor, x):
  #      chl_floor    = 0.10 ug/L (standard fluorometric/spectrophotometric Chl a MDL)
  #      doc_floor    = 0.25 mg/L (EPA Method 415.1 / USGS NWQL DOC reporting limit)
  #      tss_floor    = 0.50 mg/L (standard gravimetric TSS reporting limit)
  #      secchi_floor = 0.01 m    (minimum Secchi reporting precision)
  trim_quantiles = NULL,
  chl_floor      = 0.10,
  doc_floor      = 0.25,
  tss_floor      = 0.50,
  secchi_floor   = 0.01,
  secchi_digits  = 2,       # secchi rounded to 2 decimals before kd = 1/secchi
  kd_digits      = 5,       # kd rounded to 5 decimals (Lee et al. 2015, 2018: Kd = 1/z_SD)

  # -- Pure-water attenuation used as a fixed offset in every model:
  #    k_w = 1 / (1.5 * max(training secchi)) = 0.0478 m^-1 (with max secchi = 13.96 m).
  #    Used CONSISTENTLY both when fitting models and when back-calculating
  #    predicted Secchi depth.
  k_w_divisor    = 1.5,

  # -- Minimum group size = 100 (matching manuscript Section 2.1 & Table 1):
  #    Both ecoregion groups (zone_min) and ecoregion x water-body type groups
  #    (type_min) require > 100 observations to receive their own coefficients.
  zone_min       = 100,
  type_min       = 100,

  # -- Application to all LAC observations:
  apply_min_nap  = 0.01,    # application rows need NAP > 0.01 mg/L
  apply_min_n    = 100,     # ecoregion needs > 100 application rows and positive slopes to use Model 7
  require_positive_coefs = TRUE, # fall back to global Model 1 if any group slope <= 0 (all have p > 0.01)

  # -- Random-subset coefficient stability (Fig. S1):
  random_divisors = c(50, 100, 250, 500, 1000, 2000, 5000, 10000),  # + 1 group (all data)
  random_reps     = 10,
  random_seed     = 20240826,

  # -- Example-site selection for Figs 3-4:
  long_term_min_years = 8,      # Delaware: DRBC sites with >= 8 distinct years
  fig4_min_obs        = 7,      # CDF: sites with > 7 observations
  fig4_s7_min_obs     = 100,    # insets 2-3 drawn from sites with > 100 observations
  dominance_threshold = 0.80    # ">80% contribution" classes for Secchi depth comparison
)

# ---- Plot palettes ----------------------------------------------------------
TYPE_COLS <- c(Estuary = 'grey0', Lake = 'blue', Stream = 'red3')
LAC_COLS  <- c(chl_kd = 'green4', doc_kd = 'black', tss_kd = 'tan')

# ---- Small helpers ----------------------------------------------------------
#' Compute-once cache: evaluate `expr` and save to `data/derived/<name>.rds`,
#' or read the saved file if it already exists. Delete the file to recompute.
cache_rds <- function(name, expr, refresh = FALSE) {
  f <- file.path(PATH$derived, paste0(name, '.rds'))
  if (file.exists(f) && !refresh) return(readRDS(f))
  out <- force(expr)
  saveRDS(out, f)
  out
}

fmt_n <- function(x) format(x, big.mark = ',', scientific = FALSE, trim = TRUE)

# Source the function library in pipeline order (starting with 00_download.R)
for (f in c('R/00_download.R', 'R/prep.R', 'R/ecoregions.R', 'R/models.R',
            'R/ternary.R', 'R/figures.R', 'R/claims.R', 'R/pipeline.R', 'R/verify.R')) {
  if (file.exists(f)) source(f)
}
