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
  stopifnot('package "digest" is needed for reproducible draws below reporting limits' = requireNamespace('digest', quietly = TRUE))
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

  # -- AquaMatch Secchi MDL floor & WQP QA/QC handling (De La Torre et al. 2025, edi.1856.1):
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

  # -- Data screen (Rules 1 & 2): reference attenuation from published coefficients.
  #    Each matchup's 1/z_SD is compared with a reference attenuation
  #        K_ref = k_w + k_NAP * NAP + k_DOC * DOC + k_Chl * Chl a
  #    whose four coefficients are published values for diffuse attenuation of PAR,
  #    used exactly as printed (reference_attenuation() in R/prep.R):
  #      k_w   = 0.1    m^-1       Gohin et al. (2005) Remote Sens. Environ. 95:29-46, Eq. 21
  #      k_NAP = 0.0625 m^2 g^-1   Gohin et al. (2005), Eq. 21 (non-living suspended matter,
  #                                 defined by the same subtraction we use, their Eqs. 15-16)
  #      k_DOC = 0.22   m^2 g^-1   lakes and streams: Morris et al. (1995) Limnol. Oceanogr.
  #                                 40:1381-1391, Table 3 (p. 1387), K_dPAR in 63 lakes
  #            = 0.047  m^2 g^-1   estuaries: Gallegos (2001) Estuaries 24:381-397, Table 5
  #                                 (p. 393), mean of four Chesapeake Bay stations
  #      k_Chl = 0.016  m^2 mg^-1  Bannister (1974) Limnol. Oceanogr. 19:1-12, p. 6 (mean of
  #                                 four estimates, 0.013-0.020, p. 10)
  #    Rule 1: r = log10[(1/z_SD) / K_ref]; z_r = (r - mean(r)) / SD(r); keep |z_r| <= 2.5.
  #    Published coefficients are per unit of K_d and 1/z_SD = K_d / a (Eq. 1), so a
  #    constant a shifts every r by the same amount and cancels in z_r: no conversion
  #    of the coefficients is needed, and none is applied.
  #    Rule 2 (application set): for z_r > +2.5 (Secchi much shallower than the
  #    constituents imply) only Secchi is set aside and the three constituents are
  #    kept; for z_r < -2.5 (Secchi much deeper than they imply) the record is dropped.
  lit_kw                   = 0.1,
  lit_alpha_nap            = 0.0625,
  lit_alpha_doc_inland     = 0.22,    # lakes and streams
  lit_alpha_doc_estuary    = 0.047,   # estuaries
  lit_alpha_chl            = 0.016,
  kd_discordance_sd        = 2.5,

  # -- Chl a -> algal dry biomass (Gohin et al. 2005, Remote Sens. Environ. 95:29-46, Eqs. 15-16)
  #    biomass [mg/L] = (234/1000) * chl_a[ug/L]^0.57
  #    NAP ("tss_dead") = TSS - biomass
  chl_biomass_a  = 234,
  chl_biomass_b  = 0.57,
  chl_ratio_sensitivity = c(50, 100, 200, 234),  # explored; 234*Chl^0.57 chosen

  # -- Analytical reporting floors (replacing marginal quantile trimming):
  #    Rather than dropping the 0.1% and 99.9% tails of each variable (which
  #    discards valid oligotrophic clear-water and storm/bloom extremes), values
  #    reported as zero or below standard analytical method detection limits (MDLs)
  #    are replaced by random draws below the limit (impute_lacs() in R/prep.R;
  #    Lubin et al. 2004; Helsel 2006). Secchi has no values below its limit:
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
  dominance_threshold = 0.70,   # a LAC is "dominant" where it supplies > 70% of non-water attenuation
  majority_threshold  = 0.50,   # ... the "majority contributor" where it supplies > 50%
  nap_floor           = 0.50    # mg/L: where algal dry mass >= TSS, NAP is drawn below this limit (the TSS reporting limit)
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
