# -----------------------------------------------------------------------------
# prep.R -- from harmonized AquaMatch EDI products to the training and
#           application datasets
# -----------------------------------------------------------------------------

#' Join the four AquaMatch harmonized parameter files (chla, doc, tss, sdd)
#' on (MonitoringLocationIdentifier, harmonized_utc).
#' In `tss_harmonized_final.csv` (edi.2048.2), both TSS and SSC are harmonized;
#' when both exist at the same site-timestamp, `tss` is preferred over `ssc`.
build_aquamatch_joined <- function(raw_dir = file.path(PATH$raw, 'aquamatch')) {
  download_aquamatch_data(dest_dir = raw_dir)
  chla <- data.table::fread(
    file.path(raw_dir, 'chla_harmonized_final.csv'),
    select = c('MonitoringLocationIdentifier', 'OrganizationIdentifier', 'harmonized_utc',
               'ActivityStartDate', 'ResolvedMonitoringLocationTypeName', 'lat', 'lon',
               'harmonized_value', 'tier', 'mdl_flag')
  )
  data.table::setnames(chla, c('harmonized_value', 'tier', 'mdl_flag'), c('chl_a', 'chl_tier', 'chl_mdl'))
  chla <- unique(chla, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))

  doc <- data.table::fread(
    file.path(raw_dir, 'doc_harmonized_final.csv'),
    select = c('MonitoringLocationIdentifier', 'harmonized_utc',
               'harmonized_value', 'tier', 'mdl_flag')
  )
  data.table::setnames(doc, c('harmonized_value', 'tier', 'mdl_flag'), c('doc', 'doc_tier', 'doc_mdl'))
  doc <- unique(doc, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))

  tss <- data.table::fread(
    file.path(raw_dir, 'tss_harmonized_final.csv'),
    select = c('MonitoringLocationIdentifier', 'harmonized_utc', 'parameter',
               'harmonized_value', 'tier', 'mdl_flag')
  )
  data.table::setorder(tss, MonitoringLocationIdentifier, harmonized_utc, -parameter)
  data.table::setnames(tss, c('harmonized_value', 'tier', 'mdl_flag', 'parameter'),
                       c('tss', 'tss_tier', 'tss_mdl', 'tss_param'))
  tss <- unique(tss, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))

  sdd <- data.table::fread(
    file.path(raw_dir, 'sdd_harmonized_final.csv'),
    select = c('MonitoringLocationIdentifier', 'harmonized_utc',
               'harmonized_value', 'tier', 'mdl_flag', 'approx_flag')
  )
  data.table::setnames(sdd, c('harmonized_value', 'tier', 'mdl_flag', 'approx_flag'),
                       c('secchi', 'sdd_tier', 'sdd_mdl', 'sdd_approx'))
  sdd <- unique(sdd, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))

  m <- merge(chla, doc, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))
  m <- merge(m, tss, by = c('MonitoringLocationIdentifier', 'harmonized_utc'))
  m <- merge(m, sdd, by = c('MonitoringLocationIdentifier', 'harmonized_utc'), all.x = TRUE)

  m[, type := data.table::fifelse(grepl('Lake', ResolvedMonitoringLocationTypeName), 'Lake',
                                  ResolvedMonitoringLocationTypeName)]
  data.table::setnames(m,
    c('MonitoringLocationIdentifier', 'OrganizationIdentifier', 'harmonized_utc', 'ActivityStartDate', 'lon'),
    c('SiteID', 'org', 'date_unity', 'date', 'long')
  )

  site_meta <- unique(m[, .(SiteID, type, lat, long)], by = 'SiteID')
  m[, c('type', 'lat', 'long') := NULL]
  m <- merge(m, site_meta, by = 'SiteID')

  as_tibble(m[, .(SiteID, org, date_unity, date, type, lat, long,
                  chl_a, doc, tss, secchi, tss_param,
                  chl_tier, doc_tier, tss_tier, sdd_tier,
                  chl_mdl, doc_mdl, tss_mdl, sdd_mdl, sdd_approx)])
}

#' Site inventory: one row per unique SiteID with harmonized lat/long and water-body type.
build_site_inventory <- function(am = get_aquamatch()) {
  am %>%
    filter(type != 'Facility') %>%
    distinct(SiteID, lat, long, type)
}

#' Simultaneous (same sample date-time) Chl a + DOC + TSS (+/- Secchi) records.
#'   * treat Secchi records with sdd_mdl > P$sdd_max_mdl_flag (<= 0.01 m floor
#'     placeholders in AquaMatch harmonize_sdd.R) or sdd_approx > 0 (non-Secchi
#'     tube approximations) as NA_real_
#'   * drop 'Facility' sites
#'   * hard caps: secchi < 15 m, chl_a < 1000 ug/L, tss < 1000 mg/L, doc < 50 mg/L
#'   * deduplicate exact same-day duplicate uploads at the same monitoring site
build_simultaneous <- function(am, require_secchi = TRUE) {
  out <- am %>%
    mutate(secchi = ifelse((!is.na(sdd_mdl) & sdd_mdl > P$sdd_max_mdl_flag) |
                           (isTRUE(P$exclude_sdd_approx) & !is.na(sdd_approx) & sdd_approx > 0),
                           NA_real_, secchi)) %>%
    filter(!is.na(chl_a), !is.na(doc), !is.na(tss)) %>%
    filter(type != 'Facility') %>%
    filter(chl_a < P$max_chl_ugL, tss < P$max_tss_mgL, doc < P$max_doc_mgL)
  if (require_secchi) {
    out <- out %>% filter(!is.na(secchi), secchi < P$max_secchi_m)
  }
  if (isTRUE(P$dedup_same_day)) {
    out <- out %>% distinct(SiteID, date, secchi, chl_a, doc, tss, .keep_all = TRUE)
  }
  out
}

#' Algal dry biomass from Chl a and non-algal particles (NAP, `tss_dead`).
#' biomass [mg/L] = (a/1000) * chl_a[ug/L]^b   (Gohin et al. 2005: a = 234, b = 0.57)
add_nap <- function(df, a = P$chl_biomass_a, b = if (a == 234) P$chl_biomass_b else 1) {
  df %>% mutate(chl_ratio = a, power = b,
                chl_a_biomass = exp(log(chl_ratio / 1000) + log(chl_a) * power),
                tss_dead = tss - chl_a_biomass)
}

#' Step 1 of training data: analytical reporting floors (or optional quantile trimming).
#'   1. all four variables >= 0
#'   2. sub-detection / zero values replaced via ifelse(x < floor, floor, x)
#'      using standard EPA/USGS analytical reporting floors (P$chl_floor,
#'      P$doc_floor, P$tss_floor, P$secchi_floor)
trim_simultaneous <- function(simul) {
  out <- simul %>%
    filter(if_all(c(chl_a, doc, secchi, tss), ~ !is.na(.) & . >= 0)) %>%
    mutate(
      chl_a  = ifelse(chl_a  < P$chl_floor,    P$chl_floor,    chl_a),
      doc    = ifelse(doc    < P$doc_floor,    P$doc_floor,    doc),
      tss    = ifelse(tss    < P$tss_floor,    P$tss_floor,    tss),
      secchi = ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi),
      type   = ifelse(grepl('Lake', type), 'Lake', type)
    )
  q <- P$trim_quantiles
  if (!is.null(q) && length(q) == 2 && all(is.finite(q))) {
    out <- out %>%
      filter(if_all(c(chl_a, doc, secchi, tss),
                    ~ !is.na(.) & . > quantile(., q[1]) & . < quantile(., q[2])))
  }
  out
}

#' Compute the mean and SD of the log10 Beer-Lambert residual r = log10(kd / kd_lit)
#' across the simultaneous 4-variable dataset (`trimmed`) so residuals can be
#' standardized into z-scores (|z_r| <= P$kd_discordance_sd, i.e., +/- 2.5 SD ~ 5-fold).
kd_discordance_stats <- function(trimmed = get_trimmed()) {
  d <- trimmed %>%
    add_nap() %>%
    mutate(
      secchi_r = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits),
      kd_obs   = round(1 / secchi_r, P$kd_digits),
      kd_lit   = P$lit_kw +
                 P$lit_alpha_nap * pmax(0, tss_dead) +
                 P$lit_alpha_doc * doc +
                 P$lit_alpha_chl * chl_a,
      log10_r  = log10(kd_obs / kd_lit)
    )
  mean_r <- mean(d$log10_r, na.rm = TRUE)
  sd_r   <- sd(d$log10_r, na.rm = TRUE)
  max_sd <- if (!is.null(P$kd_discordance_sd)) P$kd_discordance_sd else 2.5
  list(
    mean_r    = mean_r,
    sd_r      = sd_r,
    max_sd    = max_sd,
    tol_log10 = max_sd * sd_r,
    tol_fold  = 10^(max_sd * sd_r)
  )
}

#' A priori Beer-Lambert physical plausibility filter for WQP training records (Rule 1).
#' Standardizes the log Beer-Lambert residual r = log10(kd / kd_lit) across the
#' 4-variable dataset (mean(r) = 0.013, SD(r) = 0.270 log10 units) and retains
#' observations within +/- 2.5 SD (|z_r| <= 2.5, equivalent to a ~4.7-fold [~5x]
#' tolerance, retaining 97.51% of records and removing 2.49% discordant outliers).
clean_wqp_discordance <- function(df, max_sd = P$kd_discordance_sd) {
  if (is.null(max_sd) || !is.finite(max_sd)) return(df)
  df %>%
    mutate(
      kd_lit  = P$lit_kw +
                P$lit_alpha_nap * pmax(0, tss_dead) +
                P$lit_alpha_doc * doc +
                P$lit_alpha_chl * chl_a,
      log10_r = log10(kd / kd_lit),
      z_r     = (log10_r - mean(log10_r, na.rm = TRUE)) / sd(log10_r, na.rm = TRUE)
    ) %>%
    filter(abs(z_r) <= max_sd) %>%
    select(-kd_lit, -log10_r, -z_r)
}

#' Step 2 of training data:
#'   3. NAP from Chl a
#'   4. secchi floored at 0.01 m, rounded to 2 dp; kd = round(1/secchi, 5)
#'   5. a priori Beer-Lambert physical plausibility screen (`clean_wqp_discordance`)
build_training <- function(trimmed) {
  trimmed %>%
    add_nap() %>%
    mutate(secchi = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits),
            kd = round(1 / secchi, P$kd_digits)) %>%
    clean_wqp_discordance()
}

#' Pure-water attenuation offset used in every model: 1 / (1.5 * max Secchi).
k_w_from_training <- function(training_or_trimmed) 1 / P$k_w_divisor / max(training_or_trimmed$secchi)

#' Application data: every record with all 3 LACs (Secchi optional) whose LAC
#' values fall inside the range spanned by the training data, with NAP > 0.01 mg/L.
#' Sub-detection / zero values are floored via ifelse() to match the training set.
#' Applies Rule 2 ("Odd-One-Out" asymmetric disposition) when a measured Secchi
#' depth is present:
#'   - Upper tail (z_r > +2.5 SD): Secchi is much too shallow for the clear water
#'     chemistry (Chl a, DOC, and TSS are low and mutually consistent, e.g. a
#'     streambed bottom hit). Only `secchi` is set to NA_real_ so the valid
#'     3-constituent chemistry is retained in the ternary application dataset!
#'   - Lower tail (z_r < -2.5 SD): One of Chl a, DOC, or TSS is implausibly high
#'     relative to a clear Secchi reading (e.g. a 1000x unit error). The row is
#'     dropped because all 3 constituents are required for ternary partitioning.
build_application <- function(no_secchi, trimmed) {
  rng <- function(v) range(trimmed[[v]], na.rm = TRUE)
  r_chl <- rng('chl_a'); r_tss <- rng('tss'); r_doc <- rng('doc'); r_sec <- rng('secchi')
  stats <- kd_discordance_stats(trimmed)
  out <- no_secchi %>%
    filter(if_all(c(chl_a, doc, tss), ~ !is.na(.) & . >= 0)) %>%
    mutate(
      chl_a  = ifelse(chl_a  < P$chl_floor, P$chl_floor, chl_a),
      doc    = ifelse(doc    < P$doc_floor, P$doc_floor, doc),
      tss    = ifelse(tss    < P$tss_floor, P$tss_floor, tss),
      secchi = ifelse(!is.na(secchi) & secchi >= 0 & secchi < P$secchi_floor, P$secchi_floor, secchi)
    ) %>%
    filter(chl_a <= r_chl[2] & chl_a >= r_chl[1],
           tss   <= r_tss[2] & tss   >= r_tss[1],
           doc   <= r_doc[2] & doc   >= r_doc[1],
           (secchi <= r_sec[2] & secchi >= r_sec[1]) | is.na(secchi)) %>%
    mutate(type = ifelse(grepl('Lake', type), 'Lake', type)) %>%
    add_nap() %>%
    filter(tss_dead > P$apply_min_nap)
  if (!is.null(stats$max_sd) && is.finite(stats$max_sd)) {
    kd_lit <- P$lit_kw + P$lit_alpha_nap * pmax(0, out$tss_dead) +
              P$lit_alpha_doc * out$doc + P$lit_alpha_chl * out$chl_a
    kd_obs <- ifelse(!is.na(out$secchi) & out$secchi > 0,
                     round(1 / round(out$secchi, P$secchi_digits), P$kd_digits),
                     NA_real_)
    z_r <- (log10(kd_obs / kd_lit) - stats$mean_r) / stats$sd_r
    # Rule 2a: Upper tail (z_r > +2.5 SD) -> Secchi is the odd one out (bottom hit);
    #          drop only `secchi` (set to NA) and retain the 3-constituent chemistry!
    out$secchi <- ifelse(!is.na(z_r) & z_r > stats$max_sd, NA_real_, out$secchi)
    # Rule 2b: Lower tail (z_r < -2.5 SD) -> Chemistry constituent is the odd one out;
    #          remove the row from ternary application.
    bad_chem <- !is.na(z_r) & z_r < -stats$max_sd
    out <- out %>% filter(!bad_chem)
  }
  out
}
