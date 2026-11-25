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

# ---- Values below reporting limits ---------------------------------------------
# Substituting the limit (or any single value) for every result below a reporting
# limit distorts the low end of the distribution (Helsel 2006, Chemosphere 65:
# 2434-2439). Following the "fill-in" approach of Lubin et al. (2004, Environ.
# Health Perspect. 112:1691-1696), each value below the limit is replaced by a
# random draw from a lognormal distribution fitted to that variable by maximum
# likelihood for left-censored data and truncated at the limit. Fewer than 3% of
# values are affected, well under the ~30% at which Lubin et al. found the approach
# starts to bias variances. Draws are tied to the record (site, date, variable), so
# a record gets the same value in the training and the application set and results
# are reproducible without a random seed.

#' Reproducible uniform(0, 1) number for each record and variable.
key_unif <- function(SiteID, date_unity, variable) {
  key <- paste(variable, SiteID, format(date_unity, '%Y-%m-%d %H:%M:%S', tz = 'UTC'), sep = '|')
  (as.numeric(digest::digest2int(key, seed = P$random_seed)) + 2^31 + 0.5) / 2^32
}

#' Lognormal fitted by maximum likelihood to data left-censored at `limit`.
fit_censored_lognormal <- function(x, limit) {
  obs <- log(x[x >= limit]); n_below <- sum(x < limit)
  nll <- function(p) -(sum(dnorm(obs, p[1], exp(p[2]), log = TRUE)) +
                         n_below * pnorm(log(limit), p[1], exp(p[2]), log.p = TRUE))
  o <- optim(c(mean(obs), log(sd(obs))), nll, method = 'BFGS')
  list(mu = o$par[1], sigma = exp(o$par[2]), limit = limit, n_below = n_below, n = length(x))
}

#' Draws from the fitted lognormal, truncated to (0, limit), at uniform numbers `u`.
draw_below_limit <- function(u, fit) exp(qnorm(u * pnorm(log(fit$limit), fit$mu, fit$sigma), fit$mu, fit$sigma))

#' Fit the four censored distributions (Chl a, DOC, TSS and derived NAP) on every
#' record that has the three LACs.
censor_fits <- function(d) {
  d <- d %>% filter(if_all(c(chl_a, doc, tss), ~ !is.na(.) & . >= 0))
  fits <- list(chl_a = fit_censored_lognormal(d$chl_a, P$chl_floor),
               doc   = fit_censored_lognormal(d$doc,   P$doc_floor),
               tss   = fit_censored_lognormal(d$tss,   P$tss_floor))
  d2 <- add_nap(impute_lacs(d, fits), impute = FALSE)
  fits$nap <- fit_censored_lognormal(d2$nap_raw, P$nap_floor)
  fits
}

#' Replace Chl a, DOC and TSS values below their reporting limits by draws below the limit.
impute_lacs <- function(df, fits = get_censor_fits()) {
  for (v in c('chl_a', 'doc', 'tss')) {
    below <- !is.na(df[[v]]) & df[[v]] < fits[[v]]$limit
    if (any(below)) df[[v]][below] <- draw_below_limit(key_unif(df$SiteID[below], df$date_unity[below], v), fits[[v]])
  }
  df
}

#' NAP for records where the algal dry mass estimate equals or exceeds TSS (raw NAP <= 0):
#' the non-algal mass is below what TSS can resolve, so it is drawn below the TSS
#' reporting limit in the same way. Positive estimates are kept as calculated.
impute_nap <- function(nap_raw, SiteID, date_unity, fit = get_censor_fits()$nap) {
  below <- !is.na(nap_raw) & nap_raw <= 0
  if (any(below)) nap_raw[below] <- draw_below_limit(key_unif(SiteID[below], date_unity[below], 'nap'), fit)
  nap_raw
}

#' Algal dry biomass from Chl a and non-algal particles (NAP, `tss_dead`).
#' biomass [mg/L] = (a/1000) * chl_a[ug/L]^b   (Gohin et al. 2005: a = 234, b = 0.57)
#' `nap_raw` is TSS minus algal dry mass as calculated; `tss_dead` is the NAP used
#' everywhere else, in which raw values <= 0 are drawn below the TSS reporting limit
#' (impute_nap()). Set `impute = FALSE` to get the raw difference only.
add_nap <- function(df, a = P$chl_biomass_a, b = if (a == 234) P$chl_biomass_b else 1, impute = TRUE) {
  out <- df %>% mutate(chl_ratio = a, power = b,
                       chl_a_biomass = exp(log(chl_ratio / 1000) + log(chl_a) * power),
                       nap_raw = tss - chl_a_biomass,
                       tss_dead = nap_raw)
  if (impute) out$tss_dead <- impute_nap(out$nap_raw, out$SiteID, out$date_unity)
  out
}

#' Step 1 of training data: analytical reporting floors (or optional quantile trimming).
#'   1. all four variables >= 0
#'   2. Chl a, DOC and TSS values below their reporting limits (P$chl_floor,
#'      P$doc_floor, P$tss_floor) replaced by draws below the limit (impute_lacs())
trim_simultaneous <- function(simul) {
  out <- simul %>%
    filter(if_all(c(chl_a, doc, secchi, tss), ~ !is.na(.) & . >= 0)) %>%
    impute_lacs() %>%
    mutate(
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

#' Reference attenuation used by the data screen: published coefficients for diffuse
#' attenuation of PAR, exactly as printed in their sources (see P$lit_* in 00_setup.R).
#' The DOC coefficient depends on water type: Morris et al. (1995) for lakes, also
#' applied to streams, and Gallegos (2001) for estuaries.
reference_attenuation <- function(df) {
  k_doc <- ifelse(df$type == 'Estuary', P$lit_alpha_doc_estuary, P$lit_alpha_doc_inland)
  P$lit_kw + P$lit_alpha_nap * pmax(0, df$tss_dead) + k_doc * df$doc + P$lit_alpha_chl * df$chl_a
}

#' Compute the mean and SD of the log10 Beer-Lambert residual r = log10(kd / kd_lit)
#' across the simultaneous 4-variable dataset (`trimmed`) so residuals can be
#' standardized into z-scores (|z_r| <= P$kd_discordance_sd, i.e., +/- 2.5 SD, about 5.4-fold with the published reference).
kd_discordance_stats <- function(trimmed = get_trimmed()) {
  d <- trimmed %>%
    add_nap() %>%
    mutate(
      secchi_r = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits),
      kd_obs   = round(1 / secchi_r, P$kd_digits),
      kd_lit   = reference_attenuation(.),
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
      kd_lit  = reference_attenuation(.),
      log10_r = log10(kd / kd_lit),
      z_r     = (log10_r - mean(log10_r, na.rm = TRUE)) / sd(log10_r, na.rm = TRUE)
    ) %>%
    filter(abs(z_r) <= max_sd) %>%
    select(-kd_lit, -log10_r, -z_r)
}

# ---- Tests of the data screen ----------------------------------------------------

#' The 4-variable matchups with NAP, rounded Secchi and kd = 1/z_SD, before the screen.
screen_input <- function(trimmed = get_trimmed()) {
  trimmed %>% add_nap() %>%
    mutate(secchi = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits),
           kd = round(1 / secchi, P$kd_digits))
}

#' Apply the screen rule to a vector of kd values and a reference attenuation.
screen_keep <- function(kd, ref, max_sd = P$kd_discordance_sd) {
  r <- log10(kd / ref)
  abs((r - mean(r)) / sd(r)) <= max_sd
}

#' R2 of the global model (Model 1 form) on the rows kept.
global_r2 <- function(d, kd = d$kd, keep = rep(TRUE, nrow(d))) {
  X <- cbind(d$tss_dead, d$doc, d$chl_a)[keep, , drop = FALSE]; y <- kd[keep]
  k_w <- 1 / P$k_w_divisor / max(1 / y)
  cor(y, X %*% lm.fit(X, y - k_w)$coefficients)^2
}

#' Shuffled-Secchi null test. Secchi depth is shuffled among the matchups, which
#' destroys any real relation with the constituents; the screen is then applied and
#' the global model fitted. If the screen created fit by itself, R2 would be well
#' above zero. Returns one row per repeat.
screen_null_test <- function(trimmed = get_trimmed(), ref = NULL, n_rep = 20, seed = P$random_seed) {
  d <- screen_input(trimmed); if (is.null(ref)) ref <- reference_attenuation(d)
  set.seed(seed)
  map_dfr(seq_len(n_rep), function(i) {
    kd <- sample(d$kd); keep <- screen_keep(kd, ref)
    tibble(rep = i, pct_flagged = 100 * mean(!keep), r2_unscreened = global_r2(d, kd), r2_screened = global_r2(d, kd, keep))
  })
}

#' How the results depend on which published coefficients form the reference.
#' Refits all ten models for each alternative. `no screen` keeps every matchup.
screen_reference_sensitivity <- function(trimmed = get_trimmed(), eco = get_eco()) {
  d <- screen_input(trimmed)
  lin <- function(kw, nap, doc, chl) kw + nap * pmax(0, d$tss_dead) + doc * d$doc + chl * d$chl_a
  refs <- list(
    'Adopted: DOC 0.22 in lakes and streams (Morris), 0.047 in estuaries (Gallegos)' = reference_attenuation(d),
    'DOC 0.22 in all waters (Morris)'                                                = lin(P$lit_kw, P$lit_alpha_nap, P$lit_alpha_doc_inland, P$lit_alpha_chl),
    'DOC 0.047 in all waters (Gallegos)'                                             = lin(P$lit_kw, P$lit_alpha_nap, P$lit_alpha_doc_estuary, P$lit_alpha_chl),
    'Earlier, unsourced values (0.08, 0.040, 0.070, 0.016)'                          = lin(0.08, 0.040, 0.070, 0.016),
    'No screen'                                                                      = NULL)
  imap_dfr(refs, function(ref, nm) {
    keep <- if (is.null(ref)) rep(TRUE, nrow(d)) else screen_keep(d$kd, ref)
    tr <- d[keep, ]; k_w <- k_w_from_training(tr)
    m <- fit_all_models(build_zone_sets(tr, eco), k_w); p <- model_performance(m, tr)
    c1 <- m$m1$coefs %>% select(term, estimate) %>% deframe()
    a <- if (is.null(ref)) rep(NA_real_, nrow(d)) else ref * d$secchi
    r <- if (is.null(ref)) NA_real_ else log10(d$kd / ref)
    nul <- if (is.null(ref)) NA_real_ else mean(screen_null_test(trimmed, ref, n_rep = 5)$r2_screened)
    tibble(reference = nm, kept = sum(keep), flagged = sum(!keep), pct_flagged = 100 * mean(!keep),
           fold_band = if (is.null(ref)) NA_real_ else 10^(P$kd_discordance_sd * sd(r)),
           r2_m1 = p$r2[1], r2_m7 = p$r2[7], r2_m10 = p$r2[10],
           a_NAP = c1[['NAP']], a_DOC = c1[['DOC']], a_Chl_a = c1[['Chl_a']],
           kxsd_estuary = median(a[d$type == 'Estuary']), kxsd_lake = median(a[d$type == 'Lake']), kxsd_stream = median(a[d$type == 'Stream']),
           r2_shuffled = nul)
  })
}

#' Step 2 of training data:
#'   3. NAP from Chl a
#'   4. secchi floored at 0.01 m, rounded to 2 dp; kd = round(1/secchi, 5)
#'   5. screen against the published reference attenuation (`clean_wqp_discordance`)
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
#' values fall inside the range spanned by the training data. No record is dropped
#' for its NAP value: where the algal dry mass estimate equals or exceeds TSS, NAP is
#' drawn below the TSS reporting limit (add_nap() / impute_nap()), in training and
#' application alike.
#' Values below reporting limits are drawn below the limit exactly as in the training set.
#' Applies Rule 2 ("Odd-One-Out" asymmetric disposition) when a measured Secchi
#' depth is present:
#'   - Upper tail (z_r > +2.5 SD): Secchi is much too shallow for the clear water
#'     chemistry (Chl a, DOC, and TSS are low and mutually consistent, e.g. a
#'     streambed bottom hit). Only `secchi` is set to NA_real_ so the valid
#'     3-constituent chemistry is retained in the ternary application dataset!
#'   - Lower tail (z_r < -2.5 SD): One of Chl a, DOC, or TSS is implausibly high
#'     relative to a clear Secchi reading (e.g. a 1000x unit error). The row is
#'     dropped because all 3 constituents are required for ternary partitioning.
#' `training` supplies the concentration ranges. `stats` must be the screen statistics
#' from the full 4-variable set (kd_discordance_stats(get_trimmed())), i.e. the same
#' mean and SD used to build the training set, so a record is treated the same way in
#' both. (Taking them from the already-screened training set narrows the band from
#' ~4.7-fold to ~3.8-fold and screens the application set a second time.)
build_application <- function(no_secchi, training, stats = kd_discordance_stats(get_trimmed())) {
  rng <- function(v) range(training[[v]], na.rm = TRUE)
  r_chl <- rng('chl_a'); r_tss <- rng('tss'); r_doc <- rng('doc'); r_sec <- rng('secchi')
  out <- no_secchi %>%
    filter(if_all(c(chl_a, doc, tss), ~ !is.na(.) & . >= 0)) %>%
    impute_lacs() %>%
    mutate(
      secchi = ifelse(!is.na(secchi) & secchi >= 0 & secchi < P$secchi_floor, P$secchi_floor, secchi)
    ) %>%
    # upper bounds only: every value below a reporting limit is a draw below that limit,
    # so the lower end of the training range is the limit itself
    filter(chl_a <= r_chl[2], tss <= r_tss[2], doc <= r_doc[2],
           (secchi <= r_sec[2] & secchi >= r_sec[1]) | is.na(secchi)) %>%
    mutate(type = ifelse(grepl('Lake', type), 'Lake', type)) %>%
    add_nap()
  if (!is.null(stats$max_sd) && is.finite(stats$max_sd)) {
    kd_lit <- reference_attenuation(out)
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
