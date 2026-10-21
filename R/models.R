# -----------------------------------------------------------------------------
# models.R -- the ten light-attenuation models, their summaries, coefficient
#             stability, and Model 7 application
#
# Model form (paper Eq. 2-3), fitted by OLS with no intercept and a fixed
# pure-water offset k_w:
#     kd - k_w = a_NAP * NAP + a_DOC * DOC + a_Chl * Chl_a,   kd = 1 / Secchi
# Coefficients vary by ecoregion (`:code`) and/or water type (`:type`).
# Because interaction-only terms with no shared main effects form a block-diagonal
# design matrix, per-group OLS via `lm.fit()` is mathematically identical to
# joint `glm()` (to ~10^-15) and ~500x faster for Level 3 and Level 4 models.
# -----------------------------------------------------------------------------

MODEL_KEY <- tribble(
  ~model, ~level,   ~by_type, ~data_key, ~scope,
  1,      'USA',    FALSE,    'global',  'Continental (all waters)',
  2,      'USA',    TRUE,     'global',  'Continental × water-body type',
  3,      'Level 1', FALSE,   'L1',      'EPA Level I ecoregion',
  4,      'Level 1', TRUE,    'L1',      'EPA Level I × water-body type',
  5,      'Level 2', FALSE,   'L2',      'EPA Level II ecoregion',
  6,      'Level 2', TRUE,    'L2',      'EPA Level II × water-body type',
  7,      'Level 3', FALSE,   'L3',      'EPA Level III ecoregion',
  8,      'Level 3', TRUE,    'L3',      'EPA Level III × water-body type',
  9,      'Level 4', FALSE,   'L4',      'EPA Level IV ecoregion',
  10,     'Level 4', TRUE,    'L4',      'EPA Level IV × water-body type'
)

#' Fit one model using exact block-diagonal OLS per group.
#' `type_min`: ecoregion x type groups need > type_min rows.
fit_one_model <- function(dat, by_type, k_w, type_min = P$type_min) {
  dat <- ungroup(dat)
  regional <- 'code' %in% names(dat)
  if (regional && by_type) {
    dat <- dat %>% add_count(code, type, name = 'n_grp') %>% filter(n_grp > type_min)
  }
  grp_df <- dat %>%
    mutate(.code = if (regional) as.character(code) else 'USA',
           .type = if (by_type) as.character(type) else 'All')

  fits <- grp_df %>%
    group_by(.code, .type) %>%
    group_modify(function(d, k) {
      X <- cbind(tss_dead = d$tss_dead, doc = d$doc, chl_a = d$chl_a)
      y <- d$kd - k_w
      fit <- lm.fit(X, y)
      rss_g <- sum(fit$residuals^2)
      df_g <- nrow(d) - fit$rank
      s2_g <- if (df_g > 0) rss_g / df_g else NA_real_
      se_g <- if (!is.na(s2_g) && fit$rank == 3) {
        sqrt(diag(chol2inv(fit$qr$qr[1:3, 1:3, drop = FALSE])) * s2_g)
      } else {
        rep(NA_real_, 3)
      }
      t_g <- fit$coefficients / se_g
      p_g <- if (df_g > 0) 2 * pt(-abs(t_g), df = df_g) else rep(NA_real_, 3)
      tibble(
        term      = c('NAP', 'DOC', 'Chl_a'),
        estimate  = unname(fit$coefficients),
        std_error = unname(se_g),
        t_value   = unname(t_g),
        p_value   = unname(p_g),
        n_group   = nrow(d)
      )
    }) %>%
    ungroup() %>%
    rename(code = .code, type = .type)

  cf_wide <- fits %>%
    select(code, type, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate)

  dat_fitted <- grp_df %>%
    left_join(cf_wide, by = c('.code' = 'code', '.type' = 'type')) %>%
    mutate(.fitted = tss_dead * NAP + doc * DOC + chl_a * Chl_a)

  n <- nrow(dat)
  p <- sum(!is.na(fits$estimate))
  rss <- sum((dat$kd - k_w - dat_fitted$.fitted)^2)
  ll <- -n / 2 * (log(2 * pi) + log(rss / n) + 1)

  list(
    coefs         = fits,
    data          = dat,
    fitted.values = dat_fitted$.fitted,
    k_w           = k_w,
    rss           = rss,
    nobs          = n,
    n_coef        = p,
    aic           = -2 * ll + 2 * (p + 1),
    bic           = -2 * ll + log(n) * (p + 1),
    regional      = regional,
    by_type       = by_type
  )
}

#' Fit all ten models. Returns a named list `m1` ... `m10`.
fit_all_models <- function(zone_sets, k_w, type_min = P$type_min) {
  out <- pmap(MODEL_KEY, function(model, level, by_type, data_key, scope)
    fit_one_model(zone_sets[[data_key]], by_type, k_w, type_min))
  set_names(out, paste0('m', MODEL_KEY$model))
}

#' Long table of coefficients with standard errors, p-values, and group sample sizes.
tidy_coefs <- function(models) {
  imap_dfr(models, function(m, nm) {
    m$coefs %>%
      mutate(model = as.integer(sub('m', '', nm)), .before = 1)
  })
}

k_w_of <- function(m) m$k_w

#' Model performance as computed for Table 1.
#' R2, RMSE, and common-sample AIC/BIC are evaluated over ALL training rows:
#' rows covered by a modeled group use its fitted values; rows in groups with
#' <= 100 observations fall back to the global model of the same kind (Model 1
#' for no-type models, Model 2 for type models).
model_performance <- function(models, training) {
  full <- ungroup(training)
  kw <- models$m1$k_w
  c1 <- models$m1$coefs %>% select(term, estimate) %>% deframe()
  c2 <- models$m2$coefs %>% select(type, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate)

  pmap_dfr(MODEL_KEY, function(model, level, by_type, data_key, scope) {
    m <- models[[paste0('m', model)]]
    md <- ungroup(m$data)
    fitted_part <- md %>% mutate(kd_hat = m$fitted.values) %>% select(SiteID, date_unity, kd, kd_hat)
    if (data_key == 'global') {
      comb <- fitted_part
      sizes <- if (by_type) count(md, type)$n else nrow(md)
    } else {
      left <- anti_join(full, fitted_part, by = c('SiteID', 'date_unity'))
      if (by_type) {
        left <- left %>% left_join(c2, by = 'type') %>%
          mutate(kd_hat = tss_dead * NAP + doc * DOC + chl_a * Chl_a)
        sizes <- count(md, code, type)$n
      } else {
        left <- left %>%
          mutate(kd_hat = tss_dead * c1[['NAP']] + doc * c1[['DOC']] + chl_a * c1[['Chl_a']])
        sizes <- count(md, code)$n
      }
      comb <- bind_rows(fitted_part, left %>% select(SiteID, date_unity, kd, kd_hat))
    }
    N_all <- nrow(comb)
    rss_all <- sum((comb$kd - kw - comb$kd_hat)^2)
    p_eff <- m$n_coef + if (data_key != 'global' && N_all > nrow(md)) (if (by_type) 9L else 3L) else 0L
    ll_all <- -N_all / 2 * (log(2 * pi) + log(rss_all / N_all) + 1)
    neg_groups <- m$coefs %>% group_by(code, type) %>% summarize(any_neg = any(estimate <= 0), .groups = 'drop')

    tibble(
      model        = model,
      level        = level,
      by_type      = by_type,
      n_groups     = length(sizes),
      group_median = round(median(sizes)),
      group_min    = min(sizes),
      group_max    = max(sizes),
      r2           = cor(comb$kd, comb$kd_hat)^2,
      rmse_kd      = sqrt(rss_all / N_all),
      n_coef       = m$n_coef,
      n_neg_coefs  = sum(m$coefs$estimate <= 0),
      n_neg_groups = sum(neg_groups$any_neg),
      aic_common   = -2 * ll_all + 2 * (p_eff + 1),
      bic_common   = -2 * ll_all + log(N_all) * (p_eff + 1),
      glm_aic      = m$aic,
      glm_bic      = m$bic,
      glm_nobs     = m$nobs
    )
  })
}

#' Wide coefficient table in the layout of manuscript Table S1.
table_s1 <- function(models) {
  tidy_coefs(models) %>%
    select(model, code, type, n_group, term, estimate, std_error, p_value) %>%
    pivot_wider(names_from = term, values_from = c(estimate, std_error, p_value)) %>%
    select(model, code, type, n_group,
           NAP = estimate_NAP, DOC = estimate_DOC, Chl_a = estimate_Chl_a,
           se_NAP = std_error_NAP, se_DOC = std_error_DOC, se_Chl_a = std_error_Chl_a,
           p_NAP = p_value_NAP, p_DOC = p_value_DOC, p_Chl_a = p_value_Chl_a) %>%
    arrange(model, type, code)
}

# ---- Coefficient stability with random subsets (Fig S1 "cones") -------------

#' Exact OLS (no intercept) for many groups at once via grouped cross-products.
grouped_ols <- function(X, y, g) {
  g <- as.integer(factor(g))
  xx <- function(i, j) rowsum(X[, i] * X[, j], g, reorder = TRUE)[, 1]
  xy <- function(i) rowsum(X[, i] * y, g, reorder = TRUE)[, 1]
  S <- cbind(xx(1,1), xx(1,2), xx(1,3), xx(2,2), xx(2,3), xx(3,3))
  B <- cbind(xy(1), xy(2), xy(3))
  t(vapply(seq_len(nrow(S)), function(k) {
    A <- matrix(S[k, c(1,2,3, 2,4,5, 3,5,6)], 3, 3)
    tryCatch(solve(A, B[k, ]), error = function(e) rep(NA_real_, 3))
  }, numeric(3)))
}

#' Split the training data into G random groups (averaging n/50 ... n/10000 rows,
#' plus G = 1), fit the global model in every group, repeat `reps` times.
random_subset_coefs <- function(training, k_w, seed = P$random_seed) {
  set.seed(seed)
  n <- nrow(training)
  X <- as.matrix(training[, c('tss_dead', 'doc', 'chl_a')])
  y <- training$kd - k_w
  G <- rep(round(n / c(P$random_divisors, n)), P$random_reps)
  map_dfr(seq_along(G), function(i) {
    rg <- sample(seq_len(G[i]), size = n, replace = TRUE)
    b <- grouped_ols(X, y, rg)
    tibble(rep = i, rgroup = seq_len(nrow(b)),
           tss_dead = b[, 1], doc = b[, 2], chl_a = b[, 3]) %>%
      pivot_longer(c(tss_dead, doc, chl_a), names_to = 'term', values_to = 'estimate') %>%
      mutate(group_size = G[i], each_group = round(n / G[i]))
  })
}

#' Envelope used for the grey cones in Fig S1: 0.001% and 99.999% quantiles.
cone_envelope <- function(rand, nonneg_only = FALSE) {
  if (nonneg_only) rand <- filter(rand, estimate >= 0)
  rand %>%
    group_by(group_size, term, each_group) %>%
    summarize(q99 = quantile(estimate, .99999, na.rm = TRUE),
              q01 = quantile(estimate, .00001, na.rm = TRUE), .groups = 'drop')
}

# ---- Application of Model 7 (Figs 2-4) --------------------------------------

#' Ternary fractions and modeled Secchi depth for every application record.
#'  * Records in Level-3 ecoregions that (a) have Model 7 coefficients with all
#'    three slopes strictly positive and (b) have > apply_min_n application
#'    records use Model 7 coefficients; all others fall back to global Model 1.
#'  * fraction_i = a_i * C_i / sum_j(a_j * C_j)  (water term excluded)
#'  * secchi_pred = round(1 / (att_sum + k_w), 2) using the consistent training k_w.
apply_model7 <- function(application, eco, models, training,
                         require_positive = P$require_positive_coefs) {
  m7 <- models$m7; m1 <- models$m1
  kw <- m7$k_w
  app3 <- eco_join(application, eco$eco3, l3code) %>%
    select(all_of(names(application)), code = l3code, name = l3name)

  c7 <- m7$coefs %>%
    select(code, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate, names_prefix = 'coef_')
  if (require_positive) {
    c7 <- c7 %>% filter(coef_NAP > 0, coef_DOC > 0, coef_Chl_a > 0)
  }

  c1 <- m1$coefs %>% select(term, estimate) %>% deframe()
  matched <- app3 %>%
    inner_join(c7, by = 'code') %>%
    add_count(code) %>%
    filter(n > P$apply_min_n) %>%
    select(-n) %>%
    mutate(model = 'z3')
  unmatched <- anti_join(app3, matched, by = c('SiteID', 'date_unity')) %>%
    mutate(coef_NAP = c1[['NAP']], coef_DOC = c1[['DOC']], coef_Chl_a = c1[['Chl_a']], model = 'global')

  bind_rows(unmatched, matched) %>%
    left_join(ungroup(training) %>% select(SiteID, date_unity, kd), by = c('SiteID', 'date_unity')) %>%
    mutate(att_nap = coef_NAP * tss_dead,
           att_doc = coef_DOC * doc,
           att_chl = coef_Chl_a * chl_a,
           att_sum = att_nap + att_doc + att_chl,
           tss_kd  = att_nap / att_sum,
           doc_kd  = att_doc / att_sum,
           chl_kd  = att_chl / att_sum,
           kd_pred = att_sum + kw,
           secchi_pred = round(1 / kd_pred, 2),
           secchi_pred_consistent = 1 / kd_pred,
           secchi_status = ifelse(is.na(secchi), 'Modeled', 'Measured'))
}
