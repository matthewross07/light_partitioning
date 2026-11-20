# -----------------------------------------------------------------------------
# pipeline.R -- cached "getters" for every pipeline product
#
# Each get_*() returns a product, building it (and its upstream inputs) only if
# needed. Slow products are cached on disk in data/derived/ (delete a file to
# force a rebuild); fast or large-in-memory products are memoised per session.
# -----------------------------------------------------------------------------

.memo <- new.env(parent = emptyenv())
memo <- function(key, expr) {
  if (!exists(key, envir = .memo, inherits = FALSE)) assign(key, force(expr), envir = .memo)
  get(key, envir = .memo, inherits = FALSE)
}

get_aquamatch <- function() memo('aquamatch', cache_rds('aquamatch_joined', build_aquamatch_joined()))
get_sites     <- function() memo('sites', cache_rds('sites', build_site_inventory(get_aquamatch())))

get_simul <- function(which = c('simul', 'no_secchi')) {
  which <- match.arg(which)
  f <- file.path(PATH$derived, paste0(c('simul', 'no_secchi'), '.rds'))
  if (!all(file.exists(f))) {
    am <- get_aquamatch()
    saveRDS(build_simultaneous(am, require_secchi = TRUE), f[1])
    saveRDS(build_simultaneous(am, require_secchi = FALSE), f[2])
  }
  memo(which, readRDS(file.path(PATH$derived, paste0(which, '.rds'))))
}

get_censor_fits <- function() memo('censor_fits', censor_fits(get_simul('no_secchi')))
get_trimmed     <- function() memo('trimmed', trim_simultaneous(get_simul('simul')))
get_training    <- function() memo('training', build_training(get_trimmed()))
#' Counts for the data screen: matchups, kept, flagged, and the two tails.
get_screen_counts <- function() memo('screen_counts', {
  t <- get_trimmed() %>% add_nap() %>%
    mutate(secchi = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits), kd = round(1 / secchi, P$kd_digits),
           r = log10(kd / (P$lit_kw + P$lit_alpha_nap * pmax(0, tss_dead) + P$lit_alpha_doc * doc + P$lit_alpha_chl * chl_a)),
           z = (r - mean(r)) / sd(r))
  list(n = nrow(t), kept = sum(abs(t$z) <= P$kd_discordance_sd), flagged = sum(abs(t$z) > P$kd_discordance_sd),
       upper = sum(t$z > P$kd_discordance_sd), lower = sum(t$z < -P$kd_discordance_sd))
})
get_k_w         <- function() k_w_from_training(get_training())
get_application <- function() memo('application', build_application(get_simul('no_secchi'), get_training(), kd_discordance_stats(get_trimmed())))
get_eco         <- function() memo('eco', cache_rds('eco', read_ecoregions()))
get_zone_sets   <- function() memo('zone_sets', cache_rds('zone_sets', build_zone_sets(get_training(), get_eco())))

get_models <- function(type_min = P$type_min) {
  memo(paste0('models_', type_min), fit_all_models(get_zone_sets(), get_k_w(), type_min))
}

get_final_data <- function() {
  memo('final_data', cache_rds('final_data',
       apply_model7(get_application(), get_eco(), get_models(), get_training())))
}

get_random_subsets <- function() memo('random', cache_rds('random_subsets', random_subset_coefs(get_training(), get_k_w())))
get_fig4_areas     <- function() memo('areas', cache_rds('fig4_areas', fig4_areas(get_final_data())))
get_fig3_data      <- function() memo('fig3', cache_rds('fig3_data', fig3_data(get_final_data())))
get_algal_sensitivity <- function() memo('algal_sens', cache_rds('algal_mass_sensitivity',
       algal_mass_sensitivity(get_training(), get_application(), get_eco(), get_k_w())))

#' Save a figure as PNG (300 dpi) and PDF into figs/.
save_fig <- function(p, name, width, height) {
  ggsave(file.path(PATH$figs, paste0(name, '.png')), p, width = width, height = height, dpi = 300, bg = 'white')
  ggsave(file.path(PATH$figs, paste0(name, '.pdf')), p, width = width, height = height, bg = 'white',
         device = grDevices::cairo_pdf)
  invisible(file.path(PATH$figs, paste0(name, '.png')))
}
