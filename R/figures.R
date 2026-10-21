# -----------------------------------------------------------------------------
# figures.R -- every manuscript figure, assembled entirely in R
# -----------------------------------------------------------------------------

# ---- Figure 1 ---------------------------------------------------------------

#' Fig 1a: distributions of LACs and Secchi, by water type, for records with
#' measured vs modeled Secchi (using consistent k_w).
fig1a <- function(fd) {
  fd %>%
    mutate(Secchi = secchi_status, secchi = ifelse(is.na(secchi), secchi_pred, secchi)) %>%
    filter(secchi > 0) %>%
    select(tss, secchi, doc, chl_a, type, Secchi) %>%
    pivot_longer(tss:chl_a, names_to = 'parameter', values_to = 'value') %>%
    filter(value >= 0.01) %>%
    mutate(parameter = factor(parameter, c('chl_a', 'doc', 'secchi', 'tss'),
                              c('Chl a (\u00b5g/L)', 'DOC (mg/L)', 'Secchi (m)', 'TSS (mg/L)'))) %>%
    ggplot(aes(value, colour = type, linetype = Secchi)) +
    geom_density(linewidth = 0.8) +
    facet_wrap(~parameter, scales = 'free', ncol = 4) +
    scale_x_log10(labels = label_number(drop0trailing = TRUE)) +
    scale_colour_manual(values = TYPE_COLS) +
    labs(x = NULL, y = 'density', colour = 'type') +
    theme_bw(base_size = 10) +
    theme(panel.grid = element_blank(), strip.background = element_blank(),
          legend.key.width = unit(1.2, 'lines'))
}

#' Fig 1b: map of observation density (0.1 degree bins) over Level III ecoregions.
fig1b <- function(fd, eco, seed = 1) {
  eco3 <- eco$eco3 %>% group_by(code = l3code, name = l3name) %>% summarize(.groups = 'drop') %>%
    st_simplify(dTolerance = 3000, preserveTopology = TRUE)
  set.seed(seed)
  greys <- sample(colorRampPalette(c('#f4f0ea', '#cfc4b8', '#a39484'))(nrow(eco3)))
  sites <- fd %>%
    mutate(lat = round(lat, 1), lon = round(long, 1)) %>%
    group_by(lat, lon, type, Secchi = secchi_status) %>% summarize(n = n(), .groups = 'drop') %>%
    filter(n > 1) %>% arrange(desc(Secchi)) %>%
    st_as_sf(coords = c('lon', 'lat'), crs = 4326) %>% st_transform(st_crs(eco3))
  ggplot() +
    geom_sf(data = eco3, aes(fill = name), colour = NA, show.legend = FALSE) +
    scale_fill_manual(values = greys) +
    geom_sf(data = sites, aes(colour = type, shape = Secchi, size = n), stroke = 0.6) +
    scale_shape_manual(values = c(Measured = 0, Modeled = 16)) +
    scale_colour_manual(values = TYPE_COLS) +
    scale_size_continuous(trans = 'log10', range = c(0.4, 4), breaks = c(10, 100, 1000),
                          labels = c('10^1', '10^2', '10^3')) +
    coord_sf(datum = NA, xlim = st_bbox(eco3)[c(1, 3)], ylim = st_bbox(eco3)[c(2, 4)]) +
    guides(size = guide_legend(override.aes = list(shape = 0))) +
    theme_void(base_size = 10) + theme(legend.position = 'right')
}

fig1 <- function(fd, eco) {
  (fig1a(fd) + theme(legend.position = 'right')) / fig1b(fd, eco) +
    plot_layout(heights = c(1, 3.4)) + plot_annotation(tag_levels = 'a', tag_suffix = ')')
}

# ---- Figure 2 ---------------------------------------------------------------

tern_summary_panel <- function(df, label, pt_alpha = 0.20, pt_size = 0.04) {
  r1 <- tern_sd_region(df, 1); r2 <- tern_sd_region(df, 2)
  mu <- df %>% summarize(across(c(chl_kd, doc_kd, tss_kd), mean)) %>% tern_xy()
  tern_frame(base_size = 10) +
    geom_point(data = tern_xy(df), aes(.x, .y), colour = 'grey55', size = pt_size, alpha = pt_alpha, shape = 16) +
    geom_polygon(data = r2, aes(.x, .y), fill = NA, colour = 'black', linewidth = 0.7) +
    geom_polygon(data = r1, aes(.x, .y), fill = NA, colour = 'red', linewidth = 0.7) +
    geom_point(data = mu, aes(.x, .y), colour = 'blue', size = 2.2) +
    annotate('text', x = -0.18, y = 0.9, hjust = 0, vjust = 1, size = 3.6,
             label = sprintf('%s\nn = %sK', label, floor(nrow(df) / 1000)))
}

#' Fig 2: ternary distributions with mean and +/-1, 2 SD regions.
fig2 <- function(fd, fd_types = fd) {
  pa <- tern_summary_panel(fd, 'a) All types')
  pb <- tern_summary_panel(filter(fd_types, type == 'Stream'), 'b) Streams')
  pc <- tern_summary_panel(filter(fd_types, type == 'Lake'), 'c) Lakes')
  pd <- tern_summary_panel(filter(fd_types, type == 'Estuary'), 'd) Estuaries')
  (pa | pb) / (pc | pd)
}

# ---- Figure 3: Delaware River longitudinal profile ---------------------------

#' Sites with long records (>= long_term_min_years distinct years).
long_term_sites <- function(fd) {
  fd %>% mutate(date = as.Date(date_unity)) %>%
    group_by(lat, long, SiteID) %>% distinct(date, .keep_all = TRUE) %>%
    mutate(year = year(date), n_years = n_distinct(year)) %>%
    filter(n_years >= P$long_term_min_years) %>% ungroup()
}

#' Distance from Trenton for every Delaware River Basin Commission mainstem site:
#' snap to RiverSR flowlines of the Delaware mainstem (path length < 158 km);
#' sites in the bay get cumulative straight-line distance from the previous site.
delaware_distances <- function(lt) {
  dr <- lt %>% filter(grepl('31DELRBC_WQX', SiteID)) %>%
    group_by(lat, long, SiteID, year, type) %>% filter(n() > 1) %>% ungroup()
  dr_sf <- dr %>% distinct(lat, long, SiteID, type) %>% arrange(-lat) %>%
    st_as_sf(coords = c('long', 'lat'), crs = 4326) %>% st_transform(2163)
  river <- st_read(file.path(PATH$raw, 'riverSR/nhdplusv2_modified_v1.0.shp'), quiet = TRUE,
                   query = "SELECT * FROM \"nhdplusv2_modified_v1.0\" WHERE TermnlP = 200005438") %>%
    filter(GNIS_NA == 'Delaware River', Pthlngt < 158) %>%
    mutate(dist_km = Pthlngt) %>% st_transform(2163) %>% select(dist_km)
  d <- dr_sf %>% st_join(river, join = st_nearest_feature) %>%
    mutate(dist_km = ifelse(dist_km == min(dist_km), NA, dist_km),
           dist_trenton = abs(dist_km - max(dist_km, na.rm = TRUE)))
  step <- c(0, as.numeric(st_distance(d$geometry[-nrow(d)], d$geometry[-1], by_element = TRUE)) / 1000)
  d %>% mutate(est_dist = ifelse(is.na(dist_trenton), step, 0),
               cume_euc = cumsum(est_dist) + max(dist_km, na.rm = TRUE),
               dist_trenton = ifelse(is.na(dist_trenton), cume_euc, dist_trenton)) %>%
    select(SiteID, type, dist_trenton) %>%
    list(sites = ., river = river)
}

DELAWARE_EXEMPLARS <- tibble(
  SiteID = c('31DELRBC_WQX-332068', '31DELRBC_WQX-892070', '31DELRBC_WQX-091011',
             '31DELRBC_WQX-091002', '31DELRBC_WQX-091020', '31DELRBC_WQX-091030'),
  label = c('1) Trenton, River', '2) 41 km, River', '3) 95 km, River',
            '4) 123 km, River', '5) 180 km, Estuary', '6) 231 km, Estuary'),
  num = 1:6)

fig3_data <- function(fd) {
  lt <- long_term_sites(fd)
  dd <- delaware_distances(lt)
  with_dist <- lt %>% inner_join(st_drop_geometry(dd$sites), by = c('SiteID', 'type'))
  list(
    site_avg = with_dist %>% group_by(SiteID, dist_trenton, type) %>%
      summarize(across(c(chl_kd, doc_kd, tss_kd), mean), .groups = 'drop'),
    year_avg = with_dist %>% group_by(SiteID, dist_trenton, type, year) %>%
      summarize(across(c(chl_kd, doc_kd, tss_kd), mean), n = n(), .groups = 'drop') %>%
      inner_join(DELAWARE_EXEMPLARS, by = 'SiteID') %>%
      mutate(site_label = factor(label, DELAWARE_EXEMPLARS$label)) %>% select(-label),
    sites = dd$sites, river = dd$river,
    years = range(with_dist$year[with_dist$SiteID %in% DELAWARE_EXEMPLARS$SiteID]))
}

basemap_tiles <- function(bbox_sf, zoom = 9, provider = 'Esri.WorldGrayCanvas') {
  cache <- file.path(PATH$derived, 'tiles'); dir.create(cache, showWarnings = FALSE)
  r <- maptiles::get_tiles(st_transform(bbox_sf, 3857), provider = provider, zoom = zoom,
                           crop = TRUE, cachedir = cache)
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3:5] <- c('r', 'g', 'b')
  df %>% mutate(fill = rgb(r, g, b, maxColorValue = 255))
}

fig3 <- function(f3) {
  lims <- c(0, max(f3$sites$dist_trenton))
  vir <- scale_colour_viridis_c(limits = lims, name = 'Distance from\nTrenton (km)')
  sites3857 <- st_transform(f3$sites, 3857)
  river3857 <- st_transform(f3$river, 3857)
  bb <- st_as_sfc(st_bbox(st_buffer(sites3857, 13000)))
  tiles <- basemap_tiles(bb)
  ex <- sites3857 %>% inner_join(DELAWARE_EXEMPLARS, by = 'SiteID')
  exy <- st_coordinates(ex)
  map <- ggplot() +
    geom_raster(data = tiles, aes(x, y, fill = fill)) + scale_fill_identity() +
    geom_sf(data = river3857, colour = 'blue', linewidth = 1.6) +
    geom_sf(data = sites3857, aes(fill = NULL, colour = dist_trenton), size = 4.2, alpha = 0.75) +
    geom_sf(data = sites3857, shape = 1, size = 4.2, colour = 'grey20') +
    vir +
    annotate('point', x = exy[, 1] + 9000, y = exy[, 2], shape = 21, size = 5.5, fill = 'white', colour = 'black') +
    annotate('text', x = exy[, 1] + 9000, y = exy[, 2], label = ex$num, size = 3.2) +
    coord_sf(crs = 3857, datum = NA, expand = FALSE) + theme_void(base_size = 9) +
    theme(legend.position = 'inside', legend.position.inside = c(0.14, 0.17),
          legend.key.height = unit(0.5, 'cm'), legend.key.width = unit(0.3, 'cm'),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6))
  inset <- tern_frame(base_size = 6) +
    geom_point(data = tern_xy(f3$site_avg), aes(.x, .y, colour = dist_trenton), size = 1.3) +
    vir + theme(panel.background = element_rect(fill = alpha('white', 0.9), colour = NA))
  map_panel <- map + inset_element(inset, left = 0.36, bottom = 0.26, right = 1.04, top = 0.66)
  facets <- tern_frame(base_size = 9) +
    geom_point(data = tern_xy(f3$year_avg), aes(.x, .y, colour = dist_trenton), size = 2.6) +
    vir + facet_wrap(~site_label, ncol = 3) +
    theme(strip.background = element_rect(fill = 'grey88', colour = 'grey40'),
          strip.text = element_text(size = 9, margin = margin(3, 3, 3, 3)))
  (map_panel | facets) + plot_layout(widths = c(1, 2.6))
}

# ---- Figure 4: temporal variability (ternary area fraction CDFs) -------------

FIG4_SITES <- tibble(
  SiteID = c('USGS-07373420', '21FLFTM_WQX-CALUSA0023FTM', 'MNPCA-29-0250-00-201',
             '21FLSJWM-LGNA', 'CEDEN-KR25444_KHSA'),
  inset = c(1, 1, 1, 2, 3),
  label = c('1S', '1E', '1L', '2L', '3S'),
  site_name = c('Mississippi R. near St. Francisville, LA',
                'Caloosahatchee Estuary (Orange R.), FL',
                'Portage Lake, MN',
                'Lake Griffin (north station), FL',
                'Klamath R. below Iron Gate Dam, CA'),
  role = c('lucistatic NAP', 'lucistatic DOC', 'lucistatic Chl a',
           'multi-annual trajectory', 'seasonal cycle'))

#' Per-site ternary area fraction (sites with > fig4_min_obs records) using the
#' exact Sutherland-Hodgman simplex clipping metric (`tern_area_fraction`).
fig4_areas <- function(fd) {
  fd %>% group_by(lat, long, SiteID, type) %>% mutate(n = n()) %>%
    filter(n > P$fig4_min_obs) %>%
    group_by(lat, long, SiteID, type, n) %>% nest() %>% ungroup() %>%
    mutate(area = map_dbl(data, tern_area_fraction, k = 2)) %>%
    select(-data)
}

season_of <- function(m) factor(case_when(m %in% c(12, 1, 2) ~ 'winter', m %in% 3:5 ~ 'spring',
                                          m %in% 6:8 ~ 'summer', TRUE ~ 'fall'),
                                c('winter', 'spring', 'summer', 'fall'))

circ_num <- function(p, n, size = 3.5) {
  p + annotate('point', x = 0.22, y = 0.83, shape = 21, size = size * 2.1, fill = 'white', stroke = 0.7) +
    annotate('text', x = 0.22, y = 0.83, label = n, size = size)
}

fig4 <- function(fd, areas) {
  ecdf_at <- function(ty, x) ecdf(areas$area[areas$type == ty])(x)
  labs_df <- FIG4_SITES %>%
    inner_join(areas %>% select(SiteID, type, area), by = 'SiteID') %>%
    mutate(cdf = map2_dbl(type, area, ecdf_at),
           label = paste0(inset, substr(type, 1, 1)))
  cdf <- ggplot(areas, aes(area, colour = type)) +
    stat_ecdf(linewidth = 0.95) +
    geom_point(data = labs_df, aes(area, cdf), shape = 21, size = 8.2, fill = 'white', stroke = 0.95, show.legend = FALSE) +
    geom_text(data = labs_df, aes(area, cdf, label = label), colour = 'black', size = 3.4, show.legend = FALSE) +
    scale_colour_manual(values = TYPE_COLS, name = NULL) +
    coord_cartesian(xlim = c(0, 0.78)) +
    labs(x = 'Ternary area fraction', y = 'CDF') + theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = 'inside',
          legend.position.inside = c(0.10, 0.86),
          legend.key.height = unit(0.42, 'cm'), legend.key.width = unit(0.65, 'cm'),
          legend.background = element_blank())
  s7 <- fd %>% group_by(lat, long, SiteID, type) %>% mutate(n = n()) %>%
    filter(n > P$fig4_s7_min_obs) %>% ungroup() %>%
    mutate(year = year(date_unity), month = month(date_unity), season = season_of(month))
  i1 <- tern_frame(base_size = 8.2, compact = TRUE) +
    geom_point(data = tern_xy(filter(fd, SiteID %in% FIG4_SITES$SiteID[1:3])), aes(.x, .y, colour = type), size = 1.35) +
    scale_colour_manual(values = TYPE_COLS)
  trend <- s7 %>% filter(SiteID == FIG4_SITES$SiteID[4]) %>% group_by(year) %>%
    summarize(across(c(chl_kd, doc_kd, tss_kd), mean))
  i2 <- tern_frame(base_size = 8.0, compact = TRUE) +
    geom_point(data = tern_xy(trend), aes(.x, .y, colour = year), size = 1.55) +
    scale_colour_viridis_c(name = NULL, breaks = function(l) unique(round(pretty(l, 4)))) +
    theme(legend.position = 'inside', legend.position.inside = c(0.92, 0.55),
          legend.key.height = unit(0.42, 'cm'), legend.key.width = unit(0.22, 'cm'),
          legend.text = element_text(size = 7.0), legend.background = element_blank())
  i3 <- tern_frame(base_size = 8.0, compact = TRUE) +
    geom_point(data = tern_xy(filter(s7, SiteID == FIG4_SITES$SiteID[5])), aes(.x, .y, colour = season), size = 1.35) +
    scale_colour_manual(values = c(winter = 'black', spring = 'lightgreen', summer = 'lightblue', fall = 'brown3'), name = NULL) +
    theme(legend.position = 'inside', legend.position.inside = c(0.91, 0.68),
          legend.key.height = unit(0.34, 'cm'), legend.key.width = unit(0.26, 'cm'),
          legend.text = element_text(size = 7.2), legend.background = element_blank())
  cdf +
    inset_element(circ_num(i1, 1), left = 0.25, bottom = 0.01, right = 0.61, top = 0.45) +
    inset_element(circ_num(i2, 2), left = 0.60, bottom = 0.01, right = 0.98, top = 0.45) +
    inset_element(circ_num(i3, 3), left = 0.57, bottom = 0.47, right = 0.95, top = 0.91)
}

# ---- Figure S1: coefficient coherence cones ----------------------------------

figS1 <- function(cones, coefs7, coefs8) {
  pts <- bind_rows(coefs7, coefs8) %>% filter(n_group > P$zone_min) %>%
    mutate(term = recode(term, NAP = 'tss_dead', DOC = 'doc', Chl_a = 'chl_a'))
  ggplot() +
    geom_ribbon(data = cones, aes(x = each_group, ymin = q01, ymax = q99), fill = 'grey60', alpha = 0.5) +
    geom_point(data = pts, aes(n_group, estimate, colour = code), size = 1.8) +
    facet_grid(term ~ type, scales = 'free_y') +
    scale_x_log10(labels = label_log(digits = 1)) +
    scale_colour_viridis_d(name = 'ecoregion', guide = 'none') +
    labs(x = 'n', y = 'coefficient') + theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), strip.background = element_blank())
}

# ---- QC Diagnostic Figure: Beer-Lambert physical plausibility screen ---------

fig_qc_discordance <- function(trimmed) {
  d <- trimmed %>%
    add_nap() %>%
    mutate(
      secchi = round(ifelse(secchi < P$secchi_floor, P$secchi_floor, secchi), P$secchi_digits),
      kd_obs = round(1 / secchi, P$kd_digits),
      kd = kd_obs,
      kd_lit = P$lit_kw + P$lit_alpha_nap * pmax(0, tss_dead) +
        P$lit_alpha_doc * doc + P$lit_alpha_chl * chl_a,
      log_ratio = log10(kd_obs / kd_lit)
    )
  n_ret <- sum(abs(d$log_ratio) <= P$max_log10_kd_discordance)
  pct_ret <- 100 * n_ret / nrow(d)
  ret_lbl <- sprintf('Retained (%.2f%%, n = %s)', pct_ret, fmt_n(n_ret))
  d <- d %>%
    mutate(
      status = case_when(
        log_ratio >  P$max_log10_kd_discordance ~ 'Removed: high kd (>6.3x; e.g. shallow bottom hits)',
        log_ratio < -P$max_log10_kd_discordance ~ 'Removed: low kd (<1/6.3x; e.g. unit/bottle errors)',
        TRUE ~ ret_lbl
      )
    )

  set.seed(42)
  d_sub <- bind_rows(
    d %>% filter(status == ret_lbl) %>% slice_sample(n = 35000),
    d %>% filter(status != ret_lbl)
  )

  col_vals <- setNames(
    c('grey45', 'firebrick', 'steelblue4'),
    c(ret_lbl,
      'Removed: high kd (>6.3x; e.g. shallow bottom hits)',
      'Removed: low kd (<1/6.3x; e.g. unit/bottle errors)')
  )

  pa <- ggplot(d_sub, aes(kd_lit, kd_obs, colour = status)) +
    geom_point(size = 0.55, alpha = 0.45, shape = 16) +
    geom_abline(slope = 1, intercept = 0, colour = 'black', linewidth = 0.7) +
    geom_abline(slope = 1, intercept =  P$max_log10_kd_discordance, colour = 'firebrick', linetype = 'dashed', linewidth = 0.65) +
    geom_abline(slope = 1, intercept = -P$max_log10_kd_discordance, colour = 'steelblue4', linetype = 'dashed', linewidth = 0.65) +
    scale_x_log10(labels = label_number(drop0trailing = TRUE)) +
    scale_y_log10(labels = label_number(drop0trailing = TRUE)) +
    scale_colour_manual(values = col_vals, name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 0.9))) +
    labs(
      title = 'a) A priori Beer-Lambert physical plausibility screen',
      x = expression(Canonical~literature~k[d*","*lit]~(m^-1)),
      y = expression(Observed~k[d] == 1/z[SD]~(m^-1))
    ) +
    theme_bw(base_size = 10.5) +
    theme(panel.grid = element_blank(), legend.position = 'bottom',
          legend.direction = 'vertical', legend.margin = margin(0, 0, 0, 0))

  d_eco3 <- eco_join(d, get_eco()$eco3, l3code) %>% mutate(code = l3code) %>% select(-any_of(c('l3code', 'l3name')))
  taus <- c(0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, Inf)
  sens <- map_dfr(taus, function(tau) {
    sub <- d_eco3 %>% filter(abs(log_ratio) <= tau)
    kw <- 1 / P$k_w_divisor / max(sub$secchi)
    fit1 <- lm.fit(cbind(sub$tss_dead, sub$doc, sub$chl_a), sub$kd - kw)
    pred1 <- fit1$fitted.values + kw
    sub_m7 <- sub %>% add_count(code) %>% filter(n > P$zone_min) %>% select(-n)
    m7_sub <- fit_one_model(sub_m7, by_type = FALSE, k_w = kw)
    fitted_m7 <- sub_m7 %>% mutate(kd_hat = m7_sub$fitted.values) %>% select(SiteID, date_unity, kd, kd_hat)
    left_m7 <- anti_join(sub, fitted_m7, by = c('SiteID', 'date_unity')) %>%
      mutate(kd_hat = tss_dead * fit1$coefficients[1] + doc * fit1$coefficients[2] + chl_a * fit1$coefficients[3]) %>%
      select(SiteID, date_unity, kd, kd_hat)
    comb7 <- bind_rows(fitted_m7, left_m7)
    tibble(
      tau = ifelse(is.infinite(tau), 1.45, tau),
      label = ifelse(is.infinite(tau), 'Unscreened', sprintf('%.2f', tau)),
      pct_removed = 100 * (1 - nrow(sub) / nrow(d_eco3)),
      R2_Model7 = cor(comb7$kd, comb7$kd_hat)^2,
      R2_Model1 = cor(sub$kd, pred1)^2,
      alpha_DOC = fit1$coefficients[2],
      alpha_NAP = fit1$coefficients[1],
      alpha_Chl10 = fit1$coefficients[3] * 10
    )
  })

  x_lbls <- map2_chr(sens$tau, sens$pct_removed, function(t_val, pct) {
    if (t_val > 1.3) sprintf('None\n(0%%)')
    else if (abs(t_val - P$max_log10_kd_discordance) < 1e-6) sprintf('%.1f*\n(%.1f%%)', t_val, pct)
    else sprintf('%.1f\n(%.1f%%)', t_val, pct)
  })

  sens_long <- sens %>%
    pivot_longer(c(R2_Model7, R2_Model1, alpha_DOC, alpha_NAP, alpha_Chl10), names_to = 'metric', values_to = 'val') %>%
    mutate(metric = factor(recode(metric,
      R2_Model7   = 'Model 7 (Level III Ecoregion) R²',
      R2_Model1   = 'Model 1 (Continental Global) R²',
      alpha_DOC   = 'α_DOC (m²/g)',
      alpha_NAP   = 'α_NAP (m²/g)',
      alpha_Chl10 = '10 × α_Chl a (m²/mg)'
    ), levels = c('Model 7 (Level III Ecoregion) R²', 'Model 1 (Continental Global) R²',
                  '10 × α_Chl a (m²/mg)', 'α_DOC (m²/g)', 'α_NAP (m²/g)')))

  pb <- ggplot(sens_long, aes(tau, val, colour = metric, group = metric)) +
    geom_vline(xintercept = P$max_log10_kd_discordance, linetype = 'dashed', colour = 'firebrick', linewidth = 0.65) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.2) +
    scale_x_continuous(
      breaks = sens$tau,
      labels = x_lbls
    ) +
    scale_colour_manual(
      values = c('Model 7 (Level III Ecoregion) R²' = '#00441b',
                 'Model 1 (Continental Global) R²' = '#5aae61',
                 '10 × α_Chl a (m²/mg)' = '#762a83',
                 'α_DOC (m²/g)' = '#2166ac',
                 'α_NAP (m²/g)' = '#d6604d'),
      name = NULL
    ) +
    labs(
      title = 'b) Sensitivity of R² and fitted slopes to threshold τ',
      x = 'Plausibility threshold τ = |log10(kd / kd,lit)| (% removed)',
      y = 'Metric / coefficient value'
    ) +
    theme_bw(base_size = 10.5) +
    theme(panel.grid = element_blank(), legend.position = 'bottom',
          legend.direction = 'vertical', legend.margin = margin(0, 0, 0, 0))

  pa | pb
}
