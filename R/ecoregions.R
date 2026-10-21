# -----------------------------------------------------------------------------
# ecoregions.R -- assigning observations to EPA/CEC ecoregions (Levels 1-4)
# -----------------------------------------------------------------------------

#' Read the three ecoregion shapefiles (Omernik & Griffith 2014).
#'   Level 1 and 2 : CEC North America Level II layer (has L1 + L2 attributes)
#'   Level 3       : US Level III layer; code = NA_L3CODE (e.g. '8.3.5')
#'   Level 4       : US Level IV layer (no state boundaries); code = US_L4CODE ('63a')
#' 'WATER' polygons are dropped from L2/L3.
read_ecoregions <- function() {
  eco <- file.path(PATH$raw, 'ecoregions')
  list(
    eco2 = st_read(file.path(eco, 'na_cec_eco_l2/NA_CEC_Eco_Level2.shp'), quiet = TRUE) %>%
      select(code = NA_L2CODE, name = NA_L2NAME, l1code = NA_L1CODE, l1name = NA_L1NAME) %>%
      filter(name != 'WATER'),
    eco3 = st_read(file.path(eco, 'us_eco_l3/us_eco_l3.shp'), quiet = TRUE) %>%
      select(l3code = NA_L3CODE, l3name = NA_L3NAME, code = NA_L2CODE, name = NA_L2NAME) %>%
      filter(name != 'WATER'),
    eco4 = st_read(file.path(eco, 'us_eco_l4/us_eco_l4_no_st.shp'), quiet = TRUE) %>%
      select(l3code = NA_L3CODE, l3name = NA_L3NAME, code = NA_L2CODE, name = NA_L2NAME,
             l4code = US_L4CODE, l4name = US_L4NAME)
  )
}

#' Join unique sites to one ecoregion layer with 1-to-1 deduplication on polygon
#' boundaries and nearest-polygon fallback for offshore/coastal sites.
join_sites_to_layer <- function(sites_df, eco_reg, code_col) {
  sf_pts <- sites_df %>%
    distinct(SiteID, lat, long) %>%
    st_as_sf(coords = c('long', 'lat'), crs = 4326) %>%
    st_transform(st_crs(eco_reg))
  pass <- st_join(sf_pts, eco_reg) %>% distinct(SiteID, .keep_all = TRUE)
  nas <- pass %>%
    filter(is.na(.data[[code_col]])) %>%
    select(SiteID) %>%
    st_join(eco_reg, join = st_nearest_feature) %>%
    distinct(SiteID, .keep_all = TRUE)
  bind_rows(filter(pass, !is.na(.data[[code_col]])), nas) %>%
    distinct(SiteID, .keep_all = TRUE) %>%
    st_drop_geometry()
}

#' Build or read the cached site-to-ecoregion lookup across all three shapefile layers.
get_site_eco_lookup <- function(sites_df, eco = get_eco()) {
  f <- file.path(PATH$derived, 'site_eco_lookup.rds')
  if (file.exists(f)) {
    lkp <- readRDS(f)
    missing_sites <- sites_df %>% distinct(SiteID, lat, long) %>% anti_join(lkp$eco3, by = 'SiteID')
    if (nrow(missing_sites) == 0) return(lkp)
    lkp$eco2 <- bind_rows(lkp$eco2, join_sites_to_layer(missing_sites, eco$eco2, 'code')) %>% distinct(SiteID, .keep_all = TRUE)
    lkp$eco3 <- bind_rows(lkp$eco3, join_sites_to_layer(missing_sites, eco$eco3, 'l3code')) %>% distinct(SiteID, .keep_all = TRUE)
    lkp$eco4 <- bind_rows(lkp$eco4, join_sites_to_layer(missing_sites, eco$eco4, 'l4code')) %>% distinct(SiteID, .keep_all = TRUE)
    saveRDS(lkp, f)
    return(lkp)
  }
  lkp <- list(
    eco2 = join_sites_to_layer(sites_df, eco$eco2, 'code'),
    eco3 = join_sites_to_layer(sites_df, eco$eco3, 'l3code'),
    eco4 = join_sites_to_layer(sites_df, eco$eco4, 'l4code')
  )
  saveRDS(lkp, f)
  lkp
}

#' Point-in-polygon join of observation sites to an ecoregion layer (1-to-1 per SiteID).
eco_join <- function(df, eco_reg, code) {
  code_col <- rlang::as_name(rlang::enquo(code))
  layer_key <- case_when(
    code_col == 'l4code' ~ 'eco4',
    code_col == 'l3code' ~ 'eco3',
    TRUE ~ 'eco2'
  )
  lkp <- get_site_eco_lookup(df)
  sites <- lkp[[layer_key]]
  inner_join(sites, df, by = 'SiteID')
}

#' Build the five model datasets used in Table 1:
#'   global (no ecoregion), and Level 1-4 ecoregions with > zone_min obs each.
build_zone_sets <- function(training, eco, zone_min = P$zone_min) {
  keep <- names(training)
  j2 <- eco_join(training, eco$eco2, code)
  j3 <- eco_join(training, eco$eco3, l3code)
  j4 <- eco_join(training, eco$eco4, l4code)
  list(
    global = training,
    L1 = j2 %>% group_by(l1name, l1code) %>% add_count() %>% filter(n > zone_min) %>%
      ungroup() %>% select(all_of(keep), code = l1code, name = l1name),
    L2 = j2 %>% group_by(name, code) %>% add_count() %>% filter(n > zone_min) %>%
      ungroup() %>% select(all_of(keep), code = code, name = name),
    L3 = j3 %>% group_by(l3code, l3name) %>% add_count() %>% filter(n > zone_min) %>%
      ungroup() %>% select(all_of(keep), code = l3code, name = l3name),
    L4 = j4 %>% group_by(l4code, l4name) %>% add_count() %>% filter(n > zone_min) %>%
      ungroup() %>% select(all_of(keep), code = l4code, name = l4name)
  )
}

#' Ecoregion counts before the zone_min filter (for documenting the threshold).
zone_counts <- function(training, eco) {
  j2 <- eco_join(training, eco$eco2, code)
  tally_by <- function(df, code_col, name_col, level) {
    df %>% group_by(code = .data[[code_col]], name = .data[[name_col]]) %>%
      summarize(n = n(), .groups = 'drop') %>% mutate(level = level)
  }
  bind_rows(
    tally_by(j2, 'l1code', 'l1name', 1),
    tally_by(j2, 'code', 'name', 2),
    tally_by(eco_join(training, eco$eco3, l3code), 'l3code', 'l3name', 3),
    tally_by(eco_join(training, eco$eco4, l4code), 'l4code', 'l4name', 4)
  )
}
