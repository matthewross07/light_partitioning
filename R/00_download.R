# -----------------------------------------------------------------------------
# 00_download.R -- Step 0: Download harmonized AquaMatch data via {AquaMatchr}
#                  and spatial reference layers (EPA/CEC Ecoregions & RiverSR)
#
# Uses AquaMatchr::download_parameters() to retrieve the harmonized in-situ WQP
# products from the Environmental Data Initiative (EDI):
#   - sdd  (Secchi disk depth,                 edi.1856.1; Brousil et al. 2025)
#   - chla (Chlorophyll a,                     edi.1756.2; Brousil et al. 2024)
#   - doc  (Dissolved organic carbon,          edi.1809.1; Brousil et al. 2024)
#   - tss  (Total suspended solids & SSC,      edi.2048.2; De La Torre et al. 2025)
#
# Note on EDI PASTA authentication (July 2026+):
#   EDI now requires an EDI_API_KEY header on pasta.lternet.edu REST endpoints
#   (used internally by EDIutils 2.x). When EDI_API_KEY is not set in the R
#   environment, `enable_edirepository_gmn_fallback()` transparently routes the
#   four EDIutils calls inside `AquaMatchr::download_parameters()` through EDI's
#   public DataONE Generic Member Node (`gmn.edirepository.org`), so
#   `AquaMatchr::download_parameters()` works out of the box without a key.
# -----------------------------------------------------------------------------

AQUAMATCH_EDI_MANIFEST <- dplyr::tribble(
  ~param, ~package_id,   ~revision, ~entity_name,            ~entity_id,                         ~bytes,     ~citation,
  "sdd",  "edi.1856.1",  1L,        "sdd_harmonized_final",  "723fe7603c501a8e68b0f523a",        635344442,  "Brousil, M.R., Meyer, M.F., Willi, K., Steele, B.G., Topp, S.N., & Ross, M.R.V. (2025). AquaMatch Secchi Disk Depth from Water Quality Portal ~1970-2024 (edi.1856.1). Environmental Data Initiative.",
  "chla", "edi.1756.2",  2L,        "chla_harmonized_final", "c97458622bc2126fa6d6c84e183101cf", 815049179,  "Brousil, M.R., Meyer, M.F., Willi, K., Steele, B.G., Topp, S.N., & Ross, M.R.V. (2024). AquaMatch Chlorophyll a from Water Quality Portal ~1970-2024 (edi.1756.2). Environmental Data Initiative.",
  "doc",  "edi.1809.1",  1L,        "doc_harmonized_final",  "c50dc72b8ea5cb6605abb613b105429a", 419874195,  "Brousil, M.R., Meyer, M.F., Willi, K., Steele, B.G., Topp, S.N., & Ross, M.R.V. (2024). AquaMatch Dissolved Organic Carbon from Water Quality Portal ~1970-2024 (edi.1809.1). Environmental Data Initiative.",
  "tss",  "edi.2048.2",  2L,        "tss_harmonized_final",  "b258b13368563cdc63113e2a972a42db", 1030574938, "De La Torre, J., Brousil, M.R., Meyer, M.F., Willi, K., Steele, B.G., & Ross, M.R.V. (2025). AquaMatch Total Suspended Solids and Suspended Sediment Concentration from Water Quality Portal ~1970-2025 (edi.2048.2). Environmental Data Initiative."
)

#' Shim EDIutils when EDI_API_KEY is not set so AquaMatchr::download_parameters()
#' works without hitting HTTP 403 on pasta.lternet.edu.
enable_edirepository_gmn_fallback <- function() {
  if (nzchar(Sys.getenv("EDI_API_KEY"))) return(invisible(FALSE))
  if (!requireNamespace("EDIutils", quietly = TRUE)) return(invisible(FALSE))

  shim_revisions <- function(scope, identifier, filter = NULL, ...) {
    row <- AQUAMATCH_EDI_MANIFEST[grepl(paste0("^", scope, "\\.", identifier, "\\."),
                                        AQUAMATCH_EDI_MANIFEST$package_id), ]
    if (nrow(row) > 0) return(row$revision[1])
    1L
  }
  shim_citation <- function(packageId, ...) {
    row <- AQUAMATCH_EDI_MANIFEST[AQUAMATCH_EDI_MANIFEST$package_id == packageId, ]
    if (nrow(row) > 0) return(row$citation[1])
    paste0("AquaMatch Harmonized Data Product (", packageId, "), Environmental Data Initiative (EDI).")
  }
  shim_entities <- function(packageId, ...) {
    row <- AQUAMATCH_EDI_MANIFEST[AQUAMATCH_EDI_MANIFEST$package_id == packageId, ]
    tibble::tibble(entityId = row$entity_id, entityName = row$entity_name)
  }
  shim_read_entity <- function(packageId, entityId, ...) {
    parts <- strsplit(packageId, "\\.")[[1]]
    pasta_url <- sprintf("https://pasta.lternet.edu/package/data/eml/%s/%s/%s/%s",
                         parts[1], parts[2], parts[3], entityId)
    gmn_url <- paste0("https://gmn.edirepository.org/mn/v2/object/",
                      utils::URLencode(pasta_url, reserved = TRUE))
    h <- curl::new_handle(timeout = 1800, connecttimeout = 60)
    con <- curl::curl(gmn_url, handle = h)
    on.exit(try(close(con), silent = TRUE))
    readr::read_file_raw(con)
  }

  utils::assignInNamespace("list_data_package_revisions", shim_revisions,   ns = "EDIutils")
  utils::assignInNamespace("read_data_package_citation",  shim_citation,    ns = "EDIutils")
  utils::assignInNamespace("read_data_entity_names",      shim_entities,    ns = "EDIutils")
  utils::assignInNamespace("read_data_entity",            shim_read_entity, ns = "EDIutils")
  invisible(TRUE)
}

#' Download AquaMatch harmonized parameter datasets using AquaMatchr::download_parameters()
#' and cache them locally in `data/raw/aquamatch/<param>_harmonized_final.csv`.
download_aquamatch_data <- function(parameters = c("sdd", "chla", "doc", "tss"),
                                    dest_dir = file.path(PATH$raw, "aquamatch"),
                                    force = FALSE) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  enable_edirepository_gmn_fallback()

  for (p in parameters) {
    csv_path <- file.path(dest_dir, paste0(p, "_harmonized_final.csv"))
    if (file.exists(csv_path) && file.size(csv_path) > 1e8 && !force) {
      next
    }
    message("Downloading AquaMatch parameter '", p, "' via AquaMatchr::download_parameters()...")
    dl <- AquaMatchr::download_parameters(parameters = p, version = "newest")
    readr::write_csv(dl[[p]], csv_path)
    rm(dl); gc(verbose = FALSE)
  }
  invisible(file.path(dest_dir, paste0(parameters, "_harmonized_final.csv")))
}

#' Registry of raw data sources used by the pipeline (AquaMatch EDI packages,
#' EPA/CEC ecoregion shapefiles, and RiverSR flowlines).
data_sources <- function() {
  tribble(
    ~id, ~dest, ~url, ~bytes, ~used_for,
    # --- In-situ water quality: AquaMatch harmonized WQP via {AquaMatchr} (EDI)
    'aquamatch_sdd', 'aquamatch/sdd_harmonized_final.csv',
    'https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1856.1', 635344442,
    'AquaMatch harmonized Secchi disk depth (sdd, edi.1856.1; 2.74M records, 1970-2024) via AquaMatchr::download_parameters()',
    'aquamatch_chla', 'aquamatch/chla_harmonized_final.csv',
    'https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1756.2', 815049179,
    'AquaMatch harmonized Chlorophyll a (chla, edi.1756.2; 3.39M records, 1970-2024) via AquaMatchr::download_parameters()',
    'aquamatch_doc', 'aquamatch/doc_harmonized_final.csv',
    'https://portal.edirepository.org/nis/mapbrowse?packageid=edi.1809.1', 419874195,
    'AquaMatch harmonized Dissolved Organic Carbon (doc, edi.1809.1; 1.84M records, 1970-2024) via AquaMatchr::download_parameters()',
    'aquamatch_tss', 'aquamatch/tss_harmonized_final.csv',
    'https://portal.edirepository.org/nis/mapbrowse?packageid=edi.2048.2', 1030574938,
    'AquaMatch harmonized Total Suspended Solids & SSC (tss, edi.2048.2; 4.61M records, 1970-2025) via AquaMatchr::download_parameters()',
    # --- US EPA / CEC ecoregions (Omernik & Griffith 2014) ------------------
    'eco_l2', 'ecoregions/na_cec_eco_l2.zip',
    'https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l2.zip', 31388641,
    'CEC North America Level I and II ecoregions (Omernik & Griffith 2014)',
    'eco_l3', 'ecoregions/us_eco_l3.zip',
    'https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l3.zip', 28424315,
    'US Level III ecoregions (NA_L3CODE used as the Level 3 code)',
    'eco_l4', 'ecoregions/us_eco_l4.zip',
    'https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/us/us_eco_l4.zip', 66167640,
    'US Level IV ecoregions (no state boundaries)',
    # --- RiverSR / modified NHDPlusV2 flowlines (Gardner et al. 2021, Zenodo 4304567)
    'riversr_shp', 'riverSR/nhdplusv2_modified_v1.0.shp',
    'https://zenodo.org/api/records/4304567/files/nhdplusv2_modified_v1.0.shp/content', 16987912,
    'River centerlines with path length (Delaware longitudinal distances, Fig 3)',
    'riversr_shx', 'riverSR/nhdplusv2_modified_v1.0.shx',
    'https://zenodo.org/api/records/4304567/files/nhdplusv2_modified_v1.0.shx/content', 271652, '',
    'riversr_dbf', 'riverSR/nhdplusv2_modified_v1.0.dbf',
    'https://zenodo.org/api/records/4304567/files/nhdplusv2_modified_v1.0.dbf/content', 132318322, '',
    'riversr_prj', 'riverSR/nhdplusv2_modified_v1.0.prj',
    'https://zenodo.org/api/records/4304567/files/nhdplusv2_modified_v1.0.prj/content', 143, ''
  )
}

#' Download one source if it is not already present.
fetch_source <- function(id, sources = data_sources(), timeout = 3600) {
  if (grepl('^aquamatch_', id)) {
    p <- sub('^aquamatch_', '', id)
    return(download_aquamatch_data(parameters = p))
  }
  s <- sources[sources$id == id, ]
  dest <- file.path(PATH$raw, s$dest)
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(dest) && file.size(dest) == s$bytes) return(invisible(dest))
  message('Downloading ', id, ' -> ', dest)
  old <- options(timeout = timeout); on.exit(options(old))
  utils::download.file(s$url, destfile = dest, mode = 'wb', quiet = TRUE)
  if (file.size(dest) != s$bytes)
    warning(id, ': size ', file.size(dest), ' != expected ', s$bytes,
            ' (source may have been updated)')
  invisible(dest)
}

#' Fetch all AquaMatch datasets via {AquaMatchr} and all spatial layers.
fetch_all <- function() {
  download_aquamatch_data()
  src <- data_sources() %>% filter(!grepl('^aquamatch_', id))
  walk(src$id, fetch_source, sources = src)
  eco <- file.path(PATH$raw, 'ecoregions')
  for (z in c('na_cec_eco_l2', 'us_eco_l3', 'us_eco_l4')) {
    if (!dir.exists(file.path(eco, z))) unzip(file.path(eco, paste0(z, '.zip')), exdir = file.path(eco, z))
  }
  invisible(TRUE)
}

#' Status table of raw inputs (for Chapter 1).
raw_status <- function() {
  data_sources() %>%
    mutate(path = file.path(PATH$raw, dest),
           present = file.exists(path),
           size_ok = present & file.size(path) == bytes) %>%
    select(id, dest, bytes, present, size_ok, url, used_for)
}
