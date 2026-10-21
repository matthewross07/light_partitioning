# -----------------------------------------------------------------------------
# sonicray_lake.R -- Native {sonicray} 3D Water Volume (`sr_upload_mesh(water=...)`)
#                    with Smooth 3-Hill Lake Basin, Cutaway Glass Front,
#                    Submerged 3D Fish & Secchi Spheres, and Single Dynamic
#                    Ternary Diagram transitioning from Clear -> Dark Water.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sonicray)
  library(tidyverse)
  library(png)
})

smooth01 <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  x * x * (3 - 2 * x)
}

smoothstep_sr <- function(e0, e1, x) {
  smooth01((x - e0) / (e1 - e0))
}

#' Build a 3D streamlined trout/fish mesh to append via `sonicray::sr_mesh_append()`
#' Colored via `sonicray::sr_paint_albedo()` at `uv_body`, `uv_stripe`, and `uv_back`.
build_sr_fish_mesh <- function(center, length = 1.15, yaw = 0.0,
                               uv_back = c(0.010, 0.990),
                               uv_stripe = c(0.025, 0.990),
                               uv_belly = c(0.040, 0.990)) {
  n_len <- 20L
  n_ang <- 14L
  s_seq <- seq(0, 1, length.out = n_len)
  th_seq <- seq(0, 2 * pi, length.out = n_ang + 1L)[-(n_ang + 1L)]

  r_y <- ifelse(
    s_seq < 0.76,
    0.20 * sin(pi * (s_seq / 0.76)^0.72),
    0.048 + 0.17 * ((s_seq - 0.76) / 0.24)^1.1
  )
  r_z <- ifelse(
    s_seq < 0.76,
    0.12 * sin(pi * (s_seq / 0.76)^0.68),
    0.022 * (1 - 0.65 * (s_seq - 0.76) / 0.24)
  )
  dorsal <- ifelse(s_seq >= 0.34 & s_seq <= 0.56, 0.08 * sin(pi * (s_seq - 0.34) / 0.22), 0.0)

  verts <- matrix(0.0, nrow = n_len * n_ang + 2L, ncol = 3L)
  uvs   <- matrix(0.0, nrow = n_len * n_ang + 2L, ncol = 2L)
  idx_m <- matrix(0L, nrow = n_len, ncol = n_ang)

  k <- 1L
  for (i in seq_len(n_len)) {
    s <- s_seq[i]
    lx <- (0.5 - s) * length
    for (j in seq_len(n_ang)) {
      th <- th_seq[j]
      sy <- sin(th)
      cz <- cos(th)
      ly <- r_y[i] * length * sy + (if (sy > 0.6) dorsal[i] * length * ((sy - 0.6) / 0.4)^1.3 else 0.0)
      lz <- r_z[i] * length * cz
      verts[k, ] <- c(lx, ly, lz)
      uvs[k, ] <- if (sy > 0.28) uv_back else if (sy > -0.18) uv_stripe else uv_belly
      idx_m[i, j] <- k
      k <- k + 1L
    }
  }
  i_nose <- k; verts[i_nose, ] <- c(0.51 * length, 0, 0); uvs[i_nose, ] <- uv_back; k <- k + 1L
  i_tail <- k; verts[i_tail, ] <- c(-0.48 * length, 0, 0); uvs[i_tail, ] <- uv_stripe

  co <- cos(yaw); si <- sin(yaw)
  vx <- verts[, 1] * co - verts[, 3] * si + center[1]
  vy <- verts[, 2] + center[2]
  vz <- verts[, 1] * si + verts[, 3] * co + center[3]
  verts <- cbind(vx, vy, vz)

  tris <- list()
  for (i in seq_len(n_len - 1L)) {
    for (j in seq_len(n_ang)) {
      jn <- if (j == n_ang) 1L else j + 1L
      a <- idx_m[i, j]; b <- idx_m[i, jn]; c <- idx_m[i + 1L, j]; d <- idx_m[i + 1L, jn]
      tris[[length(tris) + 1L]] <- rbind(c(a, c, b), c(b, c, d))
    }
  }
  for (j in seq_len(n_ang)) {
    jn <- if (j == n_ang) 1L else j + 1L
    tris[[length(tris) + 1L]] <- rbind(
      c(i_nose, idx_m[1L, j], idx_m[1L, jn]),
      c(idx_m[n_len, j], i_tail, idx_m[n_len, jn])
    )
  }
  tri <- do.call(rbind, tris) - 1L
  storage.mode(tri) <- "integer"
  cpp_normals <- getFromNamespace("cpp_vertex_normals", "sonicray")
  list(vertices = unname(verts), normals = unname(cpp_normals(verts, tri)), uvs = unname(uvs), indices = tri, mesh_dims = c(0L, 0L))
}

#' Build the smooth 3-hill lake basin scene + submerged fish & Secchi spheres using {sonicray}
build_sonicray_water_scene <- function(n = 200L, tex = 1024L) {
  set.seed(42L)
  fbm  <- getFromNamespace("fbm_noise", "sonicray")
  rmat <- getFromNamespace("resize_mat", "sonicray")

  xs <- seq(-10, 10, length.out = n)
  zs <- seq(-10, 10, length.out = n) # row 1 = -10 (far back), row n = +10 (front cutaway edge)
  X <- matrix(rep(xs, each = n), n, n)
  Z <- matrix(rep(zs, times = n), n, n)

  # Smooth bowl lakebed opening toward the front-center cutaway wall (Z > 0),
  # surrounded by 3 gentle, naturalistic rolling hills on the back-left, back-center, and back-right!
  # Water surface level = 2.60 m; deepest lakebed = 0.25 m (so max water depth = 2.35 m).
  level <- 2.60

  # Low-frequency smooth FBM (hurst = 2.8 => very smooth rolling hills, zero jagged needles!)
  fb_smooth <- fbm(n, 2.8)
  fb_tex    <- fbm(n, 2.2)

  # Keep the lake basin deep across the center & foreground (Z > 1.5) so the entire foreground water column
  # has full 2+ m optical depth (vivid sapphire blue in Clear State and authentic cutaway depth at Z = +10!)
  dz_north <- pmin(0, Z - 1.8) / 6.8
  dz_south <- pmax(0, Z - 1.8) / 15.5
  r_basin <- sqrt((X / 8.0)^2 + dz_north^2 + dz_south^2) + 0.06 * (fb_smooth - 0.5)
  basin_profile <- 0.20 + 2.40 * smoothstep_sr(0.0, 1.02, r_basin)^1.65

  # Three gentle rolling hills around the northern/lateral rim (z < 3) with soft natural valley contours
  hill_left   <- 3.2 * exp(-((X + 6.4)^2 / 19.0 + (Z + 4.2)^2 / 18.0))
  hill_center <- 3.8 * exp(-((X - 0.5)^2 / 24.0 + (Z + 6.5)^2 / 15.0))
  hill_right  <- 2.9 * exp(-((X - 6.8)^2 / 18.0 + (Z + 3.5)^2 / 19.0))
  gentle_relief <- (0.85 + 0.30 * fb_smooth) * (hill_left + hill_center + hill_right)

  H <- basin_profile + gentle_relief * smoothstep_sr(0.62, 1.28, r_basin) +
       0.10 * sin(0.65 * X) * cos(0.55 * Z)

  # Build heightfield mesh using sonicray::sr_mesh_from_heightfield (vexag = 1000 => 1 scene unit = 1 metre)
  mesh <- sr_mesh_from_heightfield(H, c(-10, 10), c(-10, 10), y0 = 0, vexag = 1000)

  # Naturalistic shoreline, forest, and sunlit meadow hill coloration:
  H_tex <- rmat(H, tex, tex)
  fb_t  <- rmat(fb_tex, tex, tex)
  gy <- H_tex[c(2:tex, tex), ] - H_tex[c(1, 1:(tex - 1)), ]
  gx <- H_tex[, c(2:tex, tex)] - H_tex[, c(1, 1:(tex - 1))]
  slope_t <- pmin(1, sqrt(gx^2 + gy^2) * 18)

  w_beach <- smoothstep_sr(level - 0.12, level + 0.18, H_tex) * (1 - smoothstep_sr(level + 0.18, level + 0.55, H_tex))
  w_land  <- smoothstep_sr(level + 0.05, level + 0.45, H_tex)
  w_high  <- smoothstep_sr(level + 0.60, level + 3.60, H_tex)

  alb_srgb <- array(0, c(tex, tex, 3))
  bed_col   <- c(0.42, 0.40, 0.34) # Natural mineral lakebed so Clear Water is a deep, vivid sapphire blue!
  beach_col <- c(0.66, 0.58, 0.44) # Sandy shoreline bank
  tree_low  <- c(0.16, 0.31, 0.14) # Riparian forest green
  hill_top  <- c(0.36, 0.41, 0.24) # Warm grassy/rocky hill crest
  rock_col  <- c(0.45, 0.41, 0.35)

  for (k in 1:3) {
    land_k <- (1 - w_high) * tree_low[k] + w_high * hill_top[k]
    land_k <- land_k * (1 - 0.45 * slope_t) + rock_col[k] * (0.45 * slope_t)
    land_k <- land_k * (0.86 + 0.28 * fb_t)
    base_k <- (1 - w_land) * bed_col[k] * (0.90 + 0.20 * fb_t) + w_land * land_k
    alb_srgb[, , k] <- (1 - w_beach) * base_k + w_beach * beach_col[k]
  }
  alb <- pmin(pmax(alb_srgb, 0.02), 0.96)^2.2

  # Paint dedicated UV swatches using `sonicray::sr_paint_albedo()` for:
  #   1. White Secchi marker spheres (`uv_secchi`)
  #   2. 3D Fish dorsal olive (`uv_back`), vibrant coral-orange lateral stripe (`uv_stripe`), and bright silver belly (`uv_belly`)
  uv_secchi <- c(0.006, 0.994)
  uv_back   <- c(0.018, 0.994)
  uv_stripe <- c(0.030, 0.994)
  uv_belly  <- c(0.042, 0.994)

  alb <- sr_paint_albedo(alb, uv_secchi, c(0.96, 0.96, 0.96)^2.2, half_px = 3L)
  alb <- sr_paint_albedo(alb, uv_back,   c(0.25, 0.40, 0.30)^2.2, half_px = 3L)
  alb <- sr_paint_albedo(alb, uv_stripe, c(0.98, 0.38, 0.22)^2.2, half_px = 3L)
  alb <- sr_paint_albedo(alb, uv_belly,  c(0.95, 0.92, 0.88)^2.2, half_px = 3L)

  # Append 4 depth tiers of submerged objects (at depths 0.40 m, 0.90 m, 1.45 m, and 2.05 m below `level = 2.60 m`):
  # Each tier has a 3D swimming trout AND a white Secchi sphere right beside it!
  tiers <- list(
    list(depth = 0.40, fish_c = c(-3.3, level - 0.40, 0.2), sphere_c = c(-4.6, level - 0.40, 0.5), yaw =  0.20, len = 1.65),
    list(depth = 0.90, fish_c = c(-0.7, level - 0.90, 1.4), sphere_c = c(-2.0, level - 0.90, 1.7), yaw =  0.14, len = 1.65),
    list(depth = 1.45, fish_c = c( 1.9, level - 1.45, 2.4), sphere_c = c( 0.6, level - 1.45, 2.7), yaw = -0.16, len = 1.68),
    list(depth = 2.05, fish_c = c( 4.3, level - 2.05, 3.4), sphere_c = c( 3.0, level - 2.05, 3.7), yaw = -0.22, len = 1.72)
  )

  for (tr in tiers) {
    mesh <- sr_mesh_append(mesh, sr_sphere_mesh(tr$sphere_c, radius = 0.28, uv = uv_secchi))
    mesh <- sr_mesh_append(mesh, build_sr_fish_mesh(tr$fish_c, length = tr$len, yaw = tr$yaw,
                                                    uv_back = uv_back, uv_stripe = uv_stripe, uv_belly = uv_belly))
  }

  # Elevated camera and warm side-angled sun (avoids blinding white sun-ball glare inside the water basin)
  cam <- list(from = c(0.0, 10.2, 14.5), at = c(0.0, 1.1, 0.0), fov = 44.0)
  el <- 48 * pi / 180
  az <- 215 * pi / 180
  sun <- data.frame(
    name = "sun",
    x = 300 * cos(el) * cos(az),
    y = 300 * sin(el),
    z = 300 * cos(el) * sin(az),
    r = 12,
    intensity = 2300
  )

  list(mesh = mesh, alb = alb, level = level, tiers = tiers, cam = cam, sun = sun)
}

#' Pre-render the optical water states natively on the Metal GPU via `sr_upload_mesh(water = ...)`
#' and `sr_render_light_layers()`, using `sr_water_preset()` as the starting point.
render_native_water_states <- function(scene, width = 1040L, height = 900L,
                                       amb_spp = 256L, light_spp = 1536L,
                                       cache_dir = "data/derived/sonicray_lake_layers") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_rds <- file.path(cache_dir, sprintf("native_sr_water_states_%dx%d_v3.rds", width, height))
  if (file.exists(cache_rds)) {
    message("[sonicray_lake] Loading cached native {sonicray} water renders from ", cache_rds)
    return(readRDS(cache_rds))
  }

  p_clear  <- sr_water_preset("clear")
  p_algal  <- sr_water_preset("algal")
  p_brown  <- sr_water_preset("brown")
  p_turbid <- sr_water_preset("turbid")

  water_configs <- list(
    # Ultra-clear blue water: strong red/green absorption + blue Rayleigh scatter -> vivid sapphire blue, all fish crisp!
    clear     = list(absorption = c(0.68, 0.140, 0.012), scattering = c(0.002, 0.015, 0.055)),
    # Clear-regime constituent pulses (low OAC load -> water shifts hue gently while keeping fish visible!)
    clear_chl = list(absorption = c(0.65, 0.110, 0.180), scattering = c(0.018, 0.065, 0.035)),
    clear_doc = list(absorption = c(0.55, 0.220, 0.480), scattering = c(0.018, 0.014, 0.006)),
    clear_tss = list(absorption = c(0.48, 0.180, 0.280), scattering = c(0.130, 0.100, 0.045)),
    # Dark/turbid-regime constituent blooms (high OAC load -> water becomes opaque, hiding deep fish & spheres!)
    algal     = list(absorption = p_algal$absorption * c(1.30, 0.88, 1.50), scattering = c(0.05, 0.40, 0.10)),
    brown     = list(absorption = c(0.42, 1.10, 2.90),                      scattering = c(0.15, 0.06, 0.01)),
    turbid    = list(absorption = c(0.38, 0.72, 1.65),                      scattering = c(0.95, 0.56, 0.14))
  )

  render_one <- function(nm, cfg) {
    message(sprintf("[sonicray_lake] Native Metal GPU water render for state '%s'...", nm))
    # Smooth capillary roughness (0.009) avoids grainy normal noise while keeping realistic Fresnel reflection
    w <- c(cfg, list(level = scene$level, cutaway = TRUE, roughness = 0.009, metres_per_unit = 1))
    gpu <- sr_upload_mesh(scene$mesh, scene$alb, water = w)
    r <- sr_render_light_layers(
      gpu, scene$sun,
      lookfrom = scene$cam$from, lookat = scene$cam$at, fov = scene$cam$fov,
      width = as.integer(width), height = as.integer(height),
      amb_spp = as.integer(amb_spp), light_spp = as.integer(light_spp),
      max_bounces = 4L, clamp = 60, batch_size = 1L
    )
    # Balance warm sun layer + sky ambient so brown DOC and tan TSS states don't get washed out by blue sky ambient
    lin <- 6.5 * array(r$aovs$amb, c(3, width, height)) + 1.35 * array(r$layers[[1]], c(3, width, height))
    lin <- aperm(lin, c(3, 2, 1)) # H x W x 3
    mm  <- t(array(r$aovs$mask, c(width, height)))
    wd  <- t(array(r$aovs$water_depth, c(width, height)))
    list(lin = lin, mask = mm, water_depth = wd)
  }

  states_lin <- lapply(names(water_configs), function(nm) render_one(nm, water_configs[[nm]]))
  names(states_lin) <- names(water_configs)

  # Project the 4 fish/sphere depth tiers to 2D pixel coordinates via `sonicray::sr_project()`
  pts3d <- do.call(rbind, lapply(scene$tiers, `[[`, "fish_c"))
  px <- sr_project(pts3d, scene$cam$from, scene$cam$at, scene$cam$fov, width, height)
  tier_df <- tibble(
    depth_m = vapply(scene$tiers, `[[`, 0.0, "depth"),
    x_px = px[, 1],
    y_px = height - px[, 2]
  )

  out <- list(width = width, height = height, states = states_lin, tier_df = tier_df)
  saveRDS(out, cache_rds)
  out
}

#' Convert blended linear HDR radiance (`lin`, H x W x 3) to sRGB using {sonicray}'s ACES tone-mapper (`tools/water_demo.R`)
tonemap_sr_water <- function(lin, mask) {
  H <- dim(lin)[1]; W <- dim(lin)[2]
  x <- pmax(lin, 0)
  aces <- (x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14)
  img <- pmin(pmax(aces, 0), 1)^(1 / 2.2)
  sky <- array(rep(c(0.05, 0.08, 0.16), each = H * W), c(H, W, 3))
  for (k in 1:3) img[, , k] <- sky[, , k] * (1 - mask) + img[, , k] * mask
  img
}

#' Bio-optical state calculator
bio_optical_state <- function(chl_a, doc, nap,
                              kw = 0.0478, a_nap = 0.0432, a_doc = 0.0732, a_chl = 0.0086) {
  k_nap <- a_nap * nap
  k_doc <- a_doc * doc
  k_chl <- a_chl * chl_a
  k_p   <- pmax(1e-6, k_nap + k_doc + k_chl)
  kd    <- kw + k_p
  z_sd  <- 1 / kd

  f_nap <- k_nap / k_p
  f_doc <- k_doc / k_p
  f_chl <- k_chl / k_p

  # Smooth transition parameter `darkness` in [0, 1] based on total non-water attenuation k_p:
  #   - When k_p < 0.18 (z_SD > 4.4 m): darkness ~ 0 (Clear Water Regime)
  #   - When k_p > 1.20 (z_SD < 0.8 m): darkness ~ 1 (Dark / Turbid Water Regime)
  darkness <- 1 - exp(-0.85 * pmax(0, k_p - 0.08))
  w_water  <- exp(-1.35 * k_p)

  c_water <- c(0.10, 0.54, 0.94)
  c_chl   <- (1 - darkness) * c(0.14, 0.78, 0.62) + darkness * c(0.10, 0.48, 0.22)
  c_doc   <- (1 - darkness) * c(0.42, 0.46, 0.34) + darkness * c(0.24, 0.14, 0.06)
  c_nap   <- (1 - darkness) * c(0.56, 0.62, 0.52) + darkness * c(0.72, 0.56, 0.32)

  col_srgb <- w_water * c_water + (1 - w_water) * (f_chl * c_chl + f_doc * c_doc + f_nap * c_nap)
  col_srgb <- pmin(pmax(col_srgb, 0.04), 0.96)
  hex_col  <- rgb(col_srgb[1], col_srgb[2], col_srgb[3])

  list(
    chl_a = chl_a, doc = doc, nap = nap,
    k_nap = k_nap, k_doc = k_doc, k_chl = k_chl, kw = kw, kd = kd, z_sd = z_sd,
    f_nap = f_nap, f_doc = f_doc, f_chl = f_chl,
    col_srgb = col_srgb, hex = hex_col,
    darkness = darkness
  )
}

#' 28-second movie loop that progresses IN TIME from:
#'   - Act I (0 - 12 s):  CLEAR WATER STATE (low OACs, z_SD = 4 - 11 m; fish & spheres visible across the ternary!)
#'   - Act II (12 - 28 s): DARK / TURBID WATER STATE (high OACs, z_SD = 0.3 - 0.8 m; titrating Chl a, DOC, and TSS)
titration_schedule <- function(t_sec, total_s = 28.0) {
  chl0 <- 0.30; doc0 <- 0.35; nap0 <- 0.35

  if (t_sec < 4.5) {
    # Stage 1: Ultra-clear oligotrophic baseline
    phase <- "Act I (Clear State): Oligotrophic Lake"
    chl <- chl0; doc <- doc0; nap <- nap0
    w_states <- c(clear = 1.0, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = 0, brown = 0, turbid = 0)
  } else if (t_sec < 12.0) {
    # Stage 2: Sweeping around the Ternary Diagram in the CLEAR WATER regime (low total OAC concentrations!)
    u <- (t_sec - 4.5) / 7.5
    a1 <- exp(-((u - 0.22) / 0.16)^2)
    a2 <- exp(-((u - 0.52) / 0.16)^2)
    a3 <- exp(-((u - 0.82) / 0.16)^2)
    chl <- chl0 + 14.0 * a1 + 1.5 * a2 + 1.5 * a3
    doc <- doc0 + 0.3 * a1  + 2.2 * a2 + 0.3 * a3
    nap <- nap0 + 0.4 * a1  + 0.4 * a2 + 3.8 * a3
    phase <- "Act I (Clear State): Low-Load Ternary Sweep"
    w_c <- pmax(0, 1 - 0.75 * (a1 + a2 + a3))
    w_states <- c(clear = w_c, clear_chl = 0.75 * a1, clear_doc = 0.75 * a2, clear_tss = 0.75 * a3,
                  algal = 0, brown = 0, turbid = 0)
    w_states <- w_states / sum(w_states)
  } else if (t_sec < 17.5) {
    # Stage 3: Act II (Dark/Turbid Regime) — Titrating High Chlorophyll a (Algal Bloom)
    u <- (t_sec - 12.0) / 5.5
    p <- sin(pi * smooth01(u))^1.15
    chl <- chl0 + 135.0 * p
    doc <- doc0 + 0.4 * p
    nap <- nap0 + 0.5 * p
    phase <- "Act II (Dark State): Algal Bloom (High Chl a)"
    w_states <- c(clear = 1 - p, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = p, brown = 0, turbid = 0)
  } else if (t_sec < 23.0) {
    # Stage 4: Act II (Dark/Turbid Regime) — Titrating High DOC (Brownification)
    u <- (t_sec - 17.5) / 5.5
    p <- sin(pi * smooth01(u))^1.15
    chl <- chl0 + 0.6 * p
    doc <- doc0 + 28.0 * p
    nap <- nap0 + 0.4 * p
    phase <- "Act II (Dark State): Brownification (High DOC)"
    w_states <- c(clear = 1 - p, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = 0, brown = p, turbid = 0)
  } else {
    # Stage 5: Act II (Dark/Turbid Regime) — Titrating High TSS / NAP (Sediment Plume)
    u <- (t_sec - 23.0) / 5.0
    p <- sin(pi * smooth01(u))^1.15
    chl <- chl0 + 0.6 * p
    doc <- doc0 + 0.4 * p
    nap <- nap0 + 75.0 * p
    phase <- "Act II (Dark State): Sediment Plume (High TSS)"
    w_states <- c(clear = 1 - p, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = 0, brown = 0, turbid = p)
  }

  st <- bio_optical_state(chl, doc, nap)
  st$t_sec <- t_sec
  st$phase <- phase
  st$w_states <- w_states
  st
}

tern_xy_pt <- function(f_nap, f_doc, f_chl) {
  s <- f_nap + f_doc + f_chl
  f_nap <- f_nap / s; f_doc <- f_doc / s; f_chl <- f_chl / s
  x <- f_doc + 0.5 * f_chl
  y <- (sqrt(3) / 2) * f_chl
  cbind(x = x, y = y)
}

build_ternary_grid <- function(res = 28) {
  pts <- list()
  for (i in 0:res) {
    for (j in 0:(res - i)) {
      k <- res - i - j
      f_nap <- max(0.01, i / res)
      f_doc <- max(0.01, j / res)
      f_chl <- max(0.01, k / res)
      s <- f_nap + f_doc + f_chl
      f_nap <- f_nap / s; f_doc <- f_doc / s; f_chl <- f_chl / s

      # Clear-state color (low k_p = 0.08 m^-1) vs Dark-state color (high k_p = 2.6 m^-1)
      kp_clear <- 0.08
      st_c <- bio_optical_state((f_chl * kp_clear) / 0.0086, (f_doc * kp_clear) / 0.0732, (f_nap * kp_clear) / 0.0432)
      kp_dark  <- 2.60
      st_d <- bio_optical_state((f_chl * kp_dark) / 0.0086,  (f_doc * kp_dark) / 0.0732,  (f_nap * kp_dark) / 0.0432)

      xy <- tern_xy_pt(f_nap, f_doc, f_chl)
      pts[[length(pts) + 1]] <- tibble(
        f_nap = f_nap, f_doc = f_doc, f_chl = f_chl,
        x = xy[1, "x"], y = xy[1, "y"],
        r_c = st_c$col_srgb[1], g_c = st_c$col_srgb[2], b_c = st_c$col_srgb[3],
        r_d = st_d$col_srgb[1], g_d = st_d$col_srgb[2], b_d = st_d$col_srgb[3]
      )
    }
  }
  bind_rows(pts)
}

#' Render a single widescreen composite frame (1600 x 900) with ONE large dynamic Ternary Diagram
draw_composite_frame <- function(out_png, lake_img, gpu_renders, st, trail_df, tgrid, width = 1600, height = 900) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(out_png, width = width, height = height, units = "px", res = 144, background = "#0b111e")
  } else {
    png(out_png, width = width, height = height, units = "px", res = 144, bg = "#0b111e", type = "quartz")
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bg = "#0b111e", fg = "white")
  plot.new()
  plot.window(xlim = c(0, 1600), ylim = c(0, 900))

  # 1. Left Panel (x: 0 -> 1040, y: 0 -> 900): Native {sonicray} 3D Water Basin
  rasterImage(lake_img, xleft = 0, ybottom = 0, xright = 1040, ytop = 900, interpolate = TRUE)

  # Subtle depth labels above the 4 submerged fish/Secchi-sphere pairs (fading out as kd hides them)
  tdf <- gpu_renders$tier_df
  for (i in seq_len(nrow(tdf))) {
    vis <- exp(-0.95 * st$kd * tdf$depth_m[i])
    if (vis > 0.22) {
      text(tdf$x_px[i], tdf$y_px[i] + 22, sprintf("%.1f m", tdf$depth_m[i]),
           col = rgb(1, 1, 1, pmin(1, vis * 1.15)), font = 2, cex = 0.60)
    }
  }

  # Clean top title bar on the 3D lake panel
  rect(24, 830, 1016, 880, col = rgb(0.05, 0.08, 0.14, 0.78), border = rgb(1, 1, 1, 0.18), lwd = 1.2)
  text(42, 855, st$phase, adj = c(0, 0.5), col = "white", font = 2, cex = 0.98)
  rect(818, 836, 1004, 874, col = st$hex, border = "white", lwd = 1.6)
  tcol_top <- if (mean(st$col_srgb) > 0.52) "#0f172a" else "#ffffff"
  text(911, 855, sprintf("Secchi = %.1f m", st$z_sd), col = tcol_top, font = 2, cex = 0.84)

  # ===========================================================================
  # 2. Right Panel (x: 1040 -> 1600, y: 0 -> 900): Single Dynamic Ternary + OAC Bars
  # ===========================================================================
  rect(1040, 0, 1600, 900, col = "#0b111e", border = "#334155", lwd = 1.5)

  # Top: 3 Constituent Bars (y: 665 -> 875)
  text(1068, 860, "Constituent Titration", adj = c(0, 0.5), col = "white", font = 2, cex = 0.92)

  draw_clean_bar <- function(yc, name, val_txt, pct_txt, frac, col_bar) {
    text(1068, yc + 11, name, adj = c(0, 0.5), col = "#e2e8f0", font = 2, cex = 0.72)
    text(1572, yc + 11, paste0(val_txt, "  (", pct_txt, ")"), adj = c(1, 0.5), col = "#f8fafc", font = 2, cex = 0.70)
    rect(1068, yc - 11, 1572, yc - 1, col = "#1e293b", border = "#475569", lwd = 1)
    w <- pmin(1, pmax(0.015, frac)) * (1572 - 1068)
    rect(1068, yc - 11, 1068 + w, yc - 1, col = col_bar, border = NA)
  }
  draw_clean_bar(810, "Chl a (Algae)", sprintf("%.1f ug/L", st$chl_a), sprintf("%.0f%%", 100 * st$f_chl), (st$chl_a / 135)^0.60, "#10b981")
  draw_clean_bar(752, "DOC (Brown)",   sprintf("%.1f mg/L", st$doc),   sprintf("%.0f%%", 100 * st$f_doc), (st$doc   / 28)^0.60,  "#b45309")
  draw_clean_bar(694, "TSS (Turbid)",  sprintf("%.1f mg/L", st$nap),   sprintf("%.0f%%", 100 * st$f_nap), (st$nap   / 75)^0.60,  "#eab308")

  # 100% Partition Bar
  x0 <- 1068; xw <- 504
  w1 <- st$f_chl * xw; w2 <- st$f_doc * xw
  rect(x0, 638, x0 + w1, 658, col = "#10b981", border = NA)
  rect(x0 + w1, 638, x0 + w1 + w2, 658, col = "#78350f", border = NA)
  rect(x0 + w1 + w2, 638, x0 + xw, 658, col = "#eab308", border = NA)
  rect(x0, 638, x0 + xw, 658, col = NA, border = "#cbd5e1", lwd = 1.1)

  # Center: ONE Large, Dynamic Ternary Diagram that transitions in time from Clear State -> Dark State!
  segments(1064, 614, 1576, 614, col = "#334155", lwd = 1.2)
  regime_label <- if (st$darkness < 0.35) "CLEAR WATER REGIME (High Visibility)" else "DARK / TURBID REGIME (Low Visibility)"
  regime_col   <- if (st$darkness < 0.35) "#38bdf8" else "#fde047"
  text(1320, 586, regime_label, col = regime_col, font = 2, cex = 0.86)

  rect(1076, 198, 1564, 560, col = "#131c2e", border = regime_col, lwd = 2.0)

  cx_left <- 1152; cx_right <- 1488
  cy_bot  <- 236;  tw <- cx_right - cx_left
  th <- tw * (sqrt(3) / 2)
  tx <- function(x) cx_left + x * tw
  ty <- function(y) cy_bot + (y / (sqrt(3) / 2)) * th

  # Blend ternary grid colors between Clear State (`r_c, g_c, b_c`) and Dark State (`r_d, g_d, b_d`) by `st$darkness`!
  dk <- st$darkness
  r_mix <- (1 - dk) * tgrid$r_c + dk * tgrid$r_d
  g_mix <- (1 - dk) * tgrid$g_c + dk * tgrid$g_d
  b_mix <- (1 - dk) * tgrid$b_c + dk * tgrid$b_d
  cols_dyn <- rgb(r_mix, g_mix, b_mix)

  points(tx(tgrid$x), ty(tgrid$y), pch = 16, cex = 1.25, col = cols_dyn)

  s32 <- sqrt(3) / 2
  polygon(c(tx(0), tx(1), tx(0.5)), c(ty(0), ty(0), ty(s32)), col = NA, border = "white", lwd = 2.0)

  text(tx(0.5), ty(s32) + 15, "100% Chl a (Algal)", col = "#4ade80", font = 2, cex = 0.74)
  text(tx(0.0) - 6, ty(0) - 16, "100% TSS", adj = c(0, 0.5), col = "#fde047", font = 2, cex = 0.74)
  text(tx(1.0) + 6, ty(0) - 16, "100% DOC", adj = c(1, 0.5), col = "#f59e0b", font = 2, cex = 0.74)

  if (nrow(trail_df) > 1) {
    lines(tx(trail_df$x), ty(trail_df$y), col = rgb(1, 1, 1, 0.65), lwd = 2.0)
  }
  cur_xy <- tern_xy_pt(st$f_nap, st$f_doc, st$f_chl)
  points(tx(cur_xy[1, "x"]), ty(cur_xy[1, "y"]), pch = 21, bg = st$hex, col = "#fde047", cex = 2.6, lwd = 2.4)

  # Bottom: Visual Underwater Visibility Meter (0.4 m, 0.9 m, 1.5 m, 2.1 m)
  rect(1076, 26, 1564, 178, col = "#131c2e", border = "#334155", lwd = 1.2)
  text(1320, 156, "Submerged Fish & Secchi Sphere Visibility", col = "white", font = 2, cex = 0.76)

  f_depths <- c(0.4, 0.9, 1.5, 2.1)
  x_centers <- seq(1144, 1496, length.out = 4)
  for (k in seq_along(f_depths)) {
    zm <- f_depths[k]
    xc <- x_centers[k]
    vis <- exp(-0.95 * st$kd * zm)
    box_col <- if (vis > 0.25) "#1e293b" else "#0f172a"
    brd_col <- if (vis > 0.25) "#38bdf8" else "#334155"
    rect(xc - 48, 42, xc + 48, 134, col = box_col, border = brd_col, lwd = if (vis > 0.25) 1.6 else 1.0)
    t_seq <- seq(0, 2 * pi, length.out = 32)
    bx <- xc + 20 * cos(t_seq)
    by <- 98 + 6.5 * sin(t_seq) * (0.75 + 0.35 * cos(t_seq))
    tx_p <- xc - 20 * c(0.75, 1.25, 1.05, 1.25)
    ty_p <- 98 + 6.5 * c(0, 0.9, 0, -0.9)
    f_alpha <- pmax(0.08, pmin(1.0, vis))
    polygon(tx_p, ty_p, col = rgb(0.95, 0.45, 0.35, f_alpha), border = NA)
    polygon(bx, by, col = rgb(0.92, 0.82, 0.72, f_alpha), border = rgb(1, 1, 1, f_alpha * 0.6), lwd = 0.8)
    text(xc, 70, sprintf("%.1f m", zm), col = "white", font = 2, cex = 0.68)
    status_txt <- if (vis > 0.50) "Clear" else if (vis > 0.22) "Faint" else "Hidden"
    status_col <- if (vis > 0.50) "#4ade80" else if (vis > 0.22) "#fde047" else "#64748b"
    text(xc, 52, status_txt, col = status_col, font = 2, cex = 0.60)
  }
}

#' Build the 6-keyframe contact sheet and full 28-second MP4 using native {sonicray} water rendering
build_sonicray_titration_visuals <- function(
  out_mp4   = "figs/sonicray_optical_titration.mp4",
  out_sheet = "figs/sonicray_optical_titration_sheet.png",
  fps       = 15L,
  duration_s = 28.0,
  sheet_only = FALSE
) {
  t0 <- Sys.time()
  scene <- build_sonicray_water_scene(n = 200L, tex = 1024L)
  gpu_renders <- render_native_water_states(scene, width = 1040L, height = 900L, amb_spp = 256L, light_spp = 1536L)
  tgrid <- build_ternary_grid(res = 28)

  times <- seq(0, duration_s, by = 1 / fps)
  states <- lapply(times, titration_schedule, total_s = duration_s)
  xy_all <- do.call(rbind, lapply(states, function(s) tern_xy_pt(s$f_nap, s$f_doc, s$f_chl)))
  trail_all <- tibble(x = xy_all[, "x"], y = xy_all[, "y"])

  # 6 keyframes:
  #   1. Clear Oligotrophic baseline (t = 2.0 s)
  #   2. Act I Clear State Low-Load Sweep — Chl a / DOC / TSS in clear water (t = 6.2 s)
  #   3. Act I Clear State Low-Load Sweep — TSS corner in clear water (t = 10.6 s)
  #   4. Act II Dark State — Algal Bloom (t = 14.8 s)
  #   5. Act II Dark State — DOC Brownification (t = 20.2 s)
  #   6. Act II Dark State — TSS Sediment Plume (t = 25.5 s)
  key_times <- c(2.0, 6.2, 10.6, 14.8, 20.2, 25.5)
  key_idx <- vapply(key_times, function(kt) which.min(abs(times - kt)), 1L)

  frames_dir <- file.path(tempdir(), "sonicray_titration_frames_native")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

  render_indices <- if (sheet_only) key_idx else seq_along(times)
  message(sprintf("[sonicray_lake] Compositing %d frames from native {sonicray} water layers...", length(render_indices)))

  state_names <- names(gpu_renders$states)
  mask_mat <- gpu_renders$states[[1]]$mask

  render_one_idx <- function(i) {
    f_png <- file.path(frames_dir, sprintf("frame_%04d.png", i))
    st_i <- states[[i]]
    w_vec <- st_i$w_states
    lin_blend <- array(0.0, dim(gpu_renders$states[[1]]$lin))
    for (nm in state_names) {
      wt <- w_vec[[nm]]
      if (wt > 1e-4) {
        lin_blend <- lin_blend + wt * gpu_renders$states[[nm]]$lin
      }
    }
    img3d <- tonemap_sr_water(lin_blend, mask_mat)
    draw_composite_frame(f_png, img3d, gpu_renders, st_i, trail_all[1:i, , drop = FALSE], tgrid)
    f_png
  }

  if (sheet_only) {
    lapply(render_indices, render_one_idx)
  } else {
    parallel::mclapply(render_indices, render_one_idx, mc.cores = min(8L, parallel::detectCores()))
  }

  # Assemble 3x2 contact sheet of the 6 key stages
  key_imgs <- lapply(key_idx, function(idx) png::readPNG(file.path(frames_dir, sprintf("frame_%04d.png", idx))))
  fh <- dim(key_imgs[[1]])[1]; fw <- dim(key_imgs[[1]])[2]
  sh <- round(fh * 0.5); sw <- round(fw * 0.5)
  pad <- 10L
  sheet <- array(0.05, c(3 * sh + 4 * pad, 2 * sw + 3 * pad, 3))
  yi <- round(seq(1, fh, length.out = sh))
  xi <- round(seq(1, fw, length.out = sw))
  for (k in seq_along(key_imgs)) {
    r0 <- ((k - 1) %/% 2) * (sh + pad) + pad
    c0 <- ((k - 1) %% 2)  * (sw + pad) + pad
    sheet[r0 + seq_len(sh), c0 + seq_len(sw), ] <- key_imgs[[k]][yi, xi, 1:3]
  }
  png::writePNG(sheet, out_sheet)
  message("[sonicray_lake] Wrote keyframe sheet: ", out_sheet)

  if (!sheet_only) {
    ff <- Sys.which("ffmpeg")
    cmd <- sprintf(
      '"%s" -y -framerate %d -i "%s/frame_%%04d.png" -c:v libx264 -preset fast -crf 18 -pix_fmt yuv420p -movflags +faststart "%s"',
      ff, fps, frames_dir, out_mp4
    )
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    message(sprintf("[sonicray_lake] Wrote %s (%.2f MB) in %.1f s total.",
                    out_mp4, file.info(out_mp4)$size / 1e6, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  invisible(list(mp4 = out_mp4, sheet = out_sheet))
}
