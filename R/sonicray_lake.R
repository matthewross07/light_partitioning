# -----------------------------------------------------------------------------
# sonicray_lake.R -- High-Resolution Native {sonicray} 3D Water Simulation
#                    with Smooth 3-Hill Lake Basin, Submerged 3D Fish & Secchi
#                    Spheres, Edge-Preserving GPU AOV Denoising, and Exact
#                    Manuscript-Style Ternary Diagram (matching Fig. 2 & Fig. 3).
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(sonicray)
  library(tidyverse)
  library(png)
})

S3_TERN <- sqrt(3) / 2

smooth01 <- function(x) {
  x <- pmin(pmax(x, 0), 1)
  x * x * (3 - 2 * x)
}

smoothstep_sr <- function(e0, e1, x) {
  smooth01((x - e0) / (e1 - e0))
}

#' Fast edge-preserving joint bilateral filter guided by {sonicray}'s deterministic
#' first-hit AOVs (`albedo` and `water_depth`) so Monte Carlo water variance is silky
#' smooth while 3D fish, Secchi spheres, shadows, and shoreline edges stay razor-sharp.
denoise_sr_water_hdr <- function(lin, albedo_aov, wdepth_aov, mask_aov, radius = 3L) {
  H <- dim(lin)[1]; W <- dim(lin)[2]
  out <- lin
  w_sum <- matrix(1.0, H, W)

  alb_r <- albedo_aov[, , 1]
  alb_g <- albedo_aov[, , 2]
  alb_b <- albedo_aov[, , 3]
  is_wat <- (wdepth_aov > 1e-4)

  for (dy in (-radius):radius) {
    for (dx in (-radius):radius) {
      if (dx == 0L && dy == 0L) next
      # Shifted index ranges
      r_src <- pmin(pmax(seq_len(H) + dy, 1L), H)
      c_src <- pmin(pmax(seq_len(W) + dx, 1L), W)

      # Guide distance in deterministic AOV space (albedo + water_depth + mask)
      d_alb <- (alb_r - alb_r[r_src, c_src])^2 +
               (alb_g - alb_g[r_src, c_src])^2 +
               (alb_b - alb_b[r_src, c_src])^2
      d_dep <- (wdepth_aov - wdepth_aov[r_src, c_src])^2
      d_msk <- (mask_aov - mask_aov[r_src, c_src])^2
      d_spa <- dx * dx + dy * dy

      # Only smooth across pixels sharing the same underwater surface/object (sharp edges preserved!)
      w <- is_wat * is_wat[r_src, c_src] *
           exp(-d_spa / (2 * 2.2^2) - d_alb / (2 * 0.015^2) - d_dep / (2 * 0.14^2) - d_msk / (2 * 0.05^2))

      for (k in 1:3) {
        out[, , k] <- out[, , k] + w * lin[r_src, c_src, k]
      }
      w_sum <- w_sum + w
    }
  }
  for (k in 1:3) out[, , k] <- out[, , k] / w_sum
  out
}

#' Build a 3D streamlined trout/fish mesh to append via `sonicray::sr_mesh_append()`
build_sr_fish_mesh <- function(center, length = 1.65, yaw = 0.0,
                               uv_back = c(0.018, 0.994),
                               uv_stripe = c(0.030, 0.994),
                               uv_belly = c(0.042, 0.994)) {
  n_len <- 24L
  n_ang <- 18L
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
build_sonicray_water_scene <- function(n = 240L, tex = 2048L) {
  set.seed(42L)
  fbm  <- getFromNamespace("fbm_noise", "sonicray")
  rmat <- getFromNamespace("resize_mat", "sonicray")

  xs <- seq(-10, 10, length.out = n)
  zs <- seq(-10, 10, length.out = n)
  X <- matrix(rep(xs, each = n), n, n)
  Z <- matrix(rep(zs, times = n), n, n)

  level <- 2.60

  fb_smooth <- fbm(n, 2.8)
  fb_tex    <- fbm(n, 2.2)

  dz_north <- pmin(0, Z - 1.8) / 6.8
  dz_south <- pmax(0, Z - 1.8) / 15.5
  r_basin <- sqrt((X / 8.0)^2 + dz_north^2 + dz_south^2) + 0.06 * (fb_smooth - 0.5)
  basin_profile <- 0.20 + 2.40 * smoothstep_sr(0.0, 1.02, r_basin)^1.65

  hill_left   <- 3.2 * exp(-((X + 6.4)^2 / 19.0 + (Z + 4.2)^2 / 18.0))
  hill_center <- 3.8 * exp(-((X - 0.5)^2 / 24.0 + (Z + 6.5)^2 / 15.0))
  hill_right  <- 2.9 * exp(-((X - 6.8)^2 / 18.0 + (Z + 3.5)^2 / 19.0))
  gentle_relief <- (0.85 + 0.30 * fb_smooth) * (hill_left + hill_center + hill_right)

  H <- basin_profile + gentle_relief * smoothstep_sr(0.62, 1.28, r_basin) +
       0.10 * sin(0.65 * X) * cos(0.55 * Z)

  mesh <- sr_mesh_from_heightfield(H, c(-10, 10), c(-10, 10), y0 = 0, vexag = 1000)

  H_tex <- rmat(H, tex, tex)
  fb_t  <- rmat(fb_tex, tex, tex)
  gy <- H_tex[c(2:tex, tex), ] - H_tex[c(1, 1:(tex - 1)), ]
  gx <- H_tex[, c(2:tex, tex)] - H_tex[, c(1, 1:(tex - 1))]
  slope_t <- pmin(1, sqrt(gx^2 + gy^2) * 18)

  w_beach <- smoothstep_sr(level - 0.12, level + 0.18, H_tex) * (1 - smoothstep_sr(level + 0.18, level + 0.55, H_tex))
  w_land  <- smoothstep_sr(level + 0.05, level + 0.45, H_tex)
  w_high  <- smoothstep_sr(level + 0.60, level + 3.60, H_tex)

  alb_srgb <- array(0, c(tex, tex, 3))
  bed_col   <- c(0.42, 0.40, 0.34)
  beach_col <- c(0.66, 0.58, 0.44)
  tree_low  <- c(0.16, 0.31, 0.14)
  hill_top  <- c(0.36, 0.41, 0.24)
  rock_col  <- c(0.45, 0.41, 0.35)

  for (k in 1:3) {
    land_k <- (1 - w_high) * tree_low[k] + w_high * hill_top[k]
    land_k <- land_k * (1 - 0.45 * slope_t) + rock_col[k] * (0.45 * slope_t)
    land_k <- land_k * (0.86 + 0.28 * fb_t)
    base_k <- (1 - w_land) * bed_col[k] * (0.92 + 0.16 * fb_t) + w_land * land_k
    alb_srgb[, , k] <- (1 - w_beach) * base_k + w_beach * beach_col[k]
  }
  alb <- pmin(pmax(alb_srgb, 0.02), 0.96)^2.2

  uv_secchi <- c(0.006, 0.994)
  uv_back   <- c(0.018, 0.994)
  uv_stripe <- c(0.030, 0.994)
  uv_belly  <- c(0.042, 0.994)

  alb <- sr_paint_albedo(alb, uv_secchi, c(0.97, 0.97, 0.97)^2.2, half_px = 5L)
  alb <- sr_paint_albedo(alb, uv_back,   c(0.24, 0.38, 0.28)^2.2, half_px = 5L)
  alb <- sr_paint_albedo(alb, uv_stripe, c(0.98, 0.36, 0.20)^2.2, half_px = 5L)
  alb <- sr_paint_albedo(alb, uv_belly,  c(0.96, 0.93, 0.89)^2.2, half_px = 5L)

  tiers <- list(
    list(depth = 0.40, fish_c = c(-3.3, level - 0.40, 0.2), sphere_c = c(-4.6, level - 0.40, 0.5), yaw =  0.20, len = 1.65),
    list(depth = 0.90, fish_c = c(-0.7, level - 0.90, 1.4), sphere_c = c(-2.0, level - 0.90, 1.7), yaw =  0.14, len = 1.65),
    list(depth = 1.45, fish_c = c( 1.9, level - 1.45, 2.4), sphere_c = c( 0.6, level - 1.45, 2.7), yaw = -0.16, len = 1.68),
    list(depth = 2.05, fish_c = c( 4.3, level - 2.05, 3.4), sphere_c = c( 3.0, level - 2.05, 3.7), yaw = -0.22, len = 1.72)
  )

  for (tr in tiers) {
    mesh <- sr_mesh_append(mesh, sr_sphere_mesh(tr$sphere_c, radius = 0.29, uv = uv_secchi, n_lat = 24L, n_lon = 36L))
    mesh <- sr_mesh_append(mesh, build_sr_fish_mesh(tr$fish_c, length = tr$len, yaw = tr$yaw,
                                                    uv_back = uv_back, uv_stripe = uv_stripe, uv_belly = uv_belly))
  }

  cam <- list(from = c(0.0, 10.2, 14.5), at = c(0.0, 1.1, 0.0), fov = 44.0)
  el <- 48 * pi / 180
  az <- 215 * pi / 180
  sun <- data.frame(
    name = "sun",
    x = 300 * cos(el) * cos(az),
    y = 300 * sin(el),
    z = 300 * cos(el) * sin(az),
    r = 14,
    intensity = 2300
  )

  list(mesh = mesh, alb = alb, level = level, tiers = tiers, cam = cam, sun = sun)
}

#' Pre-render the optical water states natively on the Metal GPU at high resolution & sample count
#' and apply joint bilateral AOV-guided denoising for silky-smooth water and razor-sharp fish/spheres.
render_native_water_states <- function(scene, width = 1140L, height = 1080L,
                                       amb_spp = 768L, light_spp = 4096L,
                                       cache_dir = "data/derived/sonicray_lake_layers") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_rds <- file.path(cache_dir, sprintf("native_sr_water_states_%dx%d_hq_v4.rds", width, height))
  if (file.exists(cache_rds)) {
    message("[sonicray_lake] Loading cached high-quality native {sonicray} water renders from ", cache_rds)
    return(readRDS(cache_rds))
  }

  p_algal <- sr_water_preset("algal")

  water_configs <- list(
    clear     = list(absorption = c(0.68, 0.140, 0.012), scattering = c(0.002, 0.015, 0.055)),
    clear_chl = list(absorption = c(0.65, 0.110, 0.180), scattering = c(0.018, 0.065, 0.035)),
    clear_doc = list(absorption = c(0.55, 0.220, 0.480), scattering = c(0.018, 0.014, 0.006)),
    clear_tss = list(absorption = c(0.48, 0.180, 0.280), scattering = c(0.130, 0.100, 0.045)),
    algal     = list(absorption = p_algal$absorption * c(1.30, 0.88, 1.50), scattering = c(0.05, 0.40, 0.10)),
    brown     = list(absorption = c(0.42, 1.10, 2.90),                      scattering = c(0.15, 0.06, 0.01)),
    turbid    = list(absorption = c(0.38, 0.72, 1.65),                      scattering = c(0.95, 0.56, 0.14))
  )

  render_one <- function(nm, cfg) {
    message(sprintf("[sonicray_lake] High-quality Metal GPU water render for state '%s' (%dx%d @ %d spp)...",
                    nm, width, height, light_spp))
    # Mirror-smooth dielectric surface (roughness = 0.002) eliminates normal-tilt speckle
    w <- c(cfg, list(level = scene$level, cutaway = TRUE, roughness = 0.002, metres_per_unit = 1))
    gpu <- sr_upload_mesh(scene$mesh, scene$alb, water = w)
    r <- sr_render_light_layers(
      gpu, scene$sun,
      lookfrom = scene$cam$from, lookat = scene$cam$at, fov = scene$cam$fov,
      width = as.integer(width), height = as.integer(height),
      amb_spp = as.integer(amb_spp), light_spp = as.integer(light_spp),
      max_bounces = 4L, clamp = 60, batch_size = 1L
    )
    lin <- 6.5 * array(r$aovs$amb, c(3, width, height)) + 1.35 * array(r$layers[[1]], c(3, width, height))
    lin <- aperm(lin, c(3, 2, 1)) # H x W x 3
    alb_aov <- aperm(array(r$aovs$albedo, c(3, width, height)), c(3, 2, 1))
    mm  <- t(array(r$aovs$mask, c(width, height)))
    wd  <- t(array(r$aovs$water_depth, c(width, height)))

    # Edge-preserving AOV-guided bilateral filter removes residual Monte Carlo scattering grain
    lin_clean <- denoise_sr_water_hdr(lin, alb_aov, wd, mm, radius = 3L)
    list(lin = lin_clean, mask = mm, water_depth = wd)
  }

  states_lin <- lapply(names(water_configs), function(nm) render_one(nm, water_configs[[nm]]))
  names(states_lin) <- names(water_configs)

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

#' ACES filmic tone-mapper (`tools/water_demo.R`)
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

  darkness <- 1 - exp(-0.85 * pmax(0, k_p - 0.08))
  w_water  <- exp(-1.35 * k_p)

  c_water <- c(0.10, 0.54, 0.94)
  c_chl   <- (1 - darkness) * c(0.14, 0.78, 0.62) + darkness * c(0.10, 0.48, 0.22)
  c_doc   <- (1 - darkness) * c(0.42, 0.46, 0.34) + darkness * c(0.32, 0.18, 0.07)
  c_nap   <- (1 - darkness) * c(0.56, 0.62, 0.52) + darkness * c(0.76, 0.64, 0.38)

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

#' 28-second movie loop: Act I (Clear Water Regime) -> Act II (Dark/Turbid Water Regime)
titration_schedule <- function(t_sec, total_s = 28.0) {
  chl0 <- 0.30; doc0 <- 0.35; nap0 <- 0.35

  if (t_sec < 4.5) {
    phase <- "Act I (Clear State): Oligotrophic Lake"
    chl <- chl0; doc <- doc0; nap <- nap0
    w_states <- c(clear = 1.0, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = 0, brown = 0, turbid = 0)
  } else if (t_sec < 12.0) {
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
    u <- (t_sec - 12.0) / 5.5
    p <- sin(pi * smooth01(u))^1.15
    chl <- chl0 + 135.0 * p
    doc <- doc0 + 0.4 * p
    nap <- nap0 + 0.5 * p
    phase <- "Act II (Dark State): Algal Bloom (High Chl a)"
    w_states <- c(clear = 1 - p, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = p, brown = 0, turbid = 0)
  } else if (t_sec < 23.0) {
    u <- (t_sec - 17.5) / 5.5
    p <- sin(pi * smooth01(u))^1.15
    chl <- chl0 + 0.6 * p
    doc <- doc0 + 28.0 * p
    nap <- nap0 + 0.4 * p
    phase <- "Act II (Dark State): Brownification (High DOC)"
    w_states <- c(clear = 1 - p, clear_chl = 0, clear_doc = 0, clear_tss = 0, algal = 0, brown = p, turbid = 0)
  } else {
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

#' EXACT Ternary Coordinate Mapping from `R/ternary.R` (Figures 2, 3, 4 of the paper):
#'   Top vertex (0.5, sqrt(3)/2) = DOC (100%)
#'   Bottom-left vertex (0, 0)   = Chl a (100%)
#'   Bottom-right vertex (1, 0)  = NAP (100%)
#'   Cartesian: x = f_nap + 0.5 * f_doc,  y = f_doc * sqrt(3)/2
tern_xy_paper <- function(f_chl, f_doc, f_nap) {
  s <- f_chl + f_doc + f_nap
  f_chl <- f_chl / s; f_doc <- f_doc / s; f_nap <- f_nap / s
  cbind(x = f_nap + 0.5 * f_doc, y = f_doc * S3_TERN)
}

#' Build a seamless sub-triangle mesh of the ternary simplex (res = 36 => 1,296 seamless sub-triangles!)
build_paper_ternary_mesh <- function(res = 36L) {
  # Build sub-triangles on the simplex grid (i = nap index, j = doc index, k = res - i - j = chl index)
  tris <- list()
  bary_xy <- function(i_nap, j_doc) {
    f_nap <- i_nap / res
    f_doc <- j_doc / res
    c(f_nap + 0.5 * f_doc, f_doc * S3_TERN)
  }

  for (j in 0:(res - 1L)) {
    for (i in 0:(res - 1L - j)) {
      # Upward triangle: (i, j), (i+1, j), (i, j+1)
      p1 <- bary_xy(i, j); p2 <- bary_xy(i + 1L, j); p3 <- bary_xy(i, j + 1L)
      f_nap <- (i + 1 / 3) / res
      f_doc <- (j + 1 / 3) / res
      f_chl <- max(0, 1 - f_nap - f_doc)
      st_c <- bio_optical_state(max(1e-3, f_chl) * 0.45 / 0.0086, max(1e-3, f_doc) * 0.45 / 0.0732, max(1e-3, f_nap) * 0.45 / 0.0432)
      st_d <- bio_optical_state(max(1e-3, f_chl) * 3.20 / 0.0086, max(1e-3, f_doc) * 3.20 / 0.0732, max(1e-3, f_nap) * 3.20 / 0.0432)
      tris[[length(tris) + 1L]] <- list(
        x = c(p1[1], p2[1], p3[1]), y = c(p1[2], p2[2], p3[2]),
        col_c = st_c$col_srgb, col_d = st_d$col_srgb
      )
      # Downward triangle if inside simplex: (i+1, j), (i+1, j+1), (i, j+1)
      if (i + j + 2L <= res) {
        p4 <- bary_xy(i + 1L, j + 1L)
        f_nap2 <- (i + 2 / 3) / res
        f_doc2 <- (j + 2 / 3) / res
        f_chl2 <- max(0, 1 - f_nap2 - f_doc2)
        st_c2 <- bio_optical_state(max(1e-3, f_chl2) * 0.45 / 0.0086, max(1e-3, f_doc2) * 0.45 / 0.0732, max(1e-3, f_nap2) * 0.45 / 0.0432)
        st_d2 <- bio_optical_state(max(1e-3, f_chl2) * 3.20 / 0.0086, max(1e-3, f_doc2) * 3.20 / 0.0732, max(1e-3, f_nap2) * 3.20 / 0.0432)
        tris[[length(tris) + 1L]] <- list(
          x = c(p2[1], p4[1], p3[1]), y = c(p2[2], p4[2], p3[2]),
          col_c = st_c2$col_srgb, col_d = st_d2$col_srgb
        )
      }
    }
  }

  list(tris = tris)
}

#' Draw a high-resolution Full HD frame (1920 x 1080 @ 180 dpi) on the dark theme (`#0b111e`)
#' where the Ternary Diagram echoes the orientation, dashed 20/40/60/80 grid, angled ticks,
#' and angled arrow labels (`NAP ->`, `<- DOC`, `<- Chl a`) of the manuscript figures!
draw_composite_frame <- function(out_png, lake_img, gpu_renders, st, trail_df, tmesh, width = 1920, height = 1080) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(out_png, width = width, height = height, units = "px", res = 180, background = "#0b111e")
  } else {
    png(out_png, width = width, height = height, units = "px", res = 180, bg = "#0b111e", type = "quartz")
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = c(0, 1920), ylim = c(0, 1080))

  # ===========================================================================
  # 1. Left Panel (x: 0 -> 1140, y: 0 -> 1080): High-Res Native {sonicray} 3D Lake
  # ===========================================================================
  rasterImage(lake_img, xleft = 0, ybottom = 0, xright = 1140, ytop = 1080, interpolate = TRUE)

  # Crisp depth labels next to the 4 submerged fish & Secchi spheres
  tdf <- gpu_renders$tier_df
  for (i in seq_len(nrow(tdf))) {
    vis <- exp(-0.95 * st$kd * tdf$depth_m[i])
    if (vis > 0.20) {
      text(tdf$x_px[i], tdf$y_px[i] + 26, sprintf("%.1f m", tdf$depth_m[i]),
           col = rgb(1, 1, 1, pmin(1, vis * 1.15)), font = 2, cex = 0.65)
    }
  }

  # Clean top title bar on the 3D lake panel
  rect(24, 998, 1116, 1058, col = rgb(0.05, 0.08, 0.14, 0.82), border = rgb(1, 1, 1, 0.22), lwd = 1.4)
  text(46, 1028, st$phase, adj = c(0, 0.5), col = "white", font = 2, cex = 0.96)
  rect(888, 1006, 1102, 1050, col = st$hex, border = "white", lwd = 1.8)
  tcol_top <- if (mean(st$col_srgb) > 0.52) "#0f172a" else "#ffffff"
  text(995, 1028, sprintf("Secchi = %.1f m", st$z_sd), col = tcol_top, font = 2, cex = 0.86)

  # ===========================================================================
  # 2. Right Panel (x: 1140 -> 1920, y: 0 -> 1080): Dark Studio Panel (`#0b111e`)
  #    with Ternary Diagram Echoing the Manuscript Style (DOC top, Chl a left, NAP right)
  # ===========================================================================
  rect(1140, 0, 1920, 1080, col = "#0b111e", border = "#1e293b", lwd = 1.5)

  # Top: 3 Constituent Bars (y: 830 -> 1055) on dark card
  rect(1162, 780, 1898, 1058, col = "#131c2e", border = "#26354d", lwd = 1.3)
  text(1184, 1034, "Constituent Titration & Attenuation Partition", adj = c(0, 0.5), col = "#e2e8f0", font = 2, cex = 0.86)

  draw_dark_bar <- function(yc, name, val_txt, pct_txt, frac, col_bar) {
    text(1184, yc + 13, name, adj = c(0, 0.5), col = "#cbd5e1", font = 2, cex = 0.72)
    text(1876, yc + 13, paste0(val_txt, "  (", pct_txt, ")"), adj = c(1, 0.5), col = "white", font = 2, cex = 0.70)
    rect(1184, yc - 12, 1876, yc, col = "#090e17", border = "#334155", lwd = 1.0)
    w <- pmin(1, pmax(0.015, frac)) * (1876 - 1184)
    rect(1184, yc - 12, 1184 + w, yc, col = col_bar, border = NA)
  }
  draw_dark_bar(984, "Chl a (Algae)", sprintf("%.1f ug/L", st$chl_a), sprintf("%.0f%% kd", 100 * st$f_chl), (st$chl_a / 135)^0.60, "#22c55e")
  draw_dark_bar(918, "DOC (Brown)",   sprintf("%.1f mg/L", st$doc),   sprintf("%.0f%% kd", 100 * st$f_doc), (st$doc   / 28)^0.60,  "#d97706")
  draw_dark_bar(852, "NAP / TSS",     sprintf("%.1f mg/L", st$nap),   sprintf("%.0f%% kd", 100 * st$f_nap), (st$nap   / 75)^0.60,  "#facc15")

  # 100% Partition Bar
  x0 <- 1184; xw <- 692
  w1 <- st$f_chl * xw; w2 <- st$f_doc * xw
  rect(x0, 796, x0 + w1, 818, col = "#22c55e", border = NA)
  rect(x0 + w1, 796, x0 + w1 + w2, 818, col = "#b45309", border = NA)
  rect(x0 + w1 + w2, 796, x0 + xw, 818, col = "#eab308", border = NA)
  rect(x0, 796, x0 + xw, 818, col = NA, border = "white", lwd = 1.1)

  # ===========================================================================
  # Center: Ternary Diagram Echoing `tern_frame()` from `R/ternary.R`
  #   - Top vertex = DOC (100), Bottom-left = Chl a (100), Bottom-right = NAP (100)
  #   - Smooth continuous triangular color mesh (Clear -> Dark regime)
  #   - Dashed 20, 40, 60, 80 internal grid lines
  #   - Angled 20, 40, 60, 80, 100 tick labels along all 3 edges
  #   - Angled axis titles: `<- Chl a`, `<- DOC`, `NAP ->`
  # ===========================================================================
  rect(1162, 190, 1898, 764, col = "#131c2e", border = "#26354d", lwd = 1.3)

  regime_txt <- if (st$darkness < 0.35) "Clear Water Regime (High Visibility)" else "Dark / Turbid Regime (Low Visibility)"
  regime_col <- if (st$darkness < 0.35) "#38bdf8" else "#f59e0b"
  rect(1235, 712, 1825, 748, col = "#090e17", border = regime_col, lwd = 1.6)
  text(1530, 730, regime_txt, col = regime_col, font = 2, cex = 0.80)

  # Map ternary unit coordinates (x in [0, 1], y in [0, sqrt(3)/2]) to pixel box:
  cx_left <- 1285; cx_right <- 1775
  tw <- cx_right - cx_left # 490 px wide equilateral triangle!
  cy_bot <- 254
  tx <- function(x) cx_left + x * tw
  ty <- function(y) cy_bot + y * tw

  # 1. Draw seamless continuous water-color background inside the triangle (blended by `st$darkness`)
  dk <- st$darkness
  n_tri <- length(tmesh$tris)
  poly_x <- numeric(4 * n_tri)
  poly_y <- numeric(4 * n_tri)
  poly_cols <- character(n_tri)
  for (m in seq_len(n_tri)) {
    tr <- tmesh$tris[[m]]
    cm <- (1 - dk) * tr$col_c + dk * tr$col_d
    poly_cols[m] <- rgb(cm[1], cm[2], cm[3])
    idx4 <- (m - 1L) * 4L
    poly_x[idx4 + 1:3] <- tx(tr$x); poly_x[idx4 + 4L] <- NA
    poly_y[idx4 + 1:3] <- ty(tr$y); poly_y[idx4 + 4L] <- NA
  }
  polygon(poly_x[-length(poly_x)], poly_y[-length(poly_y)], col = poly_cols, border = poly_cols, lwd = 0.4)

  # 2. Draw dashed 20%, 40%, 60%, 80% internal grid lines echoing `tern_frame()` (`R/ternary.R`)
  v_grid <- c(0.2, 0.4, 0.6, 0.8)
  for (vg in v_grid) {
    # Constant NAP lines
    segments(tx(vg), ty(0), tx(vg + (1 - vg) / 2), ty((1 - vg) * S3_TERN), col = rgb(1, 1, 1, 0.36), lwd = 1.2, lty = 2)
    # Constant DOC lines (horizontal)
    segments(tx(vg / 2), ty(vg * S3_TERN), tx(1 - vg / 2), ty(vg * S3_TERN), col = rgb(1, 1, 1, 0.36), lwd = 1.2, lty = 2)
    # Constant Chl a lines
    segments(tx(1 - vg), ty(0), tx((1 - vg) / 2), ty((1 - vg) * S3_TERN), col = rgb(1, 1, 1, 0.36), lwd = 1.2, lty = 2)
  }

  # 3. Outer triangle border
  polygon(c(tx(0), tx(1), tx(0.5)), c(ty(0), ty(0), ty(S3_TERN)), col = NA, border = "#e2e8f0", lwd = 2.4)

  # 4. Angled tick labels (20, 40, 60, 80, 100) echoing `tern_frame()` (`R/ternary.R`)
  tk <- c(20, 40, 60, 80, 100)
  f_tk <- tk / 100
  for (idx in seq_along(tk)) {
    fv <- f_tk[idx]; lbl <- as.character(tk[idx])
    # NAP along bottom edge (left -> right, angle = 60, hjust = 1)
    text(tx(fv + 0.015), ty(-0.032), lbl, srt = 60, adj = c(1, 0.5), col = "#cbd5e1", cex = 0.68)
    # DOC along right edge (bottom -> top, angle = 0, hjust = 0)
    text(tx(1 - fv / 2 + 0.026), ty(fv * S3_TERN), lbl, srt = 0, adj = c(0, 0.5), col = "#cbd5e1", cex = 0.68)
    # Chl a along left edge (top -> bottom-left, angle = -60, hjust = 1)
    text(tx((1 - fv) / 2 - 0.026), ty((1 - fv) * S3_TERN + 0.01), lbl, srt = -60, adj = c(1, 0.5), col = "#cbd5e1", cex = 0.68)
  }

  # 5. Angled axis titles (`NAP ->`, `<- DOC`, `<- Chl a`) echoing `tern_frame()` (`R/ternary.R`)
  text(tx(0.50), ty(-0.145), "NAP \u2192",   srt =   0, adj = c(0.5, 0.5), col = "#facc15", font = 2, cex = 0.88)
  text(tx(0.84), ty( 0.505), "\u2190 DOC",   srt = -60, adj = c(0.5, 0.5), col = "#f59e0b", font = 2, cex = 0.88)
  text(tx(0.16), ty( 0.505), "\u2190 Chl a", srt =  60, adj = c(0.5, 0.5), col = "#4ade80", font = 2, cex = 0.88)

  # 6. Live titration trajectory & current point on the ternary diagram
  if (nrow(trail_df) > 1) {
    lines(tx(trail_df$x), ty(trail_df$y), col = rgb(1, 1, 1, 0.72), lwd = 2.4)
  }
  cur_xy <- tern_xy_paper(st$f_chl, st$f_doc, st$f_nap)
  points(tx(cur_xy[1, "x"]), ty(cur_xy[1, "y"]), pch = 21, bg = st$hex, col = "white", cex = 2.6, lwd = 2.8)
  points(tx(cur_xy[1, "x"]), ty(cur_xy[1, "y"]), pch = 1,  col = "#090e17", cex = 2.0, lwd = 1.6)

  # ===========================================================================
  # Bottom: Visual Underwater Visibility Meter (0.4 m, 0.9 m, 1.5 m, 2.1 m)
  # ===========================================================================
  rect(1162, 22, 1898, 174, col = "#131c2e", border = "#26354d", lwd = 1.3)
  text(1530, 152, "Submerged Fish & Secchi Sphere Visibility by Depth", col = "#e2e8f0", font = 2, cex = 0.76)

  f_depths <- c(0.4, 0.9, 1.5, 2.1)
  x_centers <- seq(1266, 1794, length.out = 4)
  for (k in seq_along(f_depths)) {
    zm <- f_depths[k]
    xc <- x_centers[k]
    vis <- exp(-0.95 * st$kd * zm)
    rect(xc - 66, 36, xc + 66, 130, col = "#090e17",
         border = if (vis > 0.25) "#38bdf8" else "#334155", lwd = if (vis > 0.25) 1.8 else 1.1)
    t_seq <- seq(0, 2 * pi, length.out = 32)
    bx <- xc + 24 * cos(t_seq)
    by <- 94 + 7.5 * sin(t_seq) * (0.75 + 0.35 * cos(t_seq))
    tx_p <- xc - 24 * c(0.75, 1.25, 1.05, 1.25)
    ty_p <- 94 + 7.5 * c(0, 0.9, 0, -0.9)
    f_alpha <- pmax(0.10, pmin(1.0, vis))
    polygon(tx_p, ty_p, col = rgb(0.96, 0.42, 0.28, f_alpha), border = NA)
    polygon(bx, by, col = rgb(0.94, 0.84, 0.74, f_alpha), border = rgb(1, 1, 1, f_alpha * 0.6), lwd = 0.9)
    text(xc, 66, sprintf("%.1f m", zm), col = "white", font = 2, cex = 0.68)
    status_txt <- if (vis > 0.50) "Clear" else if (vis > 0.22) "Faint" else "Hidden"
    status_col <- if (vis > 0.50) "#4ade80" else if (vis > 0.22) "#facc15" else "#64748b"
    text(xc, 46, status_txt, col = status_col, font = 2, cex = 0.62)
  }
}

#' Anti-aliased 2D box-filter downsampler so the contact sheet has zero nearest-neighbor jaggies!
downsample_aa <- function(img, sh, sw) {
  fh <- dim(img)[1]; fw <- dim(img)[2]
  # 2x2 supersampled bilinear average
  y0 <- pmin(fh, pmax(1L, floor(seq(1, fh - 0.5, length.out = sh))))
  y1 <- pmin(fh, y0 + 1L)
  x0 <- pmin(fw, pmax(1L, floor(seq(1, fw - 0.5, length.out = sw))))
  x1 <- pmin(fw, x0 + 1L)
  0.25 * (img[y0, x0, 1:3] + img[y1, x0, 1:3] + img[y0, x1, 1:3] + img[y1, x1, 1:3])
}

#' Build the high-resolution 6-keyframe contact sheet and full 28-second MP4
build_sonicray_titration_visuals <- function(
  out_mp4   = "figs/sonicray_optical_titration.mp4",
  out_sheet = "figs/sonicray_optical_titration_sheet.png",
  fps       = 12L,
  duration_s = 28.0,
  sheet_only = FALSE
) {
  t0 <- Sys.time()
  cache_rds <- file.path(tempdir(), "sonicray_gpu_renders_1140x1080_4096spp.rds")
  if (file.exists(cache_rds)) {
    message("[sonicray_lake] Loading cached GPU water renders from ", cache_rds)
    gpu_renders <- readRDS(cache_rds)
  } else {
    scene <- build_sonicray_water_scene(n = 240L, tex = 2048L)
    gpu_renders <- render_native_water_states(scene, width = 1140L, height = 1080L, amb_spp = 768L, light_spp = 4096L)
    saveRDS(gpu_renders, cache_rds)
  }
  tmesh <- build_paper_ternary_mesh(res = 36L)

  times <- seq(0, duration_s, by = 1 / fps)
  states <- lapply(times, titration_schedule, total_s = duration_s)
  xy_all <- do.call(rbind, lapply(states, function(s) tern_xy_paper(s$f_chl, s$f_doc, s$f_nap)))
  trail_all <- tibble(x = xy_all[, "x"], y = xy_all[, "y"])

  key_times <- c(2.0, 6.2, 10.6, 14.8, 20.2, 25.5)
  key_idx <- vapply(key_times, function(kt) which.min(abs(times - kt)), 1L)

  frames_dir <- file.path(tempdir(), "sonicray_titration_frames_hq")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

  render_indices <- if (sheet_only) key_idx else seq_along(times)
  message(sprintf("[sonicray_lake] Compositing %d Full-HD frames (1920x1080)...", length(render_indices)))

  state_names <- names(gpu_renders$states)
  mask_mat <- gpu_renders$states[[1]]$mask

  for (i in render_indices) {
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
    draw_composite_frame(f_png, img3d, gpu_renders, st_i, trail_all[1:i, , drop = FALSE], tmesh)
  }

  # Assemble high-resolution anti-aliased 3x2 contact sheet (2880 x 2430 px)
  key_imgs <- lapply(key_idx, function(idx) png::readPNG(file.path(frames_dir, sprintf("frame_%04d.png", idx))))
  fh <- dim(key_imgs[[1]])[1]; fw <- dim(key_imgs[[1]])[2]
  sh <- round(fh * 0.75); sw <- round(fw * 0.75) # 1440 x 810 per tile => 2916 x 2478 total sheet!
  pad <- 12L
  sheet <- array(0.05, c(3 * sh + 4 * pad, 2 * sw + 3 * pad, 3))
  for (k in seq_along(key_imgs)) {
    r0 <- ((k - 1) %/% 2) * (sh + pad) + pad
    c0 <- ((k - 1) %% 2)  * (sw + pad) + pad
    sheet[r0 + seq_len(sh), c0 + seq_len(sw), ] <- downsample_aa(key_imgs[[k]], sh, sw)
  }
  png::writePNG(sheet, out_sheet)
  message("[sonicray_lake] Wrote high-resolution keyframe sheet: ", out_sheet)

  if (!sheet_only) {
    ff <- Sys.which("ffmpeg")
    cmd <- sprintf(
      '"%s" -y -framerate %d -i "%s/frame_%%04d.png" -c:v libx264 -preset fast -crf 17 -pix_fmt yuv420p -movflags +faststart "%s"',
      ff, fps, frames_dir, out_mp4
    )
    system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    message(sprintf("[sonicray_lake] Wrote %s (%.2f MB) in %.1f s total.",
                    out_mp4, file.info(out_mp4)$size / 1e6, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  invisible(list(mp4 = out_mp4, sheet = out_sheet))
}
