# -----------------------------------------------------------------------------
# sonicray_lake.R -- {sonicray} Metal GPU 3D `render_water()` Cutaway Diorama
#                    with 3D Path-Traced Fish, Underwater Bathymetry, Secchi Disk,
#                    and Dynamic Clear-vs-Dark Ternary Diagrams.
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

#' Combine multiple 3D triangle meshes into a single {sonicray} mesh
combine_sr_meshes <- function(meshes) {
  meshes <- Filter(Negate(is.null), meshes)
  v_list <- list()
  n_list <- list()
  u_list <- list()
  i_list <- list()
  v_offset <- 0L

  for (m in meshes) {
    v_list[[length(v_list) + 1L]] <- m$vertices
    n_list[[length(n_list) + 1L]] <- m$normals
    u_list[[length(u_list) + 1L]] <- m$uvs
    i_list[[length(i_list) + 1L]] <- m$indices + v_offset
    v_offset <- v_offset + nrow(m$vertices)
  }

  tri <- do.call(rbind, i_list)
  storage.mode(tri) <- "integer"
  list(
    vertices = do.call(rbind, v_list),
    normals  = do.call(rbind, n_list),
    uvs      = do.call(rbind, u_list),
    indices  = tri
  )
}

#' Build a closed 3D fusiform trout/salmon mesh in world space (cx, cy, cz)
#' with body length `L`, swimming yaw angle `yaw_rad`, and UV coordinates mapped
#' into the fish region of the texture atlas (`u_fish`, `v in [0.05, 0.95]`).
build_3d_fish_mesh <- function(cx, cy, cz, L = 2.3, yaw_rad = 0.0, u_fish = 0.94) {
  n_len <- 22L
  n_ang <- 16L
  s_seq <- seq(0, 1, length.out = n_len) # 0 = nose, 1 = tail tip
  theta_seq <- seq(0, 2 * pi, length.out = n_ang + 1L)[-(n_ang + 1L)]

  r_y <- ifelse(
    s_seq < 0.78,
    0.19 * sin(pi * (s_seq / 0.78)^0.72),
    0.045 + 0.16 * ((s_seq - 0.78) / 0.22)^1.1
  )
  r_z <- ifelse(
    s_seq < 0.78,
    0.11 * sin(pi * (s_seq / 0.78)^0.68),
    0.020 * (1 - 0.7 * (s_seq - 0.78) / 0.22)
  )
  dorsal_boost <- ifelse(s_seq >= 0.34 & s_seq <= 0.56, 0.09 * sin(pi * (s_seq - 0.34) / 0.22), 0.0)

  verts <- matrix(0.0, nrow = n_len * n_ang + 2L, ncol = 3L)
  uvs   <- matrix(0.0, nrow = n_len * n_ang + 2L, ncol = 2L)

  idx_mat <- matrix(0L, nrow = n_len, ncol = n_ang)
  k <- 1L
  for (i in seq_len(n_len)) {
    s <- s_seq[i]
    lx <- (0.5 - s) * L
    for (j in seq_len(n_ang)) {
      th <- theta_seq[j]
      sy <- sin(th)
      cz_local <- cos(th)
      ly <- r_y[i] * L * sy + (if (sy > 0.65) dorsal_boost[i] * L * ((sy - 0.65) / 0.35)^1.4 else 0.0)
      if (s > 0.88 && abs(sy) < 0.45) {
        lx <- lx + 0.07 * L * (1 - abs(sy) / 0.45) * ((s - 0.88) / 0.12)
      }
      lz <- r_z[i] * L * cz_local

      verts[k, ] <- c(lx, ly, lz)
      v_coord <- 0.50 + 0.42 * sy
      uvs[k, ] <- c(u_fish, pmin(pmax(v_coord, 0.05), 0.95))
      idx_mat[i, j] <- k
      k <- k + 1L
    }
  }
  i_nose <- k; verts[i_nose, ] <- c(0.51 * L, 0, 0); uvs[i_nose, ] <- c(u_fish, 0.65); k <- k + 1L
  i_tail <- k; verts[i_tail, ] <- c(-0.46 * L, 0, 0); uvs[i_tail, ] <- c(u_fish, 0.55)

  co <- cos(yaw_rad); si <- sin(yaw_rad)
  vx <- verts[, 1] * co - verts[, 3] * si + cx
  vy <- verts[, 2] + cy
  vz <- verts[, 1] * si + verts[, 3] * co + cz
  verts <- cbind(vx, vy, vz)

  tris <- list()
  for (i in seq_len(n_len - 1L)) {
    for (j in seq_len(n_ang)) {
      j_next <- if (j == n_ang) 1L else j + 1L
      a <- idx_mat[i, j]
      b <- idx_mat[i, j_next]
      c <- idx_mat[i + 1L, j]
      d <- idx_mat[i + 1L, j_next]
      tris[[length(tris) + 1L]] <- rbind(c(a, c, b), c(b, c, d))
    }
  }
  for (j in seq_len(n_ang)) {
    j_next <- if (j == n_ang) 1L else j + 1L
    tris[[length(tris) + 1L]] <- rbind(
      c(i_nose, idx_mat[1L, j], idx_mat[1L, j_next]),
      c(idx_mat[n_len, j], i_tail, idx_mat[n_len, j_next])
    )
  }
  tri <- do.call(rbind, tris) - 1L
  storage.mode(tri) <- "integer"
  cpp_normals <- getFromNamespace("cpp_vertex_normals", "sonicray")
  list(vertices = verts, normals = cpp_normals(verts, tri), uvs = uvs, indices = tri)
}

#' Build vertical cutaway skirt walls around a heightfield grid `Y_top` down to `Y_bot`
#' along all 4 boundaries of the 3D diorama block (`xr`, `zr`).
build_block_side_walls <- function(Y_top, Y_bot, xr, zr,
                                   u_front = 0.86, u_right = 0.88, u_back = 0.86, u_left = 0.88,
                                   v_top = 0.80, v_bot = 0.20) {
  nr <- nrow(Y_top); nc <- ncol(Y_top)
  xs <- seq(xr[1], xr[2], length.out = nc)
  zs <- seq(zr[1], zr[2], length.out = nr)
  cpp_normals <- getFromNamespace("cpp_vertex_normals", "sonicray")

  make_strip <- function(x_vec, z_vec, yt_vec, yb_vec, u_val, flip = FALSE) {
    m <- length(x_vec)
    V <- rbind(
      cbind(x_vec, yt_vec, z_vec),
      cbind(x_vec, yb_vec, z_vec)
    )
    U <- rbind(
      cbind(rep(u_val, m), rep(v_top, m)),
      cbind(rep(u_val, m), rep(v_bot, m))
    )
    a <- seq_len(m - 1L)
    b <- a + 1L
    c <- a + m
    d <- b + m
    Tri <- if (!flip) rbind(cbind(a, c, b), cbind(b, c, d)) else rbind(cbind(a, b, c), cbind(b, d, c))
    Tri <- Tri - 1L
    storage.mode(Tri) <- "integer"
    list(vertices = V, normals = cpp_normals(V, Tri), uvs = U, indices = Tri)
  }

  w_front <- make_strip(xs, rep(zs[nr], nc), Y_top[nr, ], Y_bot[nr, ], u_val = u_front, flip = FALSE)
  w_back  <- make_strip(xs, rep(zs[1L], nc), Y_top[1L, ], Y_bot[1L, ], u_val = u_back,  flip = TRUE)
  w_right <- make_strip(rep(xs[nc], nr), zs, Y_top[, nc], Y_bot[, nc], u_val = u_right, flip = TRUE)
  w_left  <- make_strip(rep(xs[1L], nr), zs, Y_top[, 1L], Y_bot[, 1L], u_val = u_left,  flip = FALSE)

  combine_sr_meshes(list(w_front, w_back, w_right, w_left))
}

#' Build the complete 3D `rayshader::render_water()` style cutaway lake diorama scene:
#'   - Full 3D diorama block (`xr = c(-13.2, 13.2)`, `zr = c(-13.2, 13.2)`) framed cleanly in perspective
#'   - Mountain amphitheater on the back-left descending into a deep 3D sculpted lake basin
#'     that cuts cleanly through the front (+z) and right (+x) cutaway block faces!
#'   - 3D path-traced swimming fish at 4 depths (0.8 m, 2.2 m, 4.2 m, 6.8 m) casting real shadows on the lakebed
build_lake_cutaway_scene <- function(seed = 42L, n = 400L, tex = 2048L) {
  set.seed(seed)
  fbm  <- getFromNamespace("fbm_noise", "sonicray")
  rmat <- getFromNamespace("resize_mat", "sonicray")
  cpp_normals <- getFromNamespace("cpp_vertex_normals", "sonicray")

  # r = 0 (back -z) to r = 1 (front +z cutaway face)
  # cx = 0 (left -x) to cx = 1 (right +x cutaway face)
  r  <- matrix(rep(seq(0, 1, length.out = n), times = n), n, n)
  cx <- matrix(rep(seq(0, 1, length.out = n), each = n), n, n)

  fb1 <- fbm(n, 2.4)
  fb2 <- fbm(n, 1.9)
  fb3 <- fbm(n, 3.2)
  ridge <- 1 - abs(2 * fbm(n, 2.1) - 1)

  # The lake basin fills the front-right 70% of the 3D block and intersects BOTH the front (r = 1)
  # and right (cx = 1) cutaway faces, creating the signature `rayshader::render_water()` 3D cutaway slab!
  Ell <- sqrt(((cx - 0.72) / 0.62)^2 + ((r - 0.72) / 0.62)^2) +
         0.07 * (fb2 - 0.5) + 0.025 * sin(5 * pi * cx) * cos(4 * pi * r)

  in_lake <- (Ell < 1.0)

  # Bathymetric depth in meters (0 m at shore -> 2 m littoral shelf -> 4.5 m slope -> 8.8 m deep canyon at front-right corner)
  s_in <- pmax(0, 1.0 - Ell)
  boulder_field <- 0.55 * pmax(0, fb3 - 0.56)^1.4 * 7.5
  bathy_m <- matrix(0, n, n)
  bathy_m[in_lake] <- pmax(
    0.15,
    8.8 * (0.24 * smoothstep_sr(0.0, 0.25, s_in[in_lake]) +
           0.76 * smoothstep_sr(0.18, 0.85, s_in[in_lake])^1.3) -
    1.0 * boulder_field[in_lake] * smoothstep_sr(0.08, 0.45, s_in[in_lake])
  )

  # Mountain amphitheater wrapping around the back-left shoreline (Ell >= 1.0)
  s_out <- pmax(0, Ell - 1.0)
  # Taper extreme back-corner height so the highest peak frames cleanly below the top title bar
  corner_taper <- 1.0 - 0.36 * exp(-((cx / 0.28)^2 + (r / 0.28)^2))
  mtn_env <- (0.36 + 1.05 * (1 - r)^1.15 + 1.00 * (1 - cx)^1.15) * corner_taper
  stream_d <- abs((cx - 0.24) - 0.55 * (r - 0.22))
  gorge <- 0.35 + 0.65 * smoothstep_sr(0.02, 0.14, stream_d)

  # World-space Y coordinates (1 m water depth = 0.52 world units so 8.8 m = 4.58 units deep!)
  water_scale <- 0.52
  Y_mtn  <- pmin(5.2, (5.8 * mtn_env * gorge * (0.56 * fb1 + 0.44 * ridge) + 0.45 * fb2) * smoothstep_sr(0.0, 0.50, s_out)^1.12)
  Y_bed  <- matrix(Y_mtn * (!in_lake) - bathy_m * water_scale * in_lake, n, n)
  Y_surf <- matrix(Y_mtn * (!in_lake) + 0.0 * in_lake, n, n)
  Y_base <- matrix(-5.6, n, n)

  xr <- c(-13.2, 13.2)
  zr <- c(-13.2, 13.2)

  make_hf_mesh <- function(Y_mat, u_max = 0.84) {
    nr <- nrow(Y_mat); nc <- ncol(Y_mat)
    xs <- seq(xr[1], xr[2], length.out = nc)
    zs <- seq(zr[1], zr[2], length.out = nr)
    idx <- matrix(seq_len(nr * nc), nr, nc)
    V <- cbind(rep(xs, each = nr), as.vector(Y_mat), rep(zs, times = nc))
    U <- cbind(rep(seq(0.01, u_max, length.out = nc), each = nr),
               1 - rep(seq(0.01, 0.99, length.out = nr), times = nc))
    a <- as.vector(idx[-nr, -nc]); b <- as.vector(idx[-nr, -1])
    c <- as.vector(idx[-1, -nc]);  d <- as.vector(idx[-1, -1])
    Tri <- rbind(cbind(a, c, b), cbind(b, c, d)) - 1L
    storage.mode(Tri) <- "integer"
    list(vertices = V, normals = cpp_normals(V, Tri), uvs = U, indices = Tri)
  }

  # 1. Mesh for Pass 1: 3D Mountain + 3D Carved Underwater Lakebed + 3D Diorama Base Walls (from Y_bed down to Y_base)
  mesh_top_bed   <- make_hf_mesh(Y_bed, u_max = 0.84)
  mesh_base_wall <- build_block_side_walls(Y_bed, Y_base, xr, zr,
                                           u_front = 0.855, u_right = 0.865, u_back = 0.855, u_left = 0.865,
                                           v_top = 0.75, v_bot = 0.25)
  mesh_bed_no_fish <- combine_sr_meshes(list(mesh_top_bed, mesh_base_wall))

  # 2. Build 3D Fish School at 4 Depth Tiers (0.8 m, 2.2 m, 4.2 m, 6.8 m) inside the 3D cutaway lake!
  fish_specs <- list(
    # Tier 1: Shallow littoral / epilimnion (0.8 m & 1.0 m)
    list(x =  0.2, z =  1.8, z_m = 0.8, L = 2.2, yaw =  0.22, u_fish = 0.925),
    list(x =  4.5, z = -0.8, z_m = 1.0, L = 1.9, yaw = -0.30, u_fish = 0.928),
    # Tier 2: Upper pelagic (2.2 m & 2.5 m)
    list(x =  2.2, z =  5.2, z_m = 2.2, L = 2.3, yaw =  0.18, u_fish = 0.940),
    list(x =  7.2, z =  3.2, z_m = 2.5, L = 2.0, yaw = -0.22, u_fish = 0.943),
    # Tier 3: Mid-depth reef (4.2 m & 4.5 m)
    list(x =  4.2, z =  8.4, z_m = 4.2, L = 2.4, yaw =  0.15, u_fish = 0.955),
    list(x =  8.8, z =  6.8, z_m = 4.5, L = 2.1, yaw = -0.25, u_fish = 0.958),
    # Tier 4: Deep pelagic canyon near front-right cutaway corner (6.8 m & 7.1 m)
    list(x =  6.8, z = 10.6, z_m = 6.8, L = 2.5, yaw =  0.20, u_fish = 0.970),
    list(x = 10.2, z =  9.4, z_m = 7.1, L = 2.1, yaw =  0.10, u_fish = 0.973)
  )
  fish_meshes <- lapply(fish_specs, function(f) {
    build_3d_fish_mesh(cx = f$x, cy = -f$z_m * water_scale, cz = f$z, L = f$L, yaw_rad = f$yaw, u_fish = f$u_fish)
  })
  mesh_bed_with_fish <- combine_sr_meshes(c(list(mesh_top_bed, mesh_base_wall), fish_meshes))

  # 3. Mesh for Pass 3 (`gpu_water_slab`): 3D Water Surface + 3D Vertical Cutaway Water Walls (`render_water` glass block!)
  wave_h <- 0.035 * (fbm(n, 7.5) - 0.5) * in_lake * smoothstep_sr(0.0, 0.08, s_in)
  Y_water_top <- Y_surf + wave_h
  mesh_top_water <- make_hf_mesh(Y_water_top, u_max = 0.84)
  # Notice: On the vertical cutaway walls, we build the WATER cutaway wall strictly from `Y_water_top` down to `Y_bed`!
  # Where there is no lake (`!in_lake`), `Y_water_top == Y_bed`, so the strip has zero height;
  # where there IS lake (`in_lake`), the strip is the exact 3D vertical glass cutaway wall of the water column!
  mesh_water_cut_wall <- build_block_side_walls(Y_water_top, Y_bed, xr, zr,
                                                u_front = 0.885, u_right = 0.900, u_back = 0.885, u_left = 0.900,
                                                v_top = 0.95, v_bot = 0.05)
  mesh_water_slab <- combine_sr_meshes(list(mesh_top_water, mesh_water_cut_wall, mesh_base_wall))

  # ---------------------------------------------------------------------------
  # Build the 2048x2048 Albedo Texture Atlas (`albedo_bed`)
  # ---------------------------------------------------------------------------
  e <- matrix(pmax(0, Y_bed) / max(Y_bed), n, n)
  gy <- Y_bed[c(2:n, n), ] - Y_bed[c(1, 1:(n - 1)), ]
  gx <- Y_bed[, c(2:n, n)] - Y_bed[, c(1, 1:(n - 1))]
  slope <- matrix(pmin(sqrt(gx^2 + gy^2) / 0.38, 1), n, n)

  mix3 <- function(a, b, w) array(c(outer(1 - w, a) + outer(w, b)), c(dim(w), 3))
  col_land <- mix3(c(0.09, 0.24, 0.11), c(0.27, 0.28, 0.18), matrix(smoothstep_sr(0.03, 0.28, e), n, n))
  rockw <- matrix(pmax(smoothstep_sr(0.24, 0.58, e), slope * 0.85), n, n)
  for (k in 1:3) col_land[, , k] <- col_land[, , k] * (1 - rockw) + c(0.47, 0.44, 0.41)[k] * rockw * (0.82 + 0.36 * fb2)
  snow <- matrix(smoothstep_sr(0.54, 0.74, e) * (1 - 0.50 * slope), n, n)
  for (k in 1:3) col_land[, , k] <- col_land[, , k] * (1 - snow) + c(0.94, 0.96, 1.00)[k] * snow

  # Underwater 3D lakebed texture: warm littoral sand near shore -> natural cobble/boulder basin
  d_norm <- matrix(pmin(1, bathy_m / 8.8), n, n)
  col_bed <- mix3(c(0.76, 0.68, 0.52), c(0.32, 0.38, 0.44), matrix(smoothstep_sr(0.02, 0.65, d_norm), n, n))
  boulder_shade <- 0.76 + 0.48 * fb3
  for (k in 1:3) col_bed[, , k] <- pmin(pmax(col_bed[, , k] * boulder_shade, 0.06), 0.92)

  col_full <- col_land
  for (k in 1:3) col_full[, , k] <- col_land[, , k] * (!in_lake) + col_bed[, , k] * in_lake

  u_split1 <- round(0.84 * tex)
  u_split2 <- round(0.91 * tex)
  albedo_bed <- array(0, c(tex, tex, 3))
  for (k in 1:3) {
    albedo_bed[, 1:u_split1, k] <- pmin(pmax(rmat(col_full[, , k], tex, u_split1), 0), 1)^2.2
    # Diorama base wall: front face slightly brighter than right face for crisp 3D block bevel!
    c_front <- round(0.855 * tex); c_right <- round(0.865 * tex)
    albedo_bed[, (u_split1 + 1L):u_split2, k] <- c(0.18, 0.20, 0.24)[k]^2.2
    albedo_bed[, max(u_split1 + 1L, c_right - 12L):u_split2, k] <- c(0.12, 0.14, 0.17)[k]^2.2
  }
  # Fish coloration along rows of columns (u_split2 + 1L):tex
  v_seq <- seq(1, 0, length.out = tex)
  fish_r <- ifelse(v_seq > 0.62, 0.20, ifelse(v_seq > 0.42, 0.98, 0.95))
  fish_g <- ifelse(v_seq > 0.62, 0.35, ifelse(v_seq > 0.42, 0.30, 0.92))
  fish_b <- ifelse(v_seq > 0.62, 0.28, ifelse(v_seq > 0.42, 0.26, 0.88))
  for (col_i in (u_split2 + 1L):tex) {
    albedo_bed[, col_i, 1] <- fish_r^2.2
    albedo_bed[, col_i, 2] <- fish_g^2.2
    albedo_bed[, col_i, 3] <- fish_b^2.2
  }

  # Uniform neutral albedo for the water surface lighting pass (prevents any bathy_m spotlight artifact!)
  albedo_wat_neutral <- albedo_bed
  for (k in 1:3) {
    albedo_wat_neutral[, 1:u_split1, k] <- pmin(pmax(rmat((col_land[, , k] * (!in_lake) + 0.50 * in_lake), tex, u_split1), 0), 1)^2.2
  }

  # Exact Geometry/Depth Metadata Texture (`meta_tex`, 2048x2048x3):
  #   Channel 1 (R): exact bathymetric / object depth in meters (`z_m / 10.0`)
  #   Channel 2 (G): surface class:
  #     - `1.0` = top water surface / underwater lakebed (`in_lake`)
  #     - `0.0` = terrestrial land (`!in_lake`)
  #     - `0.25` = solid rock diorama base wall (`u in [0.85, 0.87]`)
  #     - `0.55` = front vertical water cutaway wall (`u ~ 0.885`)
  #     - `0.75` = right vertical water cutaway wall (`u ~ 0.900`)
  #   Channel 3 (B): `is_fish` mask (`1.0` on 3D fish mesh, `0.0` elsewhere; on cutaway wall encodes vertical frac `1 - v`)
  meta_tex <- array(0, c(tex, tex, 3))
  meta_tex[, 1:u_split1, 1] <- pmin(pmax(rmat(bathy_m / 10.0, tex, u_split1), 0), 1)
  meta_tex[, 1:u_split1, 2] <- pmin(pmax(rmat(in_lake * 1.0, tex, u_split1), 0), 1)
  meta_tex[, 1:u_split1, 3] <- 0.0

  # Base rock wall columns (`u in [0.84, 0.875]`)
  c_base_end <- round(0.875 * tex)
  meta_tex[, (u_split1 + 1L):c_base_end, 1] <- 0.0
  meta_tex[, (u_split1 + 1L):c_base_end, 2] <- 0.25
  meta_tex[, (u_split1 + 1L):c_base_end, 3] <- 0.0

  # Vertical water cutaway glass walls (`u in [0.875, 0.91]`):
  #   Front glass wall (`u ~ 0.885`) -> G = 0.55
  #   Right glass wall (`u ~ 0.900`) -> G = 0.75
  #   B channel encodes vertical position on the glass wall (`v_seq`: 0 at top water surface -> 1 at lakebed bottom)
  c_cut_mid <- round(0.892 * tex)
  meta_tex[, (c_base_end + 1L):c_cut_mid, 2] <- 0.55
  meta_tex[, (c_cut_mid + 1L):u_split2,   2] <- 0.75
  for (col_i in (c_base_end + 1L):u_split2) {
    meta_tex[, col_i, 3] <- 1 - v_seq
  }

  # Fish columns (`u_split2+1 .. tex`): encode exact fish depth `z_m / 10.0` in R and `is_fish = 1.0` in B
  for (f in fish_specs) {
    c_idx <- max(u_split2 + 1L, min(tex, round(f$u_fish * tex)))
    c_rng <- max(u_split2 + 1L, c_idx - 6L):min(tex, c_idx + 6L)
    meta_tex[, c_rng, 1] <- f$z_m / 10.0
    meta_tex[, c_rng, 2] <- 1.0
    meta_tex[, c_rng, 3] <- 1.0
  }

  # Elevated 3/4 isometric camera framing the ENTIRE 3D `render_water()` cutaway block inside the left panel!
  cam <- list(
    lookfrom = c(31.5, 24.5, 38.0),
    lookat   = c(-0.2, -0.5,  0.0),
    fov      = 38.5
  )

  list(
    mesh_bed_no_fish   = mesh_bed_no_fish,
    mesh_bed_with_fish = mesh_bed_with_fish,
    mesh_water_slab    = mesh_water_slab,
    albedo_bed         = albedo_bed,
    albedo_wat_neutral = albedo_wat_neutral,
    meta_tex           = meta_tex,
    fish_specs         = fish_specs,
    water_scale        = water_scale,
    cam                = cam
  )
}

#' Render the 3D Metal GPU passes via {sonicray} and cache them
render_cutaway_gpu_layers <- function(scene, width = 1080L, height = 900L, amb_spp = 192L, light_spp = 96L,
                                      cache_dir = "data/derived/sonicray_lake_layers") {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  meta_rds <- file.path(cache_dir, sprintf("lake_diorama_gpu_%dx%d_v5.rds", width, height))
  if (file.exists(meta_rds)) {
    message("[sonicray_lake] Loading cached Metal GPU 3D diorama layers from ", meta_rds)
    return(readRDS(meta_rds))
  }

  lights <- rbind(
    data.frame(name = "sun_key", x = -22.0, y = 46.0, z =  18.0, r = 6.0, intensity = 5600),
    data.frame(name = "sun_rim", x =  28.0, y = 30.0, z = -16.0, r = 5.5, intensity = 2100)
  )
  cpp_aov <- getFromNamespace("cpp_render_ambient_and_aovs", "sonicray")

  # Pass 1: 3D Diorama + Carved Underwater Lakebed (No Fish)
  gpu_bed0 <- sr_upload_mesh(scene$mesh_bed_no_fish, scene$albedo_bed)
  res_bed0 <- sr_render_light_layers(
    gpu_scene = gpu_bed0, light_table = lights,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov = scene$cam$fov,
    width = as.integer(width), height = as.integer(height),
    amb_spp = as.integer(amb_spp), light_spp = as.integer(light_spp),
    max_bounces = 4L, clamp = 50, batch_size = 8L, out_dir = NULL
  )
  gpu_meta0 <- sr_upload_mesh(scene$mesh_bed_no_fish, scene$meta_tex)
  aov_meta0 <- cpp_aov(
    scene_xptr = gpu_meta0$gpu_handle,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov_deg = scene$cam$fov,
    width = as.integer(width), height = as.integer(height), samples = 64L, clamp_val = 4
  )

  # Pass 2: 3D Diorama + Carved Underwater Lakebed + 3D Swimming Fish (with ray-traced 3D fish shadows!)
  gpu_bedF <- sr_upload_mesh(scene$mesh_bed_with_fish, scene$albedo_bed)
  res_bedF <- sr_render_light_layers(
    gpu_scene = gpu_bedF, light_table = lights,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov = scene$cam$fov,
    width = as.integer(width), height = as.integer(height),
    amb_spp = as.integer(amb_spp), light_spp = as.integer(light_spp),
    max_bounces = 4L, clamp = 50, batch_size = 8L, out_dir = NULL
  )
  gpu_metaF <- sr_upload_mesh(scene$mesh_bed_with_fish, scene$meta_tex)
  aov_metaF <- cpp_aov(
    scene_xptr = gpu_metaF$gpu_handle,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov_deg = scene$cam$fov,
    width = as.integer(width), height = as.integer(height), samples = 64L, clamp_val = 4
  )

  # Pass 3: 3D Water Slab Envelope (`render_water` top surface + vertical cutaway water glass walls)
  gpu_wat <- sr_upload_mesh(scene$mesh_water_slab, scene$albedo_wat_neutral)
  res_wat <- sr_render_light_layers(
    gpu_scene = gpu_wat, light_table = lights,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov = scene$cam$fov,
    width = as.integer(width), height = as.integer(height),
    amb_spp = 96L, light_spp = 64L,
    max_bounces = 2L, clamp = 50, batch_size = 8L, out_dir = NULL
  )
  gpu_metaW <- sr_upload_mesh(scene$mesh_water_slab, scene$meta_tex)
  aov_metaW <- cpp_aov(
    scene_xptr = gpu_metaW$gpu_handle,
    lookfrom = scene$cam$lookfrom, lookat = scene$cam$lookat, fov_deg = scene$cam$fov,
    width = as.integer(width), height = as.integer(height), samples = 64L, clamp_val = 4
  )

  # Project the 3D fish centers & Secchi disk pole to 2D pixel space via `sr_project`
  fish_pts3d <- do.call(rbind, lapply(scene$fish_specs, function(f) c(f$x, -f$z_m * scene$water_scale, f$z)))
  fish_px <- sr_project(fish_pts3d, scene$cam$lookfrom, scene$cam$lookat, scene$cam$fov, width, height)
  fish_df <- tibble(
    x_px = fish_px[, 1],
    y_px = height - fish_px[, 2],
    z_m  = vapply(scene$fish_specs, `[[`, 0.0, "z_m")
  )

  # Project 3D vertical Secchi pole right inside the front cutaway water corner (x = 4.2, z = 12.6)
  z_pole_m <- seq(0, 8.2, by = 0.5)
  pole_pts3d <- cbind(rep(3.5, length(z_pole_m)), -z_pole_m * scene$water_scale, rep(12.6, length(z_pole_m)))
  pole_px <- sr_project(pole_pts3d, scene$cam$lookfrom, scene$cam$lookat, scene$cam$fov, width, height)
  pole_df <- tibble(z_m = z_pole_m, x_px = pole_px[, 1], y_px = height - pole_px[, 2])

  out <- list(
    width = width, height = height,
    bed0_amb = matrix(res_bed0$aovs$amb, ncol = 3, byrow = TRUE),
    bed0_sun = matrix(res_bed0$layers$sun_key, ncol = 3, byrow = TRUE),
    bed0_rim = matrix(res_bed0$layers$sun_rim, ncol = 3, byrow = TRUE),
    bed0_mask = as.numeric(res_bed0$aovs$mask),
    meta0    = matrix(aov_meta0$albedo, ncol = 3, byrow = TRUE),
    bedF_amb = matrix(res_bedF$aovs$amb, ncol = 3, byrow = TRUE),
    bedF_sun = matrix(res_bedF$layers$sun_key, ncol = 3, byrow = TRUE),
    bedF_rim = matrix(res_bedF$layers$sun_rim, ncol = 3, byrow = TRUE),
    metaF    = matrix(aov_metaF$albedo, ncol = 3, byrow = TRUE),
    wat_amb  = matrix(res_wat$aovs$amb, ncol = 3, byrow = TRUE),
    wat_sun  = matrix(res_wat$layers$sun_key, ncol = 3, byrow = TRUE),
    wat_rim  = matrix(res_wat$layers$sun_rim, ncol = 3, byrow = TRUE),
    wat_norm = matrix(res_wat$aovs$normal, ncol = 3, byrow = TRUE),
    wat_meta = matrix(aov_metaW$albedo, ncol = 3, byrow = TRUE),
    wat_mask = as.numeric(res_wat$aovs$mask),
    fish_df  = fish_df,
    pole_df  = pole_df
  )
  saveRDS(out, meta_rds)
  out
}

#' Bio-optical model mapping (Chl_a, DOC, NAP) -> kd, z_SD, ternary fractions,
#' and rayshader::render_water() style translucent water color.
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

  clarity <- exp(-0.55 * k_p)
  w_water <- exp(-1.45 * k_p)

  # Vivid rayshader::render_water() color palette:
  #   - Pure water: vibrant, deep sapphire-cobalt blue
  c_water <- c(0.04, 0.44, 0.94)
  #   - Chl a (phytoplankton): luminous emerald-jade green (clear) -> dense eutrophic algal green (bloom)
  c_chl   <- clarity * c(0.06, 0.72, 0.44) + (1 - clarity) * c(0.09, 0.42, 0.15)
  #   - DOC (humic CDOM): warm amber-tea (clear) -> dark espresso/mahogany blackwater (high DOC)
  c_doc   <- clarity * c(0.45, 0.26, 0.11) + (1 - clarity) * c(0.15, 0.07, 0.03)
  #   - NAP (sediment): turquoise-sand (clear) -> opaque riverine ochre-tan (high TSS)
  c_nap   <- clarity * c(0.56, 0.50, 0.34) + (1 - clarity) * c(0.72, 0.52, 0.26)

  col_srgb <- w_water * c_water + (1 - w_water) * (f_chl * c_chl + f_doc * c_doc + f_nap * c_nap)
  col_srgb <- pmin(pmax(col_srgb, 0.02), 0.95)
  col_lin  <- col_srgb^2.2
  hex_col  <- rgb(col_srgb[1], col_srgb[2], col_srgb[3])

  list(
    chl_a = chl_a, doc = doc, nap = nap,
    k_nap = k_nap, k_doc = k_doc, k_chl = k_chl, kw = kw, kd = kd, z_sd = z_sd,
    f_nap = f_nap, f_doc = f_doc, f_chl = f_chl,
    col_lin = col_lin, col_srgb = col_srgb, hex = hex_col,
    clarity = clarity
  )
}

#' Composite the 3D `render_water()` cutaway diorama frame using the 3 Metal GPU passes
render_lake_3d_frame <- function(gpu_layers, state) {
  W <- gpu_layers$width
  H <- gpu_layers$height

  z_bed_m <- pmax(0, gpu_layers$meta0[, 1] * 10.0)
  z_hit_m <- pmax(0, gpu_layers$metaF[, 1] * 10.0)
  is_fish <- pmin(1, pmax(0, gpu_layers$metaF[, 3]))

  # Classify the 3D water slab regions from `wat_meta[, 2]`:
  #   wat_m2 > 0.85 -> top water surface (`top_water`)
  #   wat_m2 in [0.45, 0.65] -> front vertical glass cutaway water wall (`cut_front`)
  #   wat_m2 in [0.65, 0.85] -> right vertical glass cutaway water wall (`cut_right`)
  wat_m2 <- gpu_layers$wat_meta[, 2]
  top_water <- pmin(1, pmax(0, (wat_m2 - 0.84) / 0.14)) * gpu_layers$wat_mask
  cut_front <- pmin(1, pmax(0, (wat_m2 > 0.44 & wat_m2 <= 0.64) * 1.0)) * (1 - top_water) * gpu_layers$wat_mask
  cut_right <- pmin(1, pmax(0, (wat_m2 > 0.64 & wat_m2 <= 0.84) * 1.0)) * (1 - top_water) * gpu_layers$wat_mask
  cut_water <- pmin(1, cut_front + cut_right)
  any_water <- pmin(1, top_water + cut_water)

  # Direct sun attenuates with depth z_bed_m (softening & removing fish shadows when water becomes murky!)
  sun_atten_bed <- ifelse(any_water > 0.1, exp(-0.35 * state$kd * z_bed_m), 1.0)
  sun_atten_hit <- ifelse(any_water > 0.1, exp(-0.35 * state$kd * z_hit_m), 1.0)

  L_bed0 <- sweep(gpu_layers$bed0_sun, 2, c(1.00, 0.96, 0.88) * 0.035, `*`) * sun_atten_bed +
            sweep(gpu_layers$bed0_rim, 2, c(0.70, 0.84, 1.00) * 0.019, `*`) +
            sweep(gpu_layers$bed0_amb, 2, c(0.62, 0.78, 0.98) * 0.30,  `*`)

  L_bedF <- sweep(gpu_layers$bedF_sun, 2, c(1.00, 0.96, 0.88) * 0.035, `*`) * sun_atten_hit +
            sweep(gpu_layers$bedF_rim, 2, c(0.70, 0.84, 1.00) * 0.019, `*`) +
            sweep(gpu_layers$bedF_amb, 2, c(0.62, 0.78, 0.98) * 0.30,  `*`)

  tm <- function(x) {
    m <- x / (1 + 0.50 * x)
    pmin(pmax(m^(1 / 2.2), 0), 1)
  }
  srgb_bed0 <- tm(L_bed0)
  srgb_bedF <- tm(L_bedF)

  # Volumetric water color (shaded by uniform-albedo path-traced mountain shadows on the lake surface!)
  surf_illum <- 0.78 + 0.22 * pmin(1.25, (gpu_layers$wat_sun[, 1] * 0.045 + gpu_layers$wat_amb[, 1] * 0.45))
  # Side-wall bevel shading on the front (`cut_front`) vs right (`cut_right`) 3D cutaway glass faces
  wall_bevel <- 1.0 - 0.08 * cut_front - 0.22 * cut_right - 0.12 * cut_water * gpu_layers$wat_meta[, 3]
  vol_col <- cbind(
    state$col_srgb[1] * surf_illum * wall_bevel,
    state$col_srgb[2] * surf_illum * wall_bevel,
    state$col_srgb[3] * surf_illum * wall_bevel
  )

  #Spectral water absorption on the underwater lakebed (`rayshader::render_water` dielectric look!):
  z_eff <- pmax(0.25, ifelse(cut_water > 0.1, pmax(z_bed_m, gpu_layers$wat_meta[, 3] * 6.5), z_bed_m))
  abs_w <- pmin(0.88, 1 - exp(-0.34 * z_eff))
  bed0_tinted <- cbind(
    srgb_bed0[, 1] * (1 - abs_w) + srgb_bed0[, 1] * state$col_srgb[1] * 1.15 * abs_w,
    srgb_bed0[, 2] * (1 - abs_w) + srgb_bed0[, 2] * state$col_srgb[2] * 1.25 * abs_w,
    srgb_bed0[, 3] * (1 - abs_w) + pmax(srgb_bed0[, 3], 0.44) * state$col_srgb[3] * 1.38 * abs_w
  )

  # Apply ray-traced 3D fish shadow onto the lakebed (shadow contrast decays with depth & kd)
  vis_shadow <- exp(-0.85 * state$kd * z_bed_m)
  shadow_ratio <- ifelse(is_fish > 0.15, 1.0, pmin(1.0, pmax(0.25, (srgb_bedF[, 2] + 0.02) / (srgb_bed0[, 2] + 0.02))))
  bed_with_shadows <- bed0_tinted * (1.0 - vis_shadow * (1.0 - shadow_ratio))

  # Lakebed + shadow composite seen through water volume of depth `z_eff`
  T_bed <- exp(-0.30 * state$kd * z_eff) * (0.84 - 0.10 * cut_water)
  water_bg_srgb <- bed_with_shadows * T_bed + vol_col * (1 - T_bed)

  # Now composite the 3D path-traced fish at their ACTUAL depth `z_hit_m` (`0.8 m`, `2.2 m`, `4.2 m`, `6.8 m`)!
  # Because a fish at `z_hit_m` only has `z_hit_m` of water above it, shallow fish shine brightly with natural
  # silver & coral-red colors while deeper fish take on more sapphire water tint and disappear first as kd rises!
  T_fish <- exp(-1.15 * state$kd * z_hit_m)
  fish_water_tint <- exp(-0.18 * z_hit_m)
  fish_srgb <- cbind(
    srgb_bedF[, 1] * (0.35 + 0.65 * fish_water_tint) + state$col_srgb[1] * (1 - fish_water_tint) * 0.45,
    srgb_bedF[, 2] * (0.45 + 0.55 * fish_water_tint) + state$col_srgb[2] * (1 - fish_water_tint) * 0.55,
    srgb_bedF[, 3] * (0.55 + 0.45 * fish_water_tint) + state$col_srgb[3] * (1 - fish_water_tint) * 0.65
  )
  fish_w <- is_fish * T_fish
  water_srgb <- fish_srgb * fish_w + water_bg_srgb * (1 - fish_w)

  # Realistic Dielectric Fresnel Reflection & Subtle Wave Glints on Top Water + Cutaway Glass Bevel Edges
  wn_x <- gpu_layers$wat_norm[, 1]
  wn_y <- gpu_layers$wat_norm[, 2]
  wn_z <- gpu_layers$wat_norm[, 3]
  cos_theta <- pmax(0.08, 0.50 * wn_x + 0.46 * abs(wn_y) + 0.58 * wn_z)
  fresnel_top <- top_water * (0.02 + 0.13 * (1 - cos_theta)^3.5)

  # Crisp `rayshader::render_water()` glass rim highlights along:
  #   1. The top water-air cutaway meniscus (where `top_water` meets `cut_water`)
  #   2. The vertical front-right glass corner (where `cut_front` meets `cut_right`)
  #   3. The bottom water-lakebed cutaway interface (`wat_meta[, 3] > 0.90`)
  rim_meniscus <- ((top_water > 0.08 & top_water < 0.92) & (wat_m2 > 0.40)) * 0.36 +
                  (cut_front > 0.05 & cut_right > 0.05) * 0.42 +
                  (cut_water > 0.2 & gpu_layers$wat_meta[, 3] > 0.91) * 0.26 +
                  cut_front * 0.06 + cut_right * 0.03

  water_srgb[, 1] <- pmin(1, water_srgb[, 1] + fresnel_top * 0.55 + rim_meniscus * 0.58)
  water_srgb[, 2] <- pmin(1, water_srgb[, 2] + fresnel_top * 0.78 + rim_meniscus * 0.82)
  water_srgb[, 3] <- pmin(1, water_srgb[, 3] + fresnel_top * 1.00 + rim_meniscus * 1.00)

  # Clean dark studio background around the floating 3D `render_water()` diorama block
  yy <- rep(seq(0, 1, length.out = H), each = W)
  bg_srgb <- cbind(
    0.043 + 0.020 * yy,
    0.067 + 0.030 * yy,
    0.118 + 0.045 * yy
  )

  mask_all <- pmax(gpu_layers$bed0_mask, gpu_layers$wat_mask)
  diorama_srgb <- srgb_bed0 * (1 - any_water) + water_srgb * any_water
  final_srgb <- diorama_srgb * mask_all + bg_srgb * (1 - mask_all)

  out <- array(0, c(H, W, 3))
  out[, , 1] <- matrix(final_srgb[, 1], nrow = H, ncol = W, byrow = TRUE)
  out[, , 2] <- matrix(final_srgb[, 2], nrow = H, ncol = W, byrow = TRUE)
  out[, , 3] <- matrix(final_srgb[, 3], nrow = H, ncol = W, byrow = TRUE)
  out
}

#' Define the 28-second optical titration schedule
titration_schedule <- function(t_sec, total_s = 28.0) {
  chl0 <- 0.30; doc0 <- 0.35; nap0 <- 0.35

  if (t_sec < 3.2) {
    phase <- "1. Clear Oligotrophic Lake"
    chl <- chl0; doc <- doc0; nap <- nap0
    mode <- "clear"
  } else if (t_sec < 9.2) {
    u <- (t_sec - 3.2) / 6.0
    pulse <- sin(pi * smooth01(u))^1.2
    chl <- chl0 + 130.0 * pulse
    doc <- doc0 + 0.25 * pulse
    nap <- nap0 + 0.30 * pulse
    phase <- "2. Titrating Chlorophyll a (Algae)"
    mode <- if (pulse > 0.32) "dark" else "clear"
  } else if (t_sec < 15.2) {
    u <- (t_sec - 9.2) / 6.0
    pulse <- sin(pi * smooth01(u))^1.2
    chl <- chl0 + 0.4 * pulse
    doc <- doc0 + 28.0 * pulse
    nap <- nap0 + 0.25 * pulse
    phase <- "3. Titrating DOC (Brownification)"
    mode <- if (pulse > 0.32) "dark" else "clear"
  } else if (t_sec < 21.2) {
    u <- (t_sec - 15.2) / 6.0
    pulse <- sin(pi * smooth01(u))^1.2
    chl <- chl0 + 0.4 * pulse
    doc <- doc0 + 0.25 * pulse
    nap <- nap0 + 75.0 * pulse
    phase <- "4. Titrating TSS (Sediment Plume)"
    mode <- if (pulse > 0.32) "dark" else "clear"
  } else {
    u <- (t_sec - 21.2) / 6.8
    if (u < 0.46) {
      w <- smooth01(u / 0.46)
      chl <- 1.6 + 0.8 * sin(2 * pi * w)
      doc <- 0.22 + 0.10 * cos(2 * pi * w)
      nap <- 0.36 + 0.16 * sin(2 * pi * w + 1.0)
      phase <- "5a. Balanced Mix — CLEAR State"
      mode <- "clear"
    } else {
      w <- smooth01((u - 0.46) / 0.54)
      mult <- 1.0 + 32.0 * sin(pi * pmin(1, w * 1.15))^1.1
      chl <- 1.7 * mult
      doc <- 0.23 * mult
      nap <- 0.38 * mult
      phase <- "5b. Same Mix — DARK / TURBID State"
      mode <- "dark"
    }
  }

  st <- bio_optical_state(chl, doc, nap)
  st$t_sec <- t_sec
  st$phase <- phase
  st$mode <- mode
  st
}

tern_xy_pt <- function(f_nap, f_doc, f_chl) {
  s <- f_nap + f_doc + f_chl
  f_nap <- f_nap / s; f_doc <- f_doc / s; f_chl <- f_chl / s
  x <- f_doc + 0.5 * f_chl
  y <- (sqrt(3) / 2) * f_chl
  cbind(x = x, y = y)
}

build_ternary_state_grids <- function(res = 26) {
  pts <- list()
  for (i in 0:res) {
    for (j in 0:(res - i)) {
      k <- res - i - j
      f_nap <- max(0.01, i / res)
      f_doc <- max(0.01, j / res)
      f_chl <- max(0.01, k / res)
      s <- f_nap + f_doc + f_chl
      f_nap <- f_nap / s; f_doc <- f_doc / s; f_chl <- f_chl / s

      kp_clear <- 0.065
      st_clear <- bio_optical_state(
        chl_a = (f_chl * kp_clear) / 0.0086,
        doc   = (f_doc * kp_clear) / 0.0732,
        nap   = (f_nap * kp_clear) / 0.0432
      )

      kp_dark <- 2.75
      st_dark <- bio_optical_state(
        chl_a = (f_chl * kp_dark) / 0.0086,
        doc   = (f_doc * kp_dark) / 0.0732,
        nap   = (f_nap * kp_dark) / 0.0432
      )

      xy <- tern_xy_pt(f_nap, f_doc, f_chl)
      pts[[length(pts) + 1]] <- tibble(
        f_nap = f_nap, f_doc = f_doc, f_chl = f_chl,
        x = xy[1, "x"], y = xy[1, "y"],
        col_clear = st_clear$hex, z_sd_clear = st_clear$z_sd,
        col_dark  = st_dark$hex,  z_sd_dark  = st_dark$z_sd
      )
    }
  }
  bind_rows(pts)
}

#' Draw a single minimalist, visual-first widescreen frame (1600 x 900)
draw_composite_frame <- function(out_png, lake_img, gpu_layers, st, trail_df, tgrid, width = 1600, height = 900) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(out_png, width = width, height = height, units = "px", res = 144, background = "#0b111e")
  } else {
    png(out_png, width = width, height = height, units = "px", res = 144, bg = "#0b111e", type = "quartz")
  }
  on.exit(dev.off(), add = TRUE)

  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", bg = "#0b111e", fg = "white")
  plot.new()
  plot.window(xlim = c(0, 1600), ylim = c(0, 900))

  # 1. Left Panel (x: 0 -> 1080, y: 0 -> 900): 3D `render_water()` Diorama + Path-Traced 3D Fish
  rasterImage(lake_img, xleft = 0, ybottom = 0, xright = 1080, ytop = 900, interpolate = TRUE)

  # Subtle depth badges next to the 4 primary 3D path-traced fish (only shown when fish is visible!)
  fdf <- gpu_layers$fish_df[c(1, 3, 5, 7), ]
  for (i in seq_len(nrow(fdf))) {
    vis <- exp(-1.15 * st$kd * fdf$z_m[i])
    if (vis > 0.22) {
      text(fdf$x_px[i], fdf$y_px[i] + 14, sprintf("%.1f m", fdf$z_m[i]),
           col = rgb(1, 1, 1, pmin(1, vis * 1.15)), font = 2, cex = 0.56)
    }
  }

  # 3D-Projected Secchi Disk & Depth Scale along the front cutaway edge
  pdf <- gpu_layers$pole_df
  segments(pdf$x_px[1], pdf$y_px[1], pdf$x_px[nrow(pdf)], pdf$y_px[nrow(pdf)],
           col = rgb(1, 1, 1, 0.55), lwd = 1.4)
  for (zd in c(0, 2, 4, 6, 8)) {
    row_z <- pdf[which.min(abs(pdf$z_m - zd)), ]
    segments(row_z$x_px - 6, row_z$y_px, row_z$x_px, row_z$y_px, col = rgb(1, 1, 1, 0.78), lwd = 1.3)
    text(row_z$x_px - 10, row_z$y_px, sprintf("%dm", zd), adj = c(1, 0.5), col = "white", font = 2, cex = 0.56)
  }
  z_disk <- pmin(st$z_sd, 8.2)
  x_disk <- approx(pdf$z_m, pdf$x_px, xout = z_disk)$y + 22
  y_disk <- approx(pdf$z_m, pdf$y_px, xout = z_disk)$y
  y_top  <- pdf$y_px[1]
  segments(x_disk, y_top + 10, x_disk, y_disk, col = rgb(1, 1, 1, 0.90), lwd = 1.5)
  ang <- seq(0, 2 * pi, length.out = 48)
  rx <- 20; ry <- 8
  polygon(x_disk + rx * cos(ang), y_disk + ry * sin(ang), col = "white", border = "black", lwd = 1.1)
  a1 <- seq(0, pi / 2, length.out = 14); a2 <- seq(pi, 3 * pi / 2, length.out = 14)
  polygon(c(x_disk, x_disk + rx * cos(a1), x_disk), c(y_disk, y_disk + ry * sin(a1), y_disk), col = "#0f172a", border = NA)
  polygon(c(x_disk, x_disk + rx * cos(a2), x_disk), c(y_disk, y_disk + ry * sin(a2), y_disk), col = "#0f172a", border = NA)
  polygon(x_disk + rx * cos(ang), y_disk + ry * sin(ang),
          col = rgb(st$col_srgb[1], st$col_srgb[2], st$col_srgb[3], 0.35), border = "#fde047", lwd = 1.5)

  # Minimal, single-line header bar at top of 3D panel
  rect(24, 832, 1056, 882, col = rgb(0.05, 0.08, 0.14, 0.78), border = rgb(1, 1, 1, 0.18), lwd = 1.2)
  text(42, 857, st$phase, adj = c(0, 0.5), col = "white", font = 2, cex = 1.02)
  rect(856, 838, 1044, 876, col = st$hex, border = "white", lwd = 1.6)
  tcol_top <- if (mean(st$col_srgb) > 0.52) "#0f172a" else "#ffffff"
  text(950, 857, sprintf("Secchi = %.1f m", st$z_sd), col = tcol_top, font = 2, cex = 0.84)

  # ===========================================================================
  # 2. Right Panel (x: 1080 -> 1600, y: 0 -> 900) — Minimalist Visual HUD
  # ===========================================================================
  rect(1080, 0, 1600, 900, col = "#0b111e", border = "#334155", lwd = 1.5)

  # Top: 3 OAC Bars (y: 660 -> 875)
  text(1104, 860, "Constituent Titration", adj = c(0, 0.5), col = "white", font = 2, cex = 0.90)

  draw_clean_bar <- function(yc, name, val_txt, pct_txt, frac, col_bar) {
    text(1104, yc + 11, name, adj = c(0, 0.5), col = "#e2e8f0", font = 2, cex = 0.72)
    text(1576, yc + 11, paste0(val_txt, "  (", pct_txt, ")"), adj = c(1, 0.5), col = "#f8fafc", font = 2, cex = 0.70)
    rect(1104, yc - 11, 1576, yc - 1, col = "#1e293b", border = "#475569", lwd = 1)
    w <- pmin(1, pmax(0.015, frac)) * (1576 - 1104)
    rect(1104, yc - 11, 1104 + w, yc - 1, col = col_bar, border = NA)
  }
  draw_clean_bar(810, "Chl a",     sprintf("%.1f ug/L", st$chl_a), sprintf("%.0f%%", 100 * st$f_chl), (st$chl_a / 130)^0.62, "#10b981")
  draw_clean_bar(752, "DOC",       sprintf("%.1f mg/L", st$doc),   sprintf("%.0f%%", 100 * st$f_doc), (st$doc   / 28)^0.62,  "#b45309")
  draw_clean_bar(694, "TSS / NAP", sprintf("%.1f mg/L", st$nap),   sprintf("%.0f%%", 100 * st$f_nap), (st$nap   / 75)^0.62,  "#eab308")

  # 100% Partition Bar
  x0 <- 1104; xw <- 472
  w1 <- st$f_chl * xw; w2 <- st$f_doc * xw
  rect(x0, 636, x0 + w1, 656, col = "#10b981", border = NA)
  rect(x0 + w1, 636, x0 + w1 + w2, 656, col = "#78350f", border = NA)
  rect(x0 + w1 + w2, 636, x0 + xw, 656, col = "#eab308", border = NA)
  rect(x0, 636, x0 + xw, 656, col = NA, border = "#cbd5e1", lwd = 1.1)

  # Middle: Clear State vs. Dark State Ternary Triangles
  segments(1100, 612, 1580, 612, col = "#334155", lwd = 1.2)
  text(1340, 586, "Ternary Composition vs. Total Visibility", col = "#fde047", font = 2, cex = 0.82)

  draw_clean_ternary <- function(cx_left, cx_right, cy_bot, cy_top, title_txt, state_type = c("clear", "dark"), is_active = FALSE) {
    state_type <- match.arg(state_type)
    b_col <- if (is_active) "#fde047" else "#334155"
    b_lwd <- if (is_active) 2.4 else 1.2
    rect(cx_left - 8, cy_bot - 52, cx_right + 8, cy_top + 42, col = "#131c2e", border = b_col, lwd = b_lwd)

    cx_mid <- 0.5 * (cx_left + cx_right)
    text(cx_mid, cy_top + 22, title_txt, col = if (is_active) "#fde047" else "white", font = 2, cex = 0.78)

    tw <- cx_right - cx_left
    th <- tw * (sqrt(3) / 2)
    tx <- function(x) cx_left + x * tw
    ty <- function(y) cy_bot + (y / (sqrt(3) / 2)) * th

    reveal_p <- smoothstep_sr(1.5, 20.0, st$t_sec)
    cols_vec <- if (state_type == "clear") tgrid$col_clear else tgrid$col_dark
    cur_xy <- tern_xy_pt(st$f_nap, st$f_doc, st$f_chl)
    d_cur <- sqrt((tgrid$x - cur_xy[1, "x"])^2 + (tgrid$y - cur_xy[1, "y"])^2)
    pt_active <- (d_cur <= (0.22 + 0.85 * reveal_p))

    points(tx(tgrid$x), ty(tgrid$y), pch = 16, cex = ifelse(pt_active, 1.10, 0.75),
           col = ifelse(pt_active, cols_vec, "#1e293b"))

    s32 <- sqrt(3) / 2
    polygon(c(tx(0), tx(1), tx(0.5)), c(ty(0), ty(0), ty(s32)), col = NA, border = "white", lwd = 1.5)

    text(tx(0.5), ty(s32) - 11, "Chl a", col = "#4ade80", font = 2, cex = 0.64)
    text(tx(0.0), ty(0) - 12, "TSS", adj = c(0, 0.5), col = "#fde047", font = 2, cex = 0.64)
    text(tx(1.0), ty(0) - 12, "DOC", adj = c(1, 0.5), col = "#f59e0b", font = 2, cex = 0.64)

    if (nrow(trail_df) > 1) {
      lines(tx(trail_df$x), ty(trail_df$y), col = rgb(1, 1, 1, 0.58), lwd = 1.5)
    }
    points(tx(cur_xy[1, "x"]), ty(cur_xy[1, "y"]), pch = 21, bg = st$hex, col = "#fde047", cex = 2.1, lwd = 2.0)

    kp_ref <- if (state_type == "clear") 0.065 else 2.75
    st_ref <- bio_optical_state(
      chl_a = (st$f_chl * kp_ref) / 0.0086,
      doc   = (st$f_doc * kp_ref) / 0.0732,
      nap   = (st$f_nap * kp_ref) / 0.0432
    )
    rect(cx_left + 6, cy_bot - 44, cx_right - 6, cy_bot - 20, col = st_ref$hex, border = "white", lwd = 1.1)
    tcol <- if (mean(st_ref$col_srgb) > 0.52) "#0f172a" else "#ffffff"
    text(cx_mid, cy_bot - 32, sprintf("Secchi = %.1f m", st_ref$z_sd), col = tcol, font = 2, cex = 0.68)
  }

  draw_clean_ternary(
    cx_left = 1102, cx_right = 1318, cy_bot = 265, cy_top = 505,
    title_txt = "CLEAR STATE",
    state_type = "clear",
    is_active = (st$mode == "clear")
  )

  draw_clean_ternary(
    cx_left = 1362, cx_right = 1578, cy_bot = 265, cy_top = 505,
    title_txt = "DARK STATE",
    state_type = "dark",
    is_active = (st$mode == "dark")
  )

  # Bottom: Visual Underwater Fish Visibility Meter
  rect(1094, 28, 1586, 185, col = "#131c2e", border = "#334155", lwd = 1.2)
  text(1340, 162, "Underwater Fish Visibility by Depth", col = "white", font = 2, cex = 0.76)

  f_depths <- c(0.8, 2.2, 4.2, 6.8)
  x_centers <- seq(1156, 1524, length.out = 4)
  for (k in seq_along(f_depths)) {
    zm <- f_depths[k]
    xc <- x_centers[k]
    vis <- exp(-1.15 * st$kd * zm)
    box_col <- if (vis > 0.25) "#1e293b" else "#0f172a"
    brd_col <- if (vis > 0.25) "#38bdf8" else "#334155"
    rect(xc - 50, 48, xc + 50, 138, col = box_col, border = brd_col, lwd = if (vis > 0.25) 1.6 else 1.0)
    t_seq <- seq(0, 2 * pi, length.out = 32)
    bx <- xc + 22 * cos(t_seq)
    by <- 102 + 7 * sin(t_seq) * (0.75 + 0.35 * cos(t_seq))
    tx_p <- xc - 22 * c(0.75, 1.25, 1.05, 1.25)
    ty_p <- 102 + 7 * c(0, 0.9, 0, -0.9)
    f_alpha <- pmax(0.08, pmin(1.0, vis))
    polygon(tx_p, ty_p, col = rgb(0.95, 0.45, 0.35, f_alpha), border = NA)
    polygon(bx, by, col = rgb(0.92, 0.82, 0.72, f_alpha), border = rgb(1, 1, 1, f_alpha * 0.6), lwd = 0.8)
    text(xc, 74, sprintf("%.1f m", zm), col = "white", font = 2, cex = 0.68)
    status_txt <- if (vis > 0.50) "Clear" else if (vis > 0.20) "Faint" else "Hidden"
    status_col <- if (vis > 0.50) "#4ade80" else if (vis > 0.20) "#fde047" else "#64748b"
    text(xc, 57, status_txt, col = status_col, font = 2, cex = 0.60)
  }
}

#' Render the 6-stage contact sheet (`figs/sonicray_optical_titration_sheet.png`)
#' and optionally the full 28-second MP4 (`figs/sonicray_optical_titration.mp4`).
build_sonicray_titration_visuals <- function(
  out_mp4   = "figs/sonicray_optical_titration.mp4",
  out_sheet = "figs/sonicray_optical_titration_sheet.png",
  fps       = 15L,
  duration_s = 28.0,
  sheet_only = FALSE
) {
  t0 <- Sys.time()
  scene <- build_lake_cutaway_scene(seed = 42L, n = 400L, tex = 2048L)
  gpu_layers <- render_cutaway_gpu_layers(scene, width = 1080L, height = 900L, amb_spp = 192L, light_spp = 96L)
  tgrid <- build_ternary_state_grids(res = 26)

  times <- seq(0, duration_s, by = 1 / fps)
  states <- lapply(times, titration_schedule, total_s = duration_s)
  xy_all <- do.call(rbind, lapply(states, function(s) tern_xy_pt(s$f_nap, s$f_doc, s$f_chl)))
  trail_all <- tibble(x = xy_all[, "x"], y = xy_all[, "y"])

  key_times <- c(1.5, 6.2, 12.2, 18.2, 22.8, 26.2)
  key_idx <- vapply(key_times, function(kt) which.min(abs(times - kt)), 1L)

  frames_dir <- file.path(tempdir(), "sonicray_titration_frames_v4")
  unlink(frames_dir, recursive = TRUE)
  dir.create(frames_dir, recursive = TRUE, showWarnings = FALSE)

  render_indices <- if (sheet_only) key_idx else seq_along(times)
  message(sprintf("[sonicray_lake] Rendering %d frames...", length(render_indices)))
  for (i in render_indices) {
    f_png <- file.path(frames_dir, sprintf("frame_%04d.png", i))
    st_i <- states[[i]]
    if (sheet_only && i == key_idx[5]) st_i$t_sec <- 22.8
    if (sheet_only && i == key_idx[6]) st_i$t_sec <- 26.2
    img3d <- render_lake_3d_frame(gpu_layers, st_i)
    draw_composite_frame(f_png, img3d, gpu_layers, st_i, trail_all[1:i, , drop = FALSE], tgrid)
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
