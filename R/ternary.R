# -----------------------------------------------------------------------------
# ternary.R -- lightweight ternary diagrams and exact simplex polygon clipping
#
# Orientation (same as the paper):
#   top vertex = DOC (100%), bottom-left = Chl a (100%), bottom-right = NAP (100%)
#   Cartesian: x = NAP + DOC/2,  y = DOC * sqrt(3)/2   (fractions sum to 1)
# -----------------------------------------------------------------------------

S3 <- sqrt(3) / 2

#' Convert compositions to Cartesian coordinates (adds .x, .y).
tern_xy <- function(df, chl = chl_kd, doc = doc_kd, nap = tss_kd) {
  df %>% mutate(.tot = {{ chl }} + {{ doc }} + {{ nap }},
                .x = ({{ nap }} + {{ doc }} / 2) / .tot,
                .y = ({{ doc }} * S3) / .tot) %>%
    select(-.tot)
}
bary_to_xy <- function(chl, doc, nap) cbind(x = nap + doc / 2, y = doc * S3)

#' Empty ternary frame: triangle, dashed 20% grid, tick labels, axis titles.
tern_frame <- function(base_size = 9, grid = TRUE, labels = TRUE, titles = TRUE,
                       lab_chl = 'Chl a', lab_doc = 'DOC', lab_nap = 'NAP',
                       grid_col = 'grey80', line_col = 'grey15', compact = FALSE) {
  v <- c(0.2, 0.4, 0.6, 0.8)
  grid_df <- bind_rows(
    tibble(id = paste0('n', v), x = v,                 y = 0,          xend = v + (1 - v) / 2,  yend = (1 - v) * S3),
    tibble(id = paste0('d', v), x = v / 2,             y = v * S3,     xend = 1 - v / 2,        yend = v * S3),
    tibble(id = paste0('c', v), x = 1 - v,             y = 0,          xend = (1 - v) / 2,      yend = (1 - v) * S3))
  tri <- tibble(x = c(0, 1, 0.5, 0), y = c(0, 0, S3, 0))
  tk <- c(20, 40, 60, 80, 100); f <- tk / 100
  fs <- base_size / ggplot2::.pt
  pm <- if (compact) margin(0, 0, 0, 0) else margin(4, 14, 4, 14)
  p <- ggplot() + coord_equal(clip = 'off') + theme_void(base_size = base_size) +
    theme(plot.margin = pm, legend.position = 'none')
  if (grid) p <- p + geom_segment(data = grid_df, aes(x, y, xend = xend, yend = yend),
                                  colour = grid_col, linewidth = 0.3, linetype = 'dashed')
  p <- p + geom_path(data = tri, aes(x, y), colour = line_col, linewidth = 0.65)
  if (labels) {
    lab_df <- bind_rows(
      # NAP along bottom edge (left -> right)
      tibble(x = f + 0.015, y = -0.032, label = tk, angle = 60, hjust = 1),
      # DOC along right edge (bottom -> top)
      tibble(x = 1 - f / 2 + 0.026, y = f * S3, label = tk, angle = 0, hjust = 0),
      # Chl a along left edge (top -> bottom-left)
      tibble(x = (1 - f) / 2 - 0.026, y = (1 - f) * S3 + 0.01, label = tk, angle = -60, hjust = 1))
    p <- p + geom_text(data = lab_df, aes(x, y, label = label, angle = angle, hjust = hjust),
                       size = fs * 0.9, colour = 'grey25', inherit.aes = FALSE)
  }
  if (titles) {
    tt_df <- if (compact) {
      tibble(x = c(0.5, 0.835, 0.165), y = c(-0.145, 0.505, 0.505), angle = c(0, -60, 60),
             label = c(paste(lab_nap, '\u2192'), paste('\u2190', lab_doc), paste('\u2190', lab_chl)))
    } else {
      tibble(x = c(0.5, 0.86, 0.14), y = c(-0.17, 0.52, 0.52), angle = c(0, -60, 60),
             label = c(paste(lab_nap, '\u2192'), paste('\u2190', lab_doc), paste('\u2190', lab_chl)))
    }
    p <- p + geom_text(data = tt_df, aes(x, y, label = label, angle = angle),
                       size = fs * 1.05, inherit.aes = FALSE)
  }
  xlims <- if (compact) c(-0.10, 1.12) else c(-0.2, 1.2)
  ylims <- if (compact) c(-0.17, 0.90) else c(-0.22, 0.93)
  p + scale_x_continuous(limits = xlims, expand = c(0, 0)) +
    scale_y_continuous(limits = ylims, expand = c(0, 0))
}

# ---- Mean +/- k SD region ("hexagon") ----------------------------------------

#' Clip a convex polygon (rows = vertices, cols = chl, doc, nap) to the
#' half-space  sign * (p[, j] - bound) >= 0  (Sutherland-Hodgman).
clip_halfspace <- function(poly, j, bound, sign) {
  if (nrow(poly) == 0 || is.na(bound)) return(poly)
  inside <- function(p) sign * (p[j] - bound) >= -1e-12
  out <- matrix(numeric(0), ncol = 3)
  n <- nrow(poly)
  for (i in seq_len(n)) {
    cur <- poly[i, ]; prv <- poly[if (i == 1) n else i - 1, ]
    if (inside(cur)) {
      if (!inside(prv)) { t <- (bound - prv[j]) / (cur[j] - prv[j]); out <- rbind(out, prv + t * (cur - prv)) }
      out <- rbind(out, cur)
    } else if (inside(prv)) {
      t <- (bound - prv[j]) / (cur[j] - prv[j]); out <- rbind(out, prv + t * (cur - prv))
    }
  }
  out
}

#' Region of the simplex where every component lies within mean +/- k*SD
#' (bounds truncated to [0, 1]). Returns vertices in counter-clockwise polygon order.
tern_sd_region <- function(df, k = 1, chl = chl_kd, doc = doc_kd, nap = tss_kd) {
  M <- cbind(pull(df, {{ chl }}), pull(df, {{ doc }}), pull(df, {{ nap }}))
  M <- M[complete.cases(M), , drop = FALSE]
  if (nrow(M) < 2) return(tibble(chl = numeric(0), doc = numeric(0), nap = numeric(0), .x = numeric(0), .y = numeric(0)))
  mu <- colMeans(M); s <- apply(M, 2, sd)
  lo <- pmax(0, mu - k * s); hi <- pmin(1, mu + k * s)
  poly <- diag(3)      # simplex vertices: (1,0,0), (0,1,0), (0,0,1)
  for (j in 1:3) {
    poly <- clip_halfspace(poly, j, lo[j], +1)
    poly <- clip_halfspace(poly, j, hi[j], -1)
  }
  colnames(poly) <- c('chl', 'doc', 'nap')
  as_tibble(poly) %>% mutate(.x = nap + doc / 2, .y = doc * S3)
}

#' Exact area of the mean +/- k*SD simplex region as a fraction of the whole
#' unit-side equilateral triangle (area = sqrt(3) / 4 = 0.4330127).
tern_area_fraction <- function(df, k = 2, ...) {
  r <- tern_sd_region(df, k = k, ...)
  if (nrow(r) < 3) return(0)
  a <- abs(sum(r$.x * c(r$.y[-1], r$.y[1]) - c(r$.x[-1], r$.x[1]) * r$.y)) / 2
  a / (sqrt(3) / 4)
}
