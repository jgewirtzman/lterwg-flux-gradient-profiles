# diagnose_wind_direction_value.R
# How much identifying variation does the symmetric-wind assumption destroy?
#
# Two sources of variation in f_wet(z,t):
#   (a) footprint SIZE   — varies with stability/u*; survives rotational averaging
#   (b) footprint AIMING — which way the footprint points; averaged away entirely
#
# (b) is the useful one. Size is collinear with dArea and with the met variables
# that independently drive the flux, so a regression of dConc on Δf_wet that has
# only (a) cannot separate "the footprint got wetter" from "the footprint got
# bigger". Aiming is independent of size, so it identifies composition.
#
# Both g(r,theta) and w(r,theta) live on the same polar grid, so rotating the
# footprint is a circular shift of the angular index — we can evaluate every
# wind direction exactly, with no re-intersection.
#
# Usage:  Rscript diagnose_wind_direction_value.R [SITE ...]

source("config.R")

sites <- { a <- commandArgs(trailingOnly = TRUE)
           if (length(a)) a else c("UNDE", "HARV", "TREE", "DELA") }
n_ts  <- 30
n_r   <- 90
n_ang <- 180
R_CAP <- 5000
nwi_dir <- file.path(proj_dir, "data", "nwi")

set.seed(7)

# Mask on the polar grid (n_r x n_ang) rather than the radial average.
polar_mask <- function(geom, grid, target_crs) {
  pts <- st_sfc(lapply(seq_along(grid$x),
                       function(i) st_point(c(grid$x[i], grid$y[i]))),
                crs = target_crs)
  hit <- lengths(suppressMessages(st_intersects(pts, geom))) > 0
  matrix(hit, nrow = grid$n_r, ncol = grid$n_ang)
}

out <- list()
for (site in sites) {
  message("\n=== ", site, " ===")
  nwi_proj <- load_nwi_for_site(site, nwi_dir)
  if (is.null(nwi_proj)) { message("  no NWI"); next }
  sd_ <- tryCatch(load_all_sites(site, data_dir), error = function(e) NULL)
  if (is.null(sd_) || is.null(sd_[[site]]$aligned$CH4)) { message("  no data"); next }

  attr_df <- sd_[[site]]$attr
  lat <- as.numeric(attr_df$LatTow[1]); lon <- as.numeric(attr_df$LonTow[1])
  target_crs <- pick_local_aea(lon, lat)
  nwi_proj <- suppressWarnings(st_make_valid(st_transform(nwi_proj, target_crs)))

  hl <- attr_df %>% select(TowerPosition, DistZaxsLvlMeasTow) %>% distinct() %>%
    mutate(TowerPosition = as.numeric(TowerPosition),
           z_m = as.numeric(DistZaxsLvlMeasTow)) %>%
    filter(!is.na(z_m), z_m > 0) %>% arrange(z_m)

  df <- sd_[[site]]$aligned$CH4
  df <- df[which(is.finite(as.numeric(df$ustar_interp)) &
                 is.finite(as.numeric(df$L_obukhov)) &
                 is.finite(as.numeric(df$roughLength_calc))), , drop = FALSE]
  df <- df[sample.int(nrow(df), min(n_ts, nrow(df))), , drop = FALSE]

  contours <- list()
  for (i in seq_len(nrow(df))) {
    fps <- tryCatch(compute_footprints_all_heights(
      met_row = df[i, , drop = FALSE], heights = hl$z_m,
      tower_positions = hl$TowerPosition, lat = lat, r = 80),
      error = function(e) NULL)
    if (is.null(fps)) next
    for (k in seq_len(nrow(fps))) {
      f <- fps$footprint[[k]]
      if (isTRUE(f$flag_err == 1) || length(f$xr) < 3) next
      contours[[length(contours) + 1]] <- list(
        ts = i, z = fps$z_m[k], xr = f$xr, yr = f$yr,
        ustar = as.numeric(df$ustar_interp[i]), zoL = as.numeric(df$zoL[i]))
    }
  }
  if (!length(contours)) { message("  no footprints"); next }

  r_max <- min(R_CAP, max(vapply(contours,
              function(c) max(sqrt(c$xr^2 + c$yr^2)), 0)) * 1.05)
  grid <- make_polar_grid(r_max, n_r = n_r, n_ang = n_ang)
  G <- polar_mask(nwi_proj, grid, target_crs)
  rw <- grid$r

  rows <- list()
  for (cc in contours) {
    poly <- suppressWarnings(st_make_valid(st_sfc(
      st_polygon(list(cbind(c(cc$xr, cc$xr[1]), c(cc$yr, cc$yr[1])))),
      crs = target_crs)))
    W <- polar_mask(poly, grid, target_crs)
    denom <- sum(W * rw)
    if (!is.finite(denom) || denom <= 0) next

    # f_wet at every wind direction: circular shift of the wetland mask.
    fw <- vapply(seq_len(n_ang), function(k) {
      Gs <- if (k == 1) G else cbind(G[, k:n_ang, drop = FALSE],
                                     G[, seq_len(k - 1), drop = FALSE])
      sum(W * Gs * rw) / denom
    }, numeric(1))

    rows[[length(rows) + 1]] <- tibble(
      site = site, ts = cc$ts, z = cc$z, ustar = cc$ustar, zoL = cc$zoL,
      area = sum(W * rw),            # area proxy, same weighting as f_wet
      fwet_uniform = mean(fw),       # what the symmetric assumption gives
      fwet_sd_dir  = sd(fw),         # spread across wind directions
      fwet_min = min(fw), fwet_max = max(fw))
  }
  if (!length(rows)) next
  out[[site]] <- bind_rows(rows)

  d <- out[[site]]
  message(sprintf("  %d footprints, %d heights", nrow(d), length(unique(d$z))))
  message(sprintf("  f_wet uniform : mean %.4f, SD across footprints %.4f",
                  mean(d$fwet_uniform), sd(d$fwet_uniform)))
  message(sprintf("  f_wet by direction: mean within-footprint SD %.4f  (min %.3f, max %.3f)",
                  mean(d$fwet_sd_dir), mean(d$fwet_min), mean(d$fwet_max)))
  message(sprintf("  ratio directional/size variation: %.1fx",
                  mean(d$fwet_sd_dir) / max(sd(d$fwet_uniform), 1e-9)))

  # Is uniform f_wet just a function of footprint size?
  fit <- tryCatch(summary(lm(fwet_uniform ~ poly(log(area), 2), data = d)),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    message(sprintf("  R^2 of uniform f_wet on footprint size alone: %.3f",
                    fit$r.squared))
  }
}

if (length(out)) {
  all <- bind_rows(out)
  write.csv(all, file.path(tab_dir, "wind_direction_diagnostic.csv"),
            row.names = FALSE)
  message("\nWrote ", file.path(tab_dir, "wind_direction_diagnostic.csv"))
}
