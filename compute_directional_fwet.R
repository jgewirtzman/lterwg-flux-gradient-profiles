# compute_directional_fwet.R
#
# Recompute footprint wetland fraction using the REAL per-timestamp wind
# direction (angZaxsErth) and the MEASURED sigma_v (veloYaxsHorSd), both from
# the NEON 30-min product. Previously wind_dir was NULL — which pinned every
# footprint to one arbitrary bearing — and sigma_v was parameterised as ~1.9*u*.
#
# Why this matters: under symmetric-wind averaging, f_wet is 83-96% explained by
# footprint SIZE alone, so delta_fwet could never be separated from delta_Area,
# and footprint size is itself driven by the same stability/u* that drive the
# flux. Aiming varies independently of size, so the directional f_wet is the
# variable that can actually identify a compositional effect.
#
# Both quantities come from ONE mask per footprint: G (wetland) and W (footprint)
# live on the same polar grid, so every wind direction is a circular shift of the
# angular index. f_wet_dir is the value at the observed bearing; f_wet_sym is the
# mean over all bearings. The mean over shifts is checked against the independent
# radial-method value as a correctness test.
#
# Usage:  Rscript compute_directional_fwet.R [SITE ...]

source("config.R")
suppressPackageStartupMessages({ library(lubridate) })

NWI_KM <- 10
n_ts   <- 500   # raised from 120; the within-site null was likely underpowered
n_r    <- 90
n_ang  <- 180
R_CAP  <- NWI_KM * 1000

nwi_dir  <- file.path(proj_dir, "data", sprintf("nwi_%dkm", NWI_KM))
ck_dir   <- file.path(proj_dir, "output", "footprint", "directional")
dir.create(ck_dir, showWarnings = FALSE, recursive = TRUE)

args  <- commandArgs(trailingOnly = TRUE)
sites <- if (length(args)) args else available_fluxfoot_sites()
message("sites with FluxFoot met: ", length(sites))
set.seed(11)

polar_mask <- function(geom, grid, crs) {
  pts <- st_sfc(lapply(seq_along(grid$x),
                       function(i) st_point(c(grid$x[i], grid$y[i]))), crs = crs)
  matrix(lengths(suppressMessages(st_intersects(pts, geom))) > 0,
         nrow = grid$n_r, ncol = grid$n_ang)
}

for (site in sites) {
  ck <- file.path(ck_dir, paste0(site, ".rds"))
  if (file.exists(ck)) { message("[cached] ", site); next }

  met <- load_fluxfoot_met(site)
  if (is.null(met) || sum(!is.na(met$wind_dir)) < 200) {
    message(site, ": no usable wind direction — skipping"); next }
  nwi <- load_nwi_for_site(site, nwi_dir)
  if (is.null(nwi)) { message(site, ": no 10 km NWI — skipping"); next }
  sd_ <- tryCatch(load_all_sites(site, data_dir), error = function(e) NULL)
  if (is.null(sd_) || is.null(sd_[[site]]$aligned$CH4)) {
    message(site, ": no aligned CH4 — skipping"); next }

  attr_df <- sd_[[site]]$attr
  lat <- as.numeric(attr_df$LatTow[1]); lon <- as.numeric(attr_df$LonTow[1])
  crs <- pick_local_aea(lon, lat)
  nwi <- suppressWarnings(st_make_valid(st_transform(nwi, crs)))

  hl <- attr_df %>% select(TowerPosition, DistZaxsLvlMeasTow) %>% distinct() %>%
    mutate(TowerPosition = as.numeric(TowerPosition),
           z_m = as.numeric(DistZaxsLvlMeasTow)) %>%
    filter(!is.na(z_m), z_m > 0) %>% arrange(z_m)
  if (nrow(hl) < 2) next

  df <- sd_[[site]]$aligned$CH4 %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    inner_join(met %>% filter(!is.na(wind_dir)), by = "time_round")
  df <- df[which(is.finite(as.numeric(df$ustar_interp)) &
                 is.finite(as.numeric(df$L_obukhov)) &
                 is.finite(as.numeric(df$roughLength_calc))), , drop = FALSE]
  if (nrow(df) < 50) { message(site, ": too few joined rows"); next }
  df <- df[sample.int(nrow(df), min(n_ts, nrow(df))), , drop = FALSE]
  message("\n=== ", site, " — ", nrow(df), " timestamps with wind direction")

  # ---- pass 1: contours (measured sigma_v where available) ----
  cont <- list()
  for (i in seq_len(nrow(df))) {
    rec <- df[i, , drop = FALSE]
    sv  <- as.numeric(rec$sigma_v_meas[1])
    ust <- as.numeric(rec$ustar_interp[1]); L <- as.numeric(rec$L_obukhov[1])
    z0  <- as.numeric(rec$roughLength_calc[1]); d <- as.numeric(rec$z_displ_calc[1])
    if (!is.finite(sv)) sv <- estimate_sigma_v(ust, as.numeric(rec$zoL[1]))
    h_blh <- estimate_blh(L, ust, lat)
    for (k in seq_len(nrow(hl))) {
      z <- hl$z_m[k]; pos <- hl$TowerPosition[k]
      ub <- paste0("ubar", pos)
      um <- if (ub %in% names(rec)) as.numeric(rec[[ub]][1]) else NA_real_
      zm <- pmax(z - ifelse(is.na(d), 0, d), 0.5)
      ffp <- compute_footprint_kljun(zm = zm, z0 = z0, umean = um, h = h_blh,
                                     ol = L, sigmav = sv, ustar = ust, r = 80)
      if (isTRUE(ffp$flag_err == 1) || length(ffp$xr) < 3) next
      cont[[length(cont) + 1]] <- list(i = i, z = z, pos = pos,
        tr = rec$time_round[1],
        wd = as.numeric(rec$wind_dir[1]), sv_meas = is.finite(as.numeric(rec$sigma_v_meas[1])),
        xr = ffp$xr, yr = ffp$yr)
    }
  }
  if (!length(cont)) { message("  no valid footprints"); next }

  reach <- max(vapply(cont, function(c) max(sqrt(c$xr^2 + c$yr^2)), 0))
  grid  <- make_polar_grid(min(R_CAP, max(reach * 1.05, 200)), n_r = n_r, n_ang = n_ang)
  G     <- polar_mask(nwi, grid, crs)
  g_rad <- rowMeans(G)          # for the independent cross-check
  rw    <- grid$r
  dtheta <- 360 / n_ang

  rows <- list()
  for (cc in cont) {
    poly <- suppressWarnings(st_make_valid(st_sfc(
      st_polygon(list(cbind(c(cc$xr, cc$xr[1]), c(cc$yr, cc$yr[1])))), crs = crs)))
    W <- polar_mask(poly, grid, crs)
    den <- sum(W * rw); if (!is.finite(den) || den <= 0) next

    fw <- vapply(seq_len(n_ang), function(k) {
      Gs <- if (k == 1) G else cbind(G[, k:n_ang, drop = FALSE],
                                     G[, seq_len(k - 1), drop = FALSE])
      sum(W * Gs * rw) / den
    }, numeric(1))

    # Rotation convention — verified against a synthetic wetland patch at known
    # bearings (see the rotcheck test): the flux footprint is the SOURCE area,
    # which lies UPWIND, so Kljun's native +x axis must be rotated to point
    # toward the direction the wind comes FROM. Kljun's own rotation code maps
    # native +x to compass bearing wind_dir, confirming this. As a math angle
    # (ccw from east) that bearing is 90 - wind_dir. Rotating the footprint by
    # phi is equivalent to shifting the wetland field by +phi.
    #
    # The earlier version used phi = -90 - wind_dir with a -phi shift, which is
    # a REFLECTION about the N-S axis (error = 2*wind_dir), not a 180 deg flip:
    # a wetland at bearing 45 peaked at 305 instead of 45. That invalidated the
    # first directional run entirely.
    phi <- (90 - cc$wd) %% 360
    kk  <- (round(phi / dtheta) %% n_ang) + 1
    rows[[length(rows) + 1]] <- tibble(
      site = site, ts = cc$i, time_round = cc$tr, z = cc$z, TowerPosition = cc$pos,
      wind_dir = cc$wd, sigma_v_measured = cc$sv_meas,
      f_wet_dir = fw[kk], f_wet_sym = mean(fw),
      f_wet_sym_radial = fwet_rotational(rowMeans(W), g_rad, rw),
      area = den, fw_min = min(fw), fw_max = max(fw))
  }
  if (!length(rows)) next
  r <- bind_rows(rows)
  saveRDS(r, ck)

  chk <- suppressWarnings(cor(r$f_wet_sym, r$f_wet_sym_radial, use = "complete.obs"))
  message(sprintf("  %d footprints | sym-vs-radial cross-check r = %.4f", nrow(r), chk))
  message(sprintf("  f_wet directional: mean %.4f, sd %.4f | symmetric: mean %.4f, sd %.4f",
                  mean(r$f_wet_dir), sd(r$f_wet_dir), mean(r$f_wet_sym), sd(r$f_wet_sym)))
  fit_d <- tryCatch(summary(lm(f_wet_dir ~ poly(log(area), 2), data = r))$r.squared, error = function(e) NA)
  fit_s <- tryCatch(summary(lm(f_wet_sym ~ poly(log(area), 2), data = r))$r.squared, error = function(e) NA)
  message(sprintf("  R^2 on footprint SIZE alone -> symmetric %.3f | directional %.3f",
                  fit_s, fit_d))
}
message("\nCheckpoints: ", ck_dir)
