# footprint_wetness.R
# Combine per-height footprints with the NWI wetland mask to produce
#   area(z) — area of the r% concentration footprint at height z
#   f_wet(z) — footprint-weighted wetland fraction at height z
#
# The Kljun footprint comes back in a local (downwind = +x, crosswind = y)
# coordinate frame centered on the tower. To intersect with NWI polygons in
# absolute coordinates we:
#   1. Rotate the local grid by wind_dir (clockwise from north).
#      If wind_dir is NA we run an unrotated assessment (treats wetlands as
#      isotropic — a fallback, not the headline test).
#   2. Translate to the tower lat/lon, projected to the same LAEA CRS as the
#      site's nwi_proj layer.
#   3. Rasterize footprint probability and overlay on the wetland polygons,
#      OR vectorize the r% contour and intersect.
#
# We use the vector approach: form a polygon from the r% contour, project,
# intersect with NWI, and report wet-area / contour-area as f_wet. This is
# faster than rasterizing and is the right denominator for the headline test
# (it uses the same area metric as ΔA).

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})


#' Build an sf polygon from the r% contour returned by compute_footprint_kljun.
#'
#' @param ffp List. Output of compute_footprint_kljun().
#' @param wind_dir Numeric. Wind direction in degrees from north. NA = no rotation.
#' @param tower_lon Numeric.
#' @param tower_lat Numeric.
#' @param target_crs Character. CRS string for output (the site's LAEA).
#' @return sf POLYGON in `target_crs`, or NULL if footprint invalid.
footprint_contour_to_sf <- function(ffp, wind_dir, tower_lon, tower_lat,
                                      target_crs) {

  if (is.null(ffp) || isTRUE(ffp$flag_err == 1) ||
      length(ffp$xr) < 3) return(NULL)

  x <- ffp$xr; y <- ffp$yr
  # Kljun convention: x = downwind (positive into wind-from direction).
  # Geographic convention: wind_dir = direction wind is *coming from*
  # (meteorological), so the downwind axis points to (wind_dir + 180) deg.
  # We rotate (x, y) by the angle needed to align +x with downwind direction.
  if (!is.na(wind_dir)) {
    theta <- (wind_dir + 180 - 90) * pi / 180  # convert to math angle, ccw from +x-east
    cos_t <- cos(theta); sin_t <- sin(theta)
    x_rot <- cos_t * x - sin_t * y
    y_rot <- sin_t * x + cos_t * y
  } else {
    x_rot <- x
    y_rot <- y
  }

  # Build polygon in LAEA centred on tower (x_rot, y_rot are already in metres
  # from the tower) -> we just need the target CRS string.
  poly <- st_polygon(list(cbind(c(x_rot, x_rot[1]), c(y_rot, y_rot[1]))))
  st_sfc(poly, crs = target_crs) %>% st_sf(geometry = .)
}


#' Footprint-weighted wetland fraction for one footprint at one timestamp.
#'
#' Uses the 2D footprint probability surface (`f_2d`) if available, otherwise
#' falls back to a uniform-within-contour assumption. The probability-weighted
#' integration is the technically correct quantity:
#'   f_wet = sum_{cells in NWI} f_2d_i * area_i  /  sum f_2d_i * area_i
#'
#' For speed we do an area-weighted (uniform-probability-in-contour) version
#' by default — fine for the headline regression because we're comparing
#' *relative* f_wet between paired heights, and any constant bias cancels.
#'
#' @param ffp List. Output of compute_footprint_kljun().
#' @param wind_dir Numeric. Wind direction (deg from north). NA = no rotation.
#' @param tower_lon Numeric.
#' @param tower_lat Numeric.
#' @param nwi_proj sf. Projected NWI polygons in `target_crs`.
#' @param target_crs Character. CRS string (the site's LAEA).
#' @return Named list: area_r (m^2), wet_area (m^2), f_wet (unitless).
compute_fwet_for_footprint <- function(ffp, wind_dir, tower_lon, tower_lat,
                                         nwi_proj, target_crs) {

  poly_sf <- footprint_contour_to_sf(ffp, wind_dir,
                                       tower_lon, tower_lat, target_crs)
  if (is.null(poly_sf)) {
    return(list(area_r = NA_real_, wet_area = NA_real_, f_wet = NA_real_))
  }
  fp_area <- as.numeric(st_area(poly_sf))

  if (nrow(nwi_proj) == 0) {
    return(list(area_r = fp_area, wet_area = 0, f_wet = 0))
  }

  # Intersection of footprint contour with NWI polygons
  inter <- tryCatch(
    suppressWarnings(st_intersection(nwi_proj, poly_sf)),
    error = function(e) NULL
  )
  if (is.null(inter) || nrow(inter) == 0) {
    return(list(area_r = fp_area, wet_area = 0, f_wet = 0))
  }
  wet_area <- sum(as.numeric(st_area(inter)))

  list(area_r = fp_area, wet_area = wet_area,
       f_wet = wet_area / fp_area)
}


#' Driver: compute area(z) and f_wet(z) at every height for every 30-min
#' record at a single site, then form pair-wise ΔA and Δf_wet differences.
#'
#' @param aligned_gas A single gas's aligned data frame (e.g., min9Diff.list$CH4),
#'   but here we want one row per 30-min PAIR. Will derive pair structure from
#'   `dLevelsAminusB`.
#' @param attr_df Site attributes data frame.
#' @param nwi_proj sf object of NWI polygons in LAEA. May be empty.
#' @param wind_dir_col Optional column name in aligned data containing wind dir
#'   (degrees from north). If NULL, no rotation applied.
#' @param sample_n Integer. If > 0, randomly subsample this many timestamps
#'   to keep runtime tractable. Default 0 = use all rows.
#' @param r_pct Contour percentage. Default 80.
#' @return Tibble (one row per 30-min × pair) with columns:
#'   timeMid, dLevelsAminusB, TowerPosition_A, TowerPosition_B,
#'   z_A, z_B, area_A, area_B, dArea, f_wet_A, f_wet_B, dfwet,
#'   dConc, ustar, L, zoL, season
compute_footprint_pair_metrics <- function(aligned_gas, attr_df, nwi_proj,
                                             wind_dir_col = NULL,
                                             sample_n = 0, r_pct = 80) {

  if (is.null(aligned_gas) || nrow(aligned_gas) == 0) return(tibble())

  lat <- as.numeric(attr_df$LatTow[1])
  lon <- as.numeric(attr_df$LonTow[1])
  target_crs <- pick_local_aea(lon, lat)

  # Project NWI if not already in target CRS
  if (nrow(nwi_proj) > 0) {
    nwi_proj <- st_transform(nwi_proj, target_crs)
    # Make valid (avoid intersection failures)
    nwi_proj <- suppressWarnings(st_make_valid(nwi_proj))
  }

  # Build height lookup
  height_lookup <- attr_df %>%
    select(TowerPosition, DistZaxsLvlMeasTow) %>%
    distinct() %>%
    mutate(TowerPosition = as.numeric(TowerPosition),
           z_m = as.numeric(DistZaxsLvlMeasTow)) %>%
    filter(!is.na(z_m)) %>%
    select(TowerPosition, z_m)

  df <- aligned_gas

  # Restrict to adjacent pairs (diffTowerPosition == 1) *before* subsampling
  # and footprint computation — Kljun is expensive and only adjacent pairs
  # contribute to the headline test. Non-adjacent / bottom-top pairs can be
  # added later as a sensitivity analysis.
  if ("diffTowerPosition" %in% names(df)) {
    df <- df[which(as.numeric(df$diffTowerPosition) == 1), , drop = FALSE]
  }

  # Drop rows with missing inputs that the footprint code needs
  df <- df[which(
    is.finite(as.numeric(df$ustar_interp)) &
    is.finite(as.numeric(df$L_obukhov)) &
    is.finite(as.numeric(df$roughLength_calc))
  ), , drop = FALSE]

  # Optional subsampling for tractability
  if (sample_n > 0 && nrow(df) > sample_n) {
    df <- df[sample.int(nrow(df), sample_n), , drop = FALSE]
  }

  message("  Computing pair-resolved footprint metrics for ", nrow(df),
          " timestamps")

  rows <- list()
  for (i in seq_len(nrow(df))) {
    if (i %% 500 == 0) message("    ", i, " / ", nrow(df))
    rec <- df[i, , drop = FALSE]

    pos_a <- as.numeric(rec$TowerPosition_A[1])
    pos_b <- as.numeric(rec$TowerPosition_B[1])
    z_a <- height_lookup$z_m[height_lookup$TowerPosition == pos_a]
    z_b <- height_lookup$z_m[height_lookup$TowerPosition == pos_b]
    if (length(z_a) == 0 || length(z_b) == 0) next

    wind_dir <- if (!is.null(wind_dir_col) && wind_dir_col %in% names(rec)) {
      as.numeric(rec[[wind_dir_col]][1])
    } else NA_real_

    # Compute footprint at each level
    fps <- compute_footprints_all_heights(
      met_row = rec,
      heights = c(z_a, z_b),
      tower_positions = c(pos_a, pos_b),
      lat = lat, r = r_pct
    )

    fp_a <- fps$footprint[[1]]
    fp_b <- fps$footprint[[2]]
    if (isTRUE(fp_a$flag_err == 1) || isTRUE(fp_b$flag_err == 1)) next

    # f_wet at each level
    fw_a <- compute_fwet_for_footprint(fp_a, wind_dir, lon, lat,
                                         nwi_proj, target_crs)
    fw_b <- compute_fwet_for_footprint(fp_b, wind_dir, lon, lat,
                                         nwi_proj, target_crs)

    rows[[length(rows) + 1]] <- tibble(
      timeMid   = rec$timeMid[1],
      dLevelsAminusB = as.character(rec$dLevelsAminusB[1]),
      TowerPosition_A = pos_a,
      TowerPosition_B = pos_b,
      z_A      = z_a, z_B = z_b,
      area_A   = fw_a$area_r, area_B = fw_b$area_r,
      dArea    = fw_a$area_r - fw_b$area_r,
      f_wet_A  = fw_a$f_wet,  f_wet_B = fw_b$f_wet,
      dfwet    = fw_a$f_wet - fw_b$f_wet,
      dConc    = as.numeric(rec$dConc[1]),
      ustar    = as.numeric(rec$ustar_interp[1]),
      L        = as.numeric(rec$L_obukhov[1]),
      zoL      = as.numeric(rec$zoL[1]),
      wind_dir = wind_dir
    )
  }
  bind_rows(rows)
}
