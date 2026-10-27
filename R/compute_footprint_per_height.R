# compute_footprint_per_height.R
# Per-height Kljun (2015) footprint computation
#
# Wraps calc_footprint_FFP() (sourced from parent repo) and supplies the inputs
# the parent function expects, parameterizing what NEON aligned data does not
# directly provide:
#   sigma_v (lateral wind sd)  -> 1.9 * u* with mild stability correction
#                                  (Panofsky & Dutton 1984)
#   BLH h                      -> default 1000 m neutral/unstable,
#                                  min(0.3 * u*/f, 200) for stable
#   wind_dir                   -> optional. If NA we keep footprint in the
#                                  Kljun-native +x = downwind orientation;
#                                  rotation to absolute coordinates is handled
#                                  in footprint_wetness.R.
#
# Footprint area derives from the 80% contour. We pick 80% (not 90%) because
# 80% contours are more numerically stable in stable/low-u* conditions where
# the 90% contour can blow up.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

# Sanity bounds. Values that fall outside these are flagged and the footprint
# is not computed (returns flag_err = 1). These are deliberately looser than
# Kljun's published validity ranges so we capture more of the dataset; we'll
# filter further at the analysis stage.
KLJUN_VALID_USTAR_MIN <- 0.1
KLJUN_VALID_ZM_OVER_H_MAX <- 0.8


#' Parameterize σ_v from u* (Panofsky & Dutton; sigma_v / u* ~ 1.9 neutral,
#' rises mildly under unstable conditions).
#'
#' @param ustar Numeric vector. Friction velocity (m/s).
#' @param zoL Numeric vector. Dimensionless stability z/L.
#' @return Numeric vector. σ_v (m/s).
estimate_sigma_v <- function(ustar, zoL) {
  # Unstable: σ_v / u* ≈ 1.9 * (1 - 0.5 * zoL)^(1/3), capped at 3
  # Stable/neutral: 1.9
  out <- 1.9 * ustar
  unstable <- !is.na(zoL) & zoL < 0
  if (any(unstable)) {
    correction <- (1 - 0.5 * pmin(zoL[unstable], 0))^(1/3)
    out[unstable] <- 1.9 * ustar[unstable] * pmin(correction, 1.6)
  }
  out
}


#' Estimate boundary-layer height.
#'   stable (L > 0):  h = min(0.3 * u*/f, 200) m   (Nieuwstadt 1981 limit)
#'   neutral/unstab:  h = 1000 m                  (mid-day continental default)
#'
#' These are rough but adequate for footprint-area calculations — Kljun (2015)
#' note the result is weakly sensitive to BLH for zm/h < ~0.1, which holds at
#' all NEON sites with their <80 m towers and typical BLHs > 500 m.
#'
#' @param L Numeric vector. Obukhov length (m).
#' @param ustar Numeric vector. Friction velocity (m/s).
#' @param lat Numeric. Site latitude (degrees) — used to compute Coriolis.
#' @return Numeric vector of BLH (m).
estimate_blh <- function(L, ustar, lat) {
  f <- 2 * 7.2921e-5 * sin(lat * pi / 180)
  f <- max(abs(f), 1e-5)  # avoid blow-up at equator
  h_stable <- pmin(0.3 * ustar / f, 200)
  h_stable <- pmax(h_stable, 50)
  out <- rep(1000, length(L))
  stable <- !is.na(L) & L > 0
  out[stable] <- h_stable[stable]
  out
}


#' Single-timestamp footprint at a single measurement height.
#'
#' Wraps the parent-repo `calc_footprint_FFP()`. Returns NULL on error.
#'
#' @param zm Numeric. Measurement height above displacement (m), i.e. (z - d).
#' @param z0 Numeric. Roughness length (m).
#' @param umean Numeric. Mean wind at z (m/s). May be NA — we then rely on z0.
#' @param h Numeric. BLH (m).
#' @param ol Numeric. Obukhov length (m).
#' @param sigmav Numeric. σ_v (m/s).
#' @param ustar Numeric. u* (m/s).
#' @param r Numeric. Contour percentage. Default 80.
#' @param nx Integer. Grid resolution. Default 600.
#' @return List with x_2d, y_2d, f_2d, xr (contour x), yr (contour y),
#'         area_80, flag_err.
compute_footprint_kljun <- function(zm, z0, umean, h, ol, sigmav, ustar,
                                      r = 80, nx = 600) {

  # Validate inputs
  if (any(is.na(c(zm, h, ol, sigmav, ustar))) ||
      ustar < KLJUN_VALID_USTAR_MIN ||
      zm <= 0 || zm / h > KLJUN_VALID_ZM_OVER_H_MAX) {
    return(list(flag_err = 1L))
  }

  # rslayer = 1: per Kljun (2015), allows footprint estimation inside the
  # roughness sublayer (canopy + 2-3 canopy heights) — a rough estimate, but
  # essential for our pair-resolved test because the lower sensor on a tall
  # tower frequently sits in the roughness sublayer. We are testing relative
  # area / wetness between pairs, which is more robust than absolute footprint
  # placement.
  ffp <- tryCatch(
    suppressMessages(suppressWarnings(
      calc_footprint_FFP(zm = zm, z0 = z0, umean = umean, h = h, ol = ol,
                          sigmav = sigmav, ustar = ustar,
                          wind_dir = NULL, r = r, nx = nx, rslayer = 1,
                          crop = 1)
    )),
    error = function(e) list(flag_err = 1L)
  )
  if (is.null(ffp) || !is.list(ffp) || isTRUE(ffp$flag_err == 1)) {
    return(list(flag_err = 1L))
  }

  # Extract the r% contour and compute its area in (m^2).
  xr <- if (!is.null(ffp$xr)) ffp$xr[[1]] else NA
  yr <- if (!is.null(ffp$yr)) ffp$yr[[1]] else NA
  area_r <- NA_real_
  if (length(xr) > 2 && length(yr) > 2 &&
      all(is.finite(xr)) && all(is.finite(yr))) {
    poly <- st_polygon(list(cbind(c(xr, xr[1]), c(yr, yr[1]))))
    area_r <- as.numeric(st_area(st_sfc(poly)))
  }

  list(
    x_2d     = ffp$x_2d,
    y_2d     = ffp$y_2d,
    f_2d     = ffp$f_2d,
    xr       = xr,
    yr       = yr,
    area_r   = area_r,
    r_pct    = r,
    flag_err = 0L
  )
}


#' Compute footprints at every measurement height for a single 30-min record.
#'
#' Picks the right level wind from ubar1..ubarN and supplies derived BLH and σ_v.
#'
#' @param met_row A single-row tibble/data.frame from aligned data with at
#'        minimum: ustar_interp, L_obukhov, zoL, roughLength_calc,
#'        z_displ_calc, ubar1..ubarN.
#' @param heights Numeric vector of measurement heights z (m above ground).
#' @param tower_positions Integer vector aligned with `heights` (1..N).
#' @param lat Numeric. Site latitude (deg).
#' @param r Numeric. Contour pct. Default 80.
#' @return Tibble: TowerPosition, z_m, area_r, flag_err, plus list-column
#'         `footprint` (the full list from compute_footprint_kljun).
compute_footprints_all_heights <- function(met_row, heights, tower_positions,
                                            lat, r = 80) {

  ustar <- as.numeric(met_row$ustar_interp[1])
  L     <- as.numeric(met_row$L_obukhov[1])
  zoL_v <- as.numeric(met_row$zoL[1])
  z0    <- as.numeric(met_row$roughLength_calc[1])
  d     <- as.numeric(met_row$z_displ_calc[1])

  sigmav <- estimate_sigma_v(ustar, zoL_v)
  h_blh  <- estimate_blh(L, ustar, lat)

  out <- vector("list", length(heights))
  for (i in seq_along(heights)) {
    z      <- heights[i]
    pos    <- tower_positions[i]
    ubar_col <- paste0("ubar", pos)
    umean <- if (ubar_col %in% names(met_row)) {
      as.numeric(met_row[[ubar_col]][1])
    } else NA_real_

    # zm = effective height above displacement
    zm <- pmax(z - ifelse(is.na(d), 0, d), 0.5)

    ffp <- compute_footprint_kljun(zm = zm, z0 = z0, umean = umean,
                                     h = h_blh, ol = L,
                                     sigmav = sigmav, ustar = ustar, r = r)
    out[[i]] <- tibble(
      TowerPosition = pos,
      z_m           = z,
      zm_eff        = zm,
      area_r        = if (is.null(ffp$area_r)) NA_real_ else ffp$area_r,
      flag_err      = ffp$flag_err,
      footprint     = list(ffp)
    )
  }
  bind_rows(out)
}
