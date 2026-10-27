# footprint_radial.R
# Rotationally-averaged ("symmetric wind") footprint wetness.
#
# We have no wind direction, so the honest assumption is that wind comes from
# all directions with equal probability. Note this is NOT what leaving
# wind_dir = NA does elsewhere in the pipeline: that pins the Kljun footprint
# to one arbitrary fixed bearing, which samples one arbitrary wedge of ground.
# Averaging over all bearings is the defensible version.
#
# Doing that by brute force (rotate the contour N times, intersect each with
# NWI) costs N intersections per footprint. Instead we decompose radially, so
# the expensive wetland geometry is touched once per site:
#
#   g(r) = fraction of the circle of radius r around the tower that is wetland
#   w(r) = fraction of the circle of radius r that falls inside the footprint
#
# g(r) depends only on the site; w(r) only on the footprint. Averaged over all
# rotations, the expected wetland fraction inside the footprint is
#
#   f_wet = Σ w(r) g(r) r  /  Σ w(r) r        (the r factor is the area
#                                              element of the annulus 2πr dr)
#
# which is exact for uniformly-distributed wind direction, not an approximation
# of the rotation average.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})


#' Polar sampling grid centred on the tower, in local LAEA metres.
#'
#' @param r_max Numeric. Outer radius (m).
#' @param n_r Integer. Number of radial bins.
#' @param n_ang Integer. Number of angular samples per radius.
#' @return List with r (n_r), and x/y (n_r * n_ang, angle varying fastest).
make_polar_grid <- function(r_max, n_r = 120, n_ang = 240) {
  r   <- seq(r_max / n_r, r_max, length.out = n_r)
  ang <- seq(0, 2 * pi, length.out = n_ang + 1)[seq_len(n_ang)]
  list(r = r, n_r = n_r, n_ang = n_ang,
       x = as.vector(outer(r, cos(ang))),
       y = as.vector(outer(r, sin(ang))))
}


#' Angular fraction of each circle that is wetland: g(r).
#'
#' @param nwi_proj sf. NWI polygons in the tower-centred LAEA CRS.
#' @param grid List from make_polar_grid().
#' @param target_crs Character. CRS of the grid coordinates.
#' @return Numeric vector, length grid$n_r, each in [0, 1].
wetland_radial_profile <- function(nwi_proj, grid, target_crs) {
  if (is.null(nwi_proj) || nrow(nwi_proj) == 0) {
    return(rep(0, grid$n_r))
  }
  pts <- st_sfc(lapply(seq_along(grid$x),
                       function(i) st_point(c(grid$x[i], grid$y[i]))),
                crs = target_crs)
  nwi_proj <- suppressWarnings(st_make_valid(nwi_proj))
  hit <- lengths(suppressMessages(st_intersects(pts, nwi_proj))) > 0
  # Column-major: angle varies fastest, so each row is one radius.
  rowMeans(matrix(hit, nrow = grid$n_r, ncol = grid$n_ang))
}


#' Angular fraction of each circle that falls inside a footprint contour: w(r).
#'
#' @param ffp List from compute_footprint_kljun().
#' @param grid List from make_polar_grid().
#' @param target_crs Character.
#' @return Numeric vector length grid$n_r, or NULL if the footprint is invalid.
footprint_radial_weight <- function(ffp, grid, target_crs) {
  if (is.null(ffp) || isTRUE(ffp$flag_err == 1) || length(ffp$xr) < 3) {
    return(NULL)
  }
  poly <- st_sfc(
    st_polygon(list(cbind(c(ffp$xr, ffp$xr[1]), c(ffp$yr, ffp$yr[1])))),
    crs = target_crs)
  poly <- suppressWarnings(st_make_valid(poly))
  pts <- st_sfc(lapply(seq_along(grid$x),
                       function(i) st_point(c(grid$x[i], grid$y[i]))),
                crs = target_crs)
  hit <- lengths(suppressMessages(st_intersects(pts, poly))) > 0
  rowMeans(matrix(hit, nrow = grid$n_r, ncol = grid$n_ang))
}


#' Rotationally-averaged wetland fraction seen by a footprint.
#'
#' @param w Numeric vector. Footprint radial weight, from footprint_radial_weight().
#' @param g Numeric vector. Wetland radial profile, from wetland_radial_profile().
#' @param r Numeric vector. Radii (m).
#' @return Scalar f_wet in [0, 1], or NA if the footprint has no area.
fwet_rotational <- function(w, g, r) {
  if (is.null(w)) return(NA_real_)
  denom <- sum(w * r)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  sum(w * g * r) / denom
}


#' Radius containing a given quantile of the rotationally-averaged footprint
#' area — a single interpretable "how far out does this height see" number.
#'
#' @param w Numeric vector. Footprint radial weight.
#' @param r Numeric vector. Radii (m).
#' @param q Numeric. Quantile (default 0.5).
#' @return Radius (m), or NA.
footprint_radius_quantile <- function(w, r, q = 0.5) {
  if (is.null(w)) return(NA_real_)
  cw <- cumsum(w * r)
  if (max(cw) <= 0) return(NA_real_)
  r[which(cw >= q * max(cw))[1]]
}
