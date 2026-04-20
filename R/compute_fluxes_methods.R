# compute_fluxes_methods.R
# Head-to-head flux computation using the three non-MBR methods from the
# parent lterwg-flux-gradient repo, adapted to our 30-min per-level data.
#
# For each adjacent-sensor pair at each timestamp we compute:
#   F_KZ  — our K(z) model (MOST + Raupach attenuation with site-specific α)
#   F_AE  — aerodynamic method: K = k·u*·z_eff / φ_h, z_eff geometric
#   F_WP  — wind-profile method: K = k²·ū·z_geom / (ln(z_geom/z₀)·φ_h)
#
# Common inputs (per 30-min window, per adjacent-sensor pair):
#   C_lower, C_upper   concentrations at the two sensors
#   z_lower, z_upper   heights (m)
#   u*                 friction velocity from sonic (m/s)
#   L                  Obukhov length (m, = z_eff/zoL)
#   ū                  wind speed at the upper sensor (or top)
#   z₀                 roughness length
#   d                  displacement height (0.66·h by default)
#
# Returns a long-format tibble with F_KZ, F_AE, F_WP for each pair.


#' Stability correction for heat/scalars, following the parent repo's
#' formulation (same as MOST used elsewhere in this project).
phih_stability <- function(zeta) {
  ifelse(
    is.na(zeta) | !is.finite(zeta), 1,
    ifelse(zeta >= 0,           # stable
           1 + 5 * zeta,
           (1 - 16 * zeta)^(-0.5)  # unstable
    )
  )
}


#' Aerodynamic-method K for a sensor pair.
#'
#' K_AE = k · u* · z_eff / φ_h(zeta)
#' z_eff = geometric mean of (z_lower − d, z_upper − d)
#' zeta  = z_eff / L
K_AE <- function(z_lower, z_upper, u_star, L, disp_h) {
  z_lower_eff <- pmax(z_lower - disp_h, 0.01)
  z_upper_eff <- pmax(z_upper - disp_h, 0.01)
  z_eff       <- sqrt(z_lower_eff * z_upper_eff)
  zeta        <- ifelse(is.finite(L) & L != 0, z_eff / L, 0)
  phi         <- pmax(phih_stability(zeta), 0.1)
  pmax(0.4 * u_star * z_eff / phi, 1e-4)
}


#' Wind-profile-method K for a sensor pair.
#'
#' K_WP = k² · ū · z_geom / (ln(z_geom / z₀) · φ_h)
#' where z_geom = geometric mean of (z_lower, z_upper) above displacement,
#' ū is mean wind speed at the UPPER sensor (a convention used in parent repo).
K_WP <- function(z_lower, z_upper, u_bar, z0, L, disp_h) {
  z_lower_eff <- pmax(z_lower - disp_h, 0.01)
  z_upper_eff <- pmax(z_upper - disp_h, 0.01)
  z_geom      <- sqrt(z_lower_eff * z_upper_eff)
  zeta        <- ifelse(is.finite(L) & L != 0, z_geom / L, 0)
  phi         <- pmax(phih_stability(zeta), 0.1)
  denom       <- log(pmax(z_geom / pmax(z0, 1e-3), 1.1)) * phi
  pmax(0.4^2 * u_bar * z_geom / denom, 1e-4)
}


#' Compute AE and WP fluxes for an adjacent-sensor pair at each 30-min window.
#'
#' Uses the same input schema as compute_layer_fluxes (profile_ts + met_ts +
#' attr_df). The existing K(z) flux is already computed and stored in
#' `all_fluxes[[site]][[gas]]`; this function adds F_AE and F_WP columns
#' (joined by time_round, z_lower, z_upper).
#'
#' @param flux_df  Existing tibble from compute_layer_fluxes() for ONE site/gas.
#' @param met_ts   Output of extract_met_timeseries(): timeMid, ustar, zoL,
#'                 roughLength, ubar_top, etc.
#' @param winds_ts Output of extract_wind_profile(): timeMid (time_round),
#'                 TowerPosition, height_m, ubar.
#' @param attr_df  Site attributes (for canopy height / displacement).
#' @return The same tibble with added K_AE, K_WP, F_AE, F_WP columns.
add_ae_wp_fluxes <- function(flux_df, met_ts, winds_ts, attr_df) {

  if (is.null(flux_df) || nrow(flux_df) == 0) return(flux_df)
  canopy_h <- get_canopy_height(attr_df)
  disp_h   <- if (is.na(canopy_h) || canopy_h < 1) 0 else 0.66 * canopy_h

  # Build a (time_round, height_m) -> ubar lookup from the wind profile
  # for the UPPER sensor of each pair.
  w <- winds_ts %>%
    mutate(time_round = if ("time_round" %in% names(.)) time_round else timeMid) %>%
    select(time_round, height_m, ubar)

  # Pull z0 from met (if present) — parent repo uses roughLength_interp
  met <- met_ts
  if (!"roughLength" %in% names(met) && "roughLength_interp" %in% names(met)) {
    met <- met %>% rename(roughLength = roughLength_interp)
  }
  if (!"roughLength" %in% names(met)) met$roughLength <- 0.1 * max(canopy_h, 0.1)

  met_30 <- met %>%
    mutate(time_round = timeMid) %>%
    select(time_round, ustar, zoL, roughLength, any_of("rhoa_kgm3"))

  # Join wind speed at upper sensor of each pair
  flux_df <- flux_df %>%
    left_join(w %>% rename(z_upper = height_m, u_upper = ubar),
              by = c("time_round", "z_upper")) %>%
    left_join(met_30 %>% select(time_round, roughLength, any_of("rhoa_kgm3")),
              by = "time_round")

  # L per row: use L_obukhov from stage-4 if present (avoids z_upper/zoL error).
  # rho_mol per row: from stage-4 rhoa_kgm3 (kg/m³) / M_air (kg/mol).
  has_L   <- "L_obukhov" %in% names(flux_df)
  has_rho <- "rhoa_kgm3" %in% names(flux_df)
  flux_df <- flux_df %>%
    mutate(
      L       = if (has_L) L_obukhov else ifelse(!is.na(zoL) & zoL != 0, z_upper / zoL, NA_real_),
      rho_mol = if (has_rho) rhoa_kgm3 / 0.028964 else 42.3,
      K_AE = K_AE(z_lower, z_upper, ustar, L, disp_h),
      K_WP = K_WP(z_lower, z_upper, u_upper, roughLength, L, disp_h),
      F_AE = -K_AE * dCdz * rho_mol,
      F_WP = -K_WP * dCdz * rho_mol
    )
  flux_df
}
