# calc_storage_flux.R
# Compute layer storage flux (dC/dt) per tower level
#
# The storage flux at each level captures the rate of concentration change over
# time at that height. This is less sensitive to temporal aliasing than gradient-
# based approaches because it uses consecutive measurements at the SAME level,
# not simultaneous measurements at different levels.

#' Calculate dC/dt at each tower level
#'
#' For each TowerPosition, computes the concentration change between consecutive
#' measurements at that level. The time between consecutive same-level measurements
#' is approximately 9min * n_levels (the Picarro cycle time).
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries() with columns:
#'   timeMid, TowerPosition, concentration, vari, numSamp
#' @param gas Character. Gas name (for units labeling)
#' @param max_gap_factor Numeric. Flag gaps where dt exceeds this multiple of the
#'   median dt for that level. Default 1.5.
#' @return Tibble with columns: timeMid, TowerPosition, dCdt (units/sec),
#'   dCdt_per_hr (units/hr), dt_sec, gap_flag
calc_storage_flux_level <- function(profile_ts, gas = "CH4", max_gap_factor = 1.5) {

  if (nrow(profile_ts) == 0) return(tibble())

  profile_ts %>%
    arrange(TowerPosition, timeMid) %>%
    group_by(TowerPosition) %>%
    mutate(
      # Time difference to next measurement at same level (seconds)
      dt_sec = as.numeric(difftime(lead(timeMid), timeMid, units = "secs")),
      # Concentration difference
      dC = lead(concentration) - concentration,
      # Rate of change
      dCdt = dC / dt_sec,
      dCdt_per_hr = dCdt * 3600,
      # Midpoint time for the derivative
      timeMid_deriv = timeMid + seconds(dt_sec / 2)
    ) %>%
    # Remove last row of each group (no lead value)
    filter(!is.na(dCdt)) %>%
    # Flag gaps: compute median dt per level, flag if dt exceeds threshold
    mutate(
      median_dt = median(dt_sec, na.rm = TRUE),
      gap_flag = dt_sec > max_gap_factor * median_dt
    ) %>%
    ungroup() %>%
    transmute(
      timeMid = timeMid_deriv,
      TowerPosition,
      dCdt,
      dCdt_per_hr,
      dt_sec,
      gap_flag,
      gas = gas
    )
}


#' Summarize storage flux as diel pattern by level and season
#'
#' @param storage_flux Tibble. Output of calc_storage_flux_level()
#' @param attr_df Data frame. Site attributes for height lookup
#' @return Tibble with mean/se of dCdt_per_hr by hour, season, level/height
summarize_storage_flux_diel <- function(storage_flux, attr_df) {

  if (nrow(storage_flux) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)

  storage_flux %>%
    # Exclude flagged gaps
    filter(!gap_flag) %>%
    # Add temporal context
    add_season("timeMid") %>%
    add_tod("timeMid") %>%
    mutate(hour = hour(timeMid)) %>%
    # Join heights
    left_join(height_lookup, by = "TowerPosition") %>%
    # Aggregate
    group_by(hour, season, tod, TowerPosition, height_m, gas) %>%
    summarise(
      dCdt_mean = mean(dCdt_per_hr, na.rm = TRUE),
      dCdt_se = sd(dCdt_per_hr, na.rm = TRUE) / sqrt(sum(!is.na(dCdt_per_hr))),
      n = sum(!is.na(dCdt_per_hr)),
      .groups = "drop"
    )
}


#' Compute column-integrated storage flux between two heights
#'
#' Integrates dC/dt through the vertical column using trapezoidal integration.
#' This gives the total storage change rate for a layer (units: ppb*m/hr or ppm*m/hr).
#'
#' @param storage_flux Tibble. Output of calc_storage_flux_level()
#' @param attr_df Data frame. Site attributes
#' @param z_bottom Numeric. Bottom height of integration (m). NULL = lowest level.
#' @param z_top Numeric. Top height of integration (m). NULL = highest level.
#' @return Tibble with timeMid and integrated_storage_flux
calc_integrated_storage <- function(storage_flux, attr_df,
                                     z_bottom = NULL, z_top = NULL) {

  height_lookup <- get_height_lookup(attr_df)

  df <- storage_flux %>%
    filter(!gap_flag) %>%
    left_join(height_lookup, by = "TowerPosition") %>%
    filter(!is.na(height_m))

  # Apply height bounds

  if (!is.null(z_bottom)) df <- df %>% filter(height_m >= z_bottom)
  if (!is.null(z_top))    df <- df %>% filter(height_m <= z_top)

  # For each timestamp, do trapezoidal integration of dCdt over height
  df %>%
    group_by(timeMid) %>%
    arrange(height_m) %>%
    summarise(
      n_levels = n(),
      integrated_storage = if (n() >= 2) {
        # Trapezoidal rule
        sum(diff(height_m) * (dCdt_per_hr[-n()] + dCdt_per_hr[-1]) / 2)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    filter(!is.na(integrated_storage), n_levels >= 2)
}


#' Compute 30-min layer storage flux for each adjacent sensor pair
#'
#' Aggregates 9-min dC/dt (from calc_storage_flux_level) to 30-min windows,
#' then integrates over each inter-sensor layer using rho_mol and Δz.
#' Returns storage_flux in the same units as F_FG from compute_layer_fluxes()
#' (µmol/m²/s for CO2; nmol/m²/s for CH4 — matching concentration units).
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param met_ts Tibble or NULL. Output of extract_met_timeseries(). Used for
#'   rho_mol (rhoa_kgm3); falls back to 42.3 mol/m³ if absent.
#' @param gas Character. Gas name (passed to calc_storage_flux_level)
#' @return Tibble: time_round, pos_lower, pos_upper, z_lower, z_upper, dz,
#'   dCdt_lower, dCdt_upper, dCdt_layer, storage_flux
compute_layer_storage_flux <- function(profile_ts, attr_df, met_ts = NULL,
                                        gas = "CH4") {

  if (nrow(profile_ts) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)
  positions <- sort(unique(height_lookup$TowerPosition))
  if (length(positions) < 2) return(tibble())

  # Per-level dC/dt from 9-min data, gap-flagged
  raw_dcdt <- calc_storage_flux_level(profile_ts, gas) %>%
    filter(!gap_flag)

  if (nrow(raw_dcdt) == 0) return(tibble())

  # Aggregate to 30-min: mean dC/dt per level
  dcdt_30 <- raw_dcdt %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    group_by(time_round, TowerPosition) %>%
    summarise(
      dCdt_mean = mean(dCdt, na.rm = TRUE),
      n_obs     = sum(!is.na(dCdt)),
      .groups   = "drop"
    ) %>%
    filter(!is.na(dCdt_mean)) %>%
    left_join(height_lookup, by = "TowerPosition")

  # rho_mol per 30-min window (mol/m³)
  if (!is.null(met_ts) && "rhoa_kgm3" %in% names(met_ts)) {
    rho_30 <- met_ts %>%
      mutate(time_round = timeMid,
             rho_mol    = rhoa_kgm3 / 0.028964) %>%
      select(time_round, rho_mol)
  } else {
    rho_30 <- tibble(time_round = unique(dcdt_30$time_round), rho_mol = 42.3)
  }

  # Layer storage flux for each adjacent sensor pair
  layers <- vector("list", length(positions) - 1)
  for (i in seq_len(length(positions) - 1)) {
    pos_lo <- positions[i]
    pos_hi <- positions[i + 1]
    z_lo   <- height_lookup$height_m[height_lookup$TowerPosition == pos_lo]
    z_hi   <- height_lookup$height_m[height_lookup$TowerPosition == pos_hi]

    lo <- dcdt_30 %>% filter(TowerPosition == pos_lo) %>%
      select(time_round, dCdt_lower = dCdt_mean)
    hi <- dcdt_30 %>% filter(TowerPosition == pos_hi) %>%
      select(time_round, dCdt_upper = dCdt_mean)

    layers[[i]] <- inner_join(lo, hi, by = "time_round") %>%
      left_join(rho_30, by = "time_round") %>%
      mutate(
        pos_lower    = pos_lo,
        pos_upper    = pos_hi,
        z_lower      = z_lo,
        z_upper      = z_hi,
        dz           = z_hi - z_lo,
        dCdt_layer   = (dCdt_lower + dCdt_upper) / 2,   # trapezoidal
        # rho_mol [mol/m³] × dz [m] × dCdt [conc_unit/mol_air/s]
        # → conc_unit/m²/s  (µmol/m²/s for CO2, nmol/m²/s for CH4)
        storage_flux = rho_mol * dz * dCdt_layer
      ) %>%
      select(time_round, pos_lower, pos_upper, z_lower, z_upper, dz,
             dCdt_lower, dCdt_upper, dCdt_layer, storage_flux)
  }

  bind_rows(layers)
}
