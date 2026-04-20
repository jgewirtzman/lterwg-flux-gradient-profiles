# calc_profile_diffs.R
# Compute pairwise concentration differences from raw 9-min profile means
#
# Instead of using the aligned paired-measurement dConc (which is attenuated
# by QC filtering and temporal misalignment), this computes concentration
# differences from the mean profile at each height under different conditions.
# This gives the true mean gradient structure.

#' Compute pairwise concentration differences from raw profile data
#'
#' Bins raw 9-min concentration data by condition (season x TOD, or
#' season x stability), computes mean concentration at each height,
#' then takes all pairwise differences.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param met_ts Tibble or NULL. Output of extract_met_timeseries().
#'   If provided, uses PAR-based day/night instead of fixed hours.
#' @param group_vars Character vector. Condition variables to group by.
#'   Default: c("season", "tod")
#' @return Tibble with columns: group_vars, height_upper, height_lower,
#'   position_upper, position_lower, conc_upper, conc_lower, dConc, n_upper, n_lower
calc_profile_diffs <- function(profile_ts, attr_df, gas = "CH4",
                                met_ts = NULL,
                                group_vars = c("season", "tod")) {

  if (nrow(profile_ts) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)

  # Add heights and temporal grouping columns
  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    add_season("timeMid")

  # Use PAR-based day/night if met data available, otherwise fixed hours
  if (!is.null(met_ts) && "PAR" %in% names(met_ts)) {
    df <- df %>%
      mutate(time_round = round_date(timeMid, "30 minutes")) %>%
      left_join(met_ts %>% select(timeMid, PAR) %>% rename(time_round = timeMid),
                by = "time_round") %>%
      add_tod_par("PAR") %>%
      select(-time_round, -PAR)
  } else {
    df <- df %>% add_tod("timeMid")
  }

  # Compute mean concentration per height per condition
  height_means <- df %>%
    group_by(across(all_of(group_vars)), TowerPosition, height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      sd_conc = sd(concentration, na.rm = TRUE),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  # Generate all pairwise combinations within each condition group
  height_means %>%
    inner_join(height_means,
               by = group_vars,
               suffix = c("_upper", "_lower"),
               relationship = "many-to-many") %>%
    # Keep only pairs where upper is physically higher
    filter(height_m_upper > height_m_lower) %>%
    mutate(
      dConc = mean_conc_upper - mean_conc_lower,
      is_adjacent = abs(TowerPosition_upper - TowerPosition_lower) == 1,
      gas = gas
    ) %>%
    rename(
      height_upper = height_m_upper,
      height_lower = height_m_lower,
      position_upper = TowerPosition_upper,
      position_lower = TowerPosition_lower,
      conc_upper = mean_conc_upper,
      conc_lower = mean_conc_lower
    ) %>%
    select(all_of(group_vars),
           position_upper, position_lower,
           height_upper, height_lower,
           conc_upper, conc_lower,
           dConc, is_adjacent,
           n_upper, n_lower,
           sd_conc_upper, sd_conc_lower,
           gas)
}


#' Extract met/stability time series from aligned data
#'
#' Pulls stability and turbulence info from the aligned data as a standalone
#' time series (not tied to specific level pairs), for use as a binning variable.
#'
#' @param aligned_data Data frame. One gas from min9Diff.list
#' @return Tibble with one row per ~30-min period: timeMid, zoL, ustar, stability
extract_met_timeseries <- function(aligned_data) {

  if (is.null(aligned_data) || nrow(aligned_data) == 0) return(tibble())

  # Find stability column
  zol_col <- if ("zoL" %in% names(aligned_data)) "zoL" else
             if ("MO.param" %in% names(aligned_data)) "MO.param" else NULL

  df <- aligned_data %>%
    mutate(timeMid = if ("timeMid" %in% names(.)) {
      ymd_hms(timeMid, quiet = TRUE)
    } else if ("timeEnd_A" %in% names(.)) {
      ymd_hms(timeEnd_A, quiet = TRUE)
    } else NA_real_)

  # Check for PAR column
  has_par <- "PAR" %in% names(aligned_data)

  # Roughness length column: prefer NEON-interpolated, fall back to calc
  z0_col <- if ("roughLength_interp" %in% names(aligned_data)) "roughLength_interp" else
            if ("roughLength_calc"   %in% names(aligned_data)) "roughLength_calc"   else NULL

  # Air density: prefer pre-computed rhoa_kgm3 from Stage 4; fall back to NA.
  has_rho  <- "rhoa_kgm3"  %in% names(aligned_data)
  # Obukhov length in metres (preferred over zoL — avoids the z_tower / z_i error).
  has_L    <- "L_obukhov"  %in% names(aligned_data)

  # Take one value per ~30-min window (use top-level pair to avoid duplicates)
  result <- df %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    group_by(time_round) %>%
    summarise(
      ustar       = median(ustar_interp, na.rm = TRUE),
      zoL         = if (!is.null(zol_col)) median(.data[[zol_col]], na.rm = TRUE) else NA_real_,
      L_obukhov   = if (has_L)   median(L_obukhov,  na.rm = TRUE) else NA_real_,
      rhoa_kgm3   = if (has_rho) median(rhoa_kgm3,  na.rm = TRUE) else NA_real_,
      roughLength = if (!is.null(z0_col))
                      median(.data[[z0_col]], na.rm = TRUE) else NA_real_,
      PAR         = if (has_par) median(PAR, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    rename(timeMid = time_round) %>%
    filter(!is.na(timeMid))

  if (!is.null(zol_col)) {
    result <- result %>% add_stability_class("zoL")
  }

  result
}


#' Extract per-level wind speed profile from aligned data
#'
#' Aligned data has columns ubar1, ubar2, ..., ubarN giving horizontal
#' wind speed at each tower level. Pivots to long format and joins heights.
#'
#' @param aligned_data Data frame. One gas from min9Diff.list
#' @param attr_df Data frame. Site attributes
#' @return Tibble with columns: timeMid, TowerPosition, height_m, ubar
extract_wind_profile <- function(aligned_data, attr_df) {

  if (is.null(aligned_data) || nrow(aligned_data) == 0) return(tibble())

  wind_cols <- grep("^ubar[0-9]+$", names(aligned_data), value = TRUE)
  if (length(wind_cols) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)

  df <- aligned_data %>%
    mutate(timeMid = if ("timeMid" %in% names(.)) {
      ymd_hms(timeMid, quiet = TRUE)
    } else if ("timeEnd_A" %in% names(.)) {
      ymd_hms(timeEnd_A, quiet = TRUE)
    } else NA_real_) %>%
    mutate(time_round = round_date(timeMid, "30 minutes"))

  # Aggregate to 30-min windows, keeping one value per level
  df %>%
    select(time_round, all_of(wind_cols)) %>%
    group_by(time_round) %>%
    summarise(across(all_of(wind_cols), ~median(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    pivot_longer(cols = all_of(wind_cols),
                 names_to = "ubar_col",
                 values_to = "ubar") %>%
    mutate(TowerPosition = as.numeric(sub("ubar", "", ubar_col))) %>%
    left_join(height_lookup, by = "TowerPosition") %>%
    filter(!is.na(ubar), !is.na(height_m)) %>%
    select(timeMid = time_round, TowerPosition, height_m, ubar) %>%
    filter(!is.na(timeMid))
}


#' Compute profile diffs binned by stability class
#'
#' Like calc_profile_diffs but uses stability (from aligned data) as
#' the grouping variable instead of or in addition to season/TOD.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param met_ts Tibble. Output of extract_met_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param group_vars Character vector. Default: c("season", "stability")
#' @return Tibble with pairwise diffs grouped by stability
calc_profile_diffs_stability <- function(profile_ts, met_ts, attr_df,
                                          gas = "CH4",
                                          group_vars = c("season", "stability")) {

  if (nrow(profile_ts) == 0 || nrow(met_ts) == 0) return(tibble())

  # Round profile times and join stability
  df <- profile_ts %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    left_join(met_ts %>% select(timeMid, stability) %>% rename(time_round = timeMid),
              by = "time_round") %>%
    filter(!is.na(stability))

  # Now compute diffs using the joined stability as a grouping variable
  height_lookup <- get_height_lookup(attr_df)

  df2 <- df %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    add_season("timeMid")

  height_means <- df2 %>%
    group_by(across(all_of(group_vars)), TowerPosition, height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      sd_conc = sd(concentration, na.rm = TRUE),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  height_means %>%
    inner_join(height_means,
               by = group_vars,
               suffix = c("_upper", "_lower"),
               relationship = "many-to-many") %>%
    filter(height_m_upper > height_m_lower) %>%
    mutate(
      dConc = mean_conc_upper - mean_conc_lower,
      is_adjacent = abs(TowerPosition_upper - TowerPosition_lower) == 1,
      gas = gas
    ) %>%
    rename(
      height_upper = height_m_upper,
      height_lower = height_m_lower,
      position_upper = TowerPosition_upper,
      position_lower = TowerPosition_lower,
      conc_upper = mean_conc_upper,
      conc_lower = mean_conc_lower
    ) %>%
    select(all_of(group_vars),
           position_upper, position_lower,
           height_upper, height_lower,
           conc_upper, conc_lower,
           dConc, is_adjacent,
           n_upper, n_lower,
           sd_conc_upper, sd_conc_lower,
           gas)
}
