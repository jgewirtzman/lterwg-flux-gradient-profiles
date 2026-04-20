# prepare_profiles.R
# Reshape raw 9-min concentration data into profile formats for analysis

#' Prepare long-format time series of concentration at each tower level
#'
#' Takes raw 9-min NEON data (one row per level per timestamp) and returns
#' a clean long-format tibble with quality-filtered concentrations.
#'
#' @param raw_9min List with elements CH4, CO2, H2O (from SITE_9min.RData)
#' @param gas Character. Which gas to process: "CH4", "CO2", or "H2O"
#' @return Tibble with columns: timeMid, timeBgn, timeEnd, TowerPosition,
#'   concentration, vari, numSamp
prepare_profile_timeseries <- function(raw_9min, gas = "CH4") {

  df <- raw_9min[[gas]]
  if (is.null(df) || nrow(df) == 0) {
    warning("No data for gas: ", gas)
    return(tibble())
  }

  df %>%
    # Parse timestamps
    mutate(
      timeBgn = ymd_hms(timeBgn, quiet = TRUE),
      timeEnd = ymd_hms(timeEnd, quiet = TRUE)
    ) %>%
    # Quality filter
    filter(qfFinl == 0, !is.na(mean)) %>%
    # Compute midpoint time
    mutate(timeMid = timeBgn + (timeEnd - timeBgn) / 2) %>%
    # Select and rename
    transmute(
      timeMid,
      timeBgn,
      timeEnd,
      TowerPosition = as.numeric(TowerPosition),
      concentration = mean,
      vari,
      numSamp
    ) %>%
    arrange(timeMid, TowerPosition)
}


#' Prepare wide-format profile data (one row per timestamp, columns per height)
#'
#' Joins with site attributes to get actual measurement heights, then pivots
#' so each row is a single timestamp with concentration at each height.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes with TowerPosition and DistZaxsLvlMeasTow
#' @return Tibble with columns: timeMid, conc_<height>m for each level,
#'   plus hour, month, season, year, doy
prepare_profile_wide <- function(profile_ts, attr_df) {

  if (nrow(profile_ts) == 0) return(tibble())

  # Build tower position -> height lookup from attributes
  height_lookup <- attr_df %>%
    select(TowerPosition, DistZaxsLvlMeasTow) %>%
    distinct() %>%
    mutate(
      TowerPosition = as.numeric(TowerPosition),
      height_m = as.numeric(DistZaxsLvlMeasTow)
    ) %>%
    filter(!is.na(height_m)) %>%
    select(TowerPosition, height_m)

  # Join heights
  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition")

  if (nrow(df) == 0) {
    warning("No matching tower positions between profile data and attributes")
    return(tibble())
  }

  # Create height label for column names
  df <- df %>%
    mutate(height_label = paste0("conc_", round(height_m, 1), "m"))

  # Pivot wider: one row per timestamp, one column per height
  df_wide <- df %>%
    select(timeMid, height_label, concentration) %>%
    pivot_wider(
      names_from = height_label,
      values_from = concentration,
      values_fn = mean  # handle rare duplicates
    )

  # Add temporal columns
  df_wide <- df_wide %>%
    mutate(
      hour = hour(timeMid),
      month = month(timeMid),
      year = year(timeMid),
      doy = yday(timeMid)
    ) %>%
    add_season("timeMid") %>%
    add_tod("timeMid")

  df_wide
}


#' Get the height lookup table for a site
#'
#' @param attr_df Data frame. Site attributes
#' @return Tibble with TowerPosition and height_m, sorted by height
get_height_lookup <- function(attr_df) {
  attr_df %>%
    select(TowerPosition, DistZaxsLvlMeasTow) %>%
    distinct() %>%
    mutate(
      TowerPosition = as.numeric(TowerPosition),
      height_m = as.numeric(DistZaxsLvlMeasTow)
    ) %>%
    filter(!is.na(height_m)) %>%
    select(TowerPosition, height_m) %>%
    arrange(height_m)
}


#' Join met/stability time series to profile data
#'
#' Rounds profile timestamps to 30-min and joins PAR, stability, ustar
#' from the met time series extracted from aligned data.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param met_ts Tibble. Output of extract_met_timeseries()
#' @return profile_ts with added columns: PAR, zoL, ustar, stability, tod (PAR-based)
join_met_to_profiles <- function(profile_ts, met_ts) {

  if (is.null(met_ts) || nrow(met_ts) == 0) return(profile_ts)

  profile_ts %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    left_join(
      met_ts %>% rename(time_round = timeMid),
      by = "time_round"
    ) %>%
    select(-time_round) %>%
    # Add PAR-based day/night if PAR is available
    {if ("PAR" %in% names(.)) add_tod_par(., "PAR") else add_tod(., "timeMid")}
}


#' Get canopy height for a site from attributes
#'
#' @param attr_df Data frame. Site attributes
#' @return Numeric. Canopy height in meters (median of DistZaxsCnpy)
get_canopy_height <- function(attr_df) {
  if ("DistZaxsCnpy" %in% names(attr_df)) {
    h <- as.numeric(attr_df$DistZaxsCnpy)
    return(median(h, na.rm = TRUE))
  }
  NA_real_
}
