# utils.R
# Shared helper functions for height-resolved concentration analysis

#' Plotting scale factor for per-gas flux display.
#'
#' Internally we compute all CH4 fluxes in µmol/m²/s, but typical values
#' are 10^-3 to 10^-2 µmol/m²/s — hard to read. Multiply by 1000 to display
#' in nmol/m²/s, matching CH4 flux literature. CO2 stays in µmol/m²/s,
#' H2O in mmol/m²/s.
flux_scale_for_gas <- function(gas) {
  switch(gas, "CH4" = 1000, "CO2" = 1, "H2O" = 1, 1)
}

#' Units label for per-gas flux plots.
flux_units_for_gas <- function(gas) {
  switch(gas,
         "CH4" = "nmol/m\u00b2/s",
         "CO2" = "\u00b5mol/m\u00b2/s",
         "H2O" = "mmol/m\u00b2/s",
         "")
}

#' Add season column (DJF/MAM/JJA/SON) based on month
add_season <- function(df, time_col = "timeMid") {
  df %>%
    mutate(season = case_when(
      month(.data[[time_col]]) %in% c(12, 1, 2)  ~ "DJF",
      month(.data[[time_col]]) %in% c(3, 4, 5)   ~ "MAM",
      month(.data[[time_col]]) %in% c(6, 7, 8)   ~ "JJA",
      month(.data[[time_col]]) %in% c(9, 10, 11)  ~ "SON"
    )) %>%
    mutate(season = factor(season, levels = c("DJF", "MAM", "JJA", "SON")))
}

#' Add time-of-day column (Day/Night) using fixed hours (fallback)
add_tod <- function(df, time_col = "timeMid", day_start = 6, day_end = 18) {
  df %>%
    mutate(tod = ifelse(
      hour(.data[[time_col]]) >= day_start & hour(.data[[time_col]]) < day_end,
      "Day", "Night"
    )) %>%
    mutate(tod = factor(tod, levels = c("Day", "Night")))
}

#' Add time-of-day column using PAR (radiation-based)
#' Requires a PAR time series joined to the data.
#' @param df Data frame with a PAR column
#' @param par_col Character. Column name for PAR
#' @param threshold Numeric. PAR threshold for day/night (µmol/m²/s). Default: 10
add_tod_par <- function(df, par_col = "PAR", threshold = 10) {
  df %>%
    mutate(tod = ifelse(
      !is.na(.data[[par_col]]) & .data[[par_col]] > threshold,
      "Day", "Night"
    )) %>%
    mutate(tod = factor(tod, levels = c("Day", "Night")))
}

#' Add atmospheric stability class based on z/L
add_stability_class <- function(df, zol_col = "zoL") {
  df %>%
    mutate(stability = case_when(
      .data[[zol_col]] < -0.1  ~ "unstable",
      .data[[zol_col]] > 0.1   ~ "stable",
      TRUE                      ~ "neutral"
    )) %>%
    mutate(stability = factor(stability, levels = c("unstable", "neutral", "stable")))
}

#' Logical vector identifying non-outliers using IQR method
remove_outliers_iqr <- function(x, multiplier = 3) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  x >= (q1 - multiplier * iqr) & x <= (q3 + multiplier * iqr) & !is.na(x)
}

#' Gas-specific color scale parameters for diverging gradient2 plots
#' Returns named list: low, mid, high
gas_colors <- function(gas) {
  switch(gas,
    "CH4" = list(low = "#4682B4", mid = "white", high = "#B22222"),
    "CO2" = list(low = "forestgreen", mid = "white", high = "chocolate4"),
    "H2O" = list(low = "goldenrod", mid = "white", high = "#4682B4"),
    stop("Unknown gas: ", gas)
  )
}

#' Gas-specific units string
gas_units <- function(gas) {
  switch(gas,
    "CH4" = "ppb",
    "CO2" = "ppm",
    "H2O" = "mmol/mol",
    stop("Unknown gas: ", gas)
  )
}

#' Gas-specific display name
gas_label <- function(gas) {
  switch(gas,
    "CH4" = expression(CH[4]),
    "CO2" = expression(CO[2]),
    "H2O" = expression(H[2]*O),
    gas
  )
}
