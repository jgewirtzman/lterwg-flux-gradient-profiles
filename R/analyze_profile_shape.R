# analyze_profile_shape.R
# Classify and quantify vertical concentration profile shapes

#' Classify the shape of each concentration profile
#'
#' For each timestamp, examines the concentration values across all tower heights
#' and assigns a profile shape category. Uses the long-format profile data joined
#' with heights.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param uniform_threshold Numeric. If range < this value, profile is "uniform".
#'   Default: 2 ppb for CH4, 2 ppm for CO2, 0.5 mmol/mol for H2O.
#' @return Tibble with one row per timestamp: timeMid, profile_shape, plus metrics
classify_profile_shape <- function(profile_ts, attr_df, gas = "CH4",
                                    uniform_threshold = NULL) {

  if (nrow(profile_ts) == 0) return(tibble())

  # Set default threshold based on gas

  if (is.null(uniform_threshold)) {
    uniform_threshold <- switch(gas,
      "CH4" = 2,    # ppb
      "CO2" = 2,    # ppm
      "H2O" = 0.5,  # mmol/mol
      2
    )
  }

  height_lookup <- get_height_lookup(attr_df)

  # Join heights and work per timestamp
  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    group_by(timeMid) %>%
    filter(n() >= 3) %>%  # need at least 3 levels for shape classification
    arrange(height_m, .by_group = TRUE) %>%
    summarise(
      n_levels = n(),
      conc_range = max(concentration) - min(concentration),
      height_of_max = height_m[which.max(concentration)],
      height_of_min = height_m[which.min(concentration)],
      max_height = max(height_m),
      min_height = min(height_m),
      # Gradient: slope of linear fit conc ~ height
      gradient_strength = tryCatch(
        coef(lm(concentration ~ height_m))[2],
        error = function(e) NA_real_
      ),
      # Curvature: quadratic coefficient
      curvature = tryCatch(
        if (n() >= 3) coef(lm(concentration ~ height_m + I(height_m^2)))[3] else NA_real_,
        error = function(e) NA_real_
      ),
      # Is max/min at interior level?
      max_is_interior = height_of_max != max_height & height_of_max != min_height,
      min_is_interior = height_of_min != max_height & height_of_min != min_height,
      # Monotonicity check: are concentrations monotonically ordered with height?
      diffs = list(diff(concentration)),
      .groups = "drop"
    ) %>%
    mutate(
      all_increasing = map_lgl(diffs, ~ all(.x > 0)),
      all_decreasing = map_lgl(diffs, ~ all(.x < 0)),
      # Classify
      profile_shape = case_when(
        conc_range < uniform_threshold         ~ "uniform",
        all_decreasing                         ~ "monotonic_decreasing",
        all_increasing                         ~ "monotonic_increasing",
        max_is_interior                        ~ "mid_canopy_max",
        min_is_interior                        ~ "mid_canopy_min",
        TRUE                                   ~ "complex"
      ),
      profile_shape = factor(profile_shape, levels = c(
        "monotonic_decreasing", "monotonic_increasing",
        "mid_canopy_max", "mid_canopy_min",
        "uniform", "complex"
      ))
    ) %>%
    select(-diffs, -all_increasing, -all_decreasing, -max_is_interior, -min_is_interior)

  df
}


#' Compute detailed profile metrics per timestamp
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @return Tibble with metrics per timestamp
calc_profile_metrics <- function(profile_ts, attr_df, gas = "CH4") {

  if (nrow(profile_ts) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)
  canopy_h <- get_canopy_height(attr_df)

  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    group_by(timeMid) %>%
    filter(n() >= 2) %>%
    arrange(height_m, .by_group = TRUE) %>%
    summarise(
      n_levels = n(),
      conc_range = max(concentration) - min(concentration),
      conc_mean = mean(concentration),
      height_of_max = height_m[which.max(concentration)],
      height_of_min = height_m[which.min(concentration)],
      gradient_strength = tryCatch(
        coef(lm(concentration ~ height_m))[2],
        error = function(e) NA_real_
      ),
      # Below vs above canopy means (if canopy height is known)
      below_canopy_mean = if (!is.na(canopy_h) && any(height_m < canopy_h)) {
        mean(concentration[height_m < canopy_h])
      } else NA_real_,
      above_canopy_mean = if (!is.na(canopy_h) && any(height_m >= canopy_h)) {
        mean(concentration[height_m >= canopy_h])
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      # Crown enrichment: below-canopy minus above-canopy
      below_above_diff = below_canopy_mean - above_canopy_mean,
      gas = gas
    )

  # Add temporal context
  df %>%
    add_season("timeMid") %>%
    add_tod("timeMid") %>%
    mutate(hour = hour(timeMid))
}


#' Summarize profile shape frequencies by condition
#'
#' @param profile_shapes Tibble. Output of classify_profile_shape()
#' @return Tibble with shape frequency by season x TOD
summarize_profile_shapes <- function(profile_shapes) {

  profile_shapes %>%
    add_season("timeMid") %>%
    add_tod("timeMid") %>%
    group_by(season, tod, profile_shape) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(season, tod) %>%
    mutate(
      total = sum(n),
      fraction = n / total
    ) %>%
    ungroup()
}
