# compare_tracers.R
# Multi-tracer comparison: CH4 vs CO2 profile shapes
#
# Comparing the vertical structure of CH4 and CO2 profiles helps distinguish:
#  - Shared soil source + trapping (both gases accumulate near ground at night)
#  - Independent canopy-level CH4 process (CH4 shows mid-canopy enrichment
#    that CO2 does not, or vice versa)

#' Compare CH4 and CO2 profile shapes at each timestamp
#'
#' @param profile_ts_ch4 Tibble. Output of prepare_profile_timeseries() for CH4
#' @param profile_ts_co2 Tibble. Output of prepare_profile_timeseries() for CO2
#' @param attr_df Data frame. Site attributes
#' @return Tibble with per-timestamp comparison metrics
compare_tracer_profiles <- function(profile_ts_ch4, profile_ts_co2, attr_df) {

  if (nrow(profile_ts_ch4) == 0 || nrow(profile_ts_co2) == 0) return(tibble())

  height_lookup <- get_height_lookup(attr_df)

  # Prepare each gas: join heights, compute per-timestamp gradient
  prep_gas <- function(prof_ts, gas_name) {
    prof_ts %>%
      inner_join(height_lookup, by = "TowerPosition") %>%
      group_by(timeMid) %>%
      filter(n() >= 3) %>%
      arrange(height_m, .by_group = TRUE) %>%
      summarise(
        n_levels = n(),
        gradient = tryCatch(coef(lm(concentration ~ height_m))[2], error = function(e) NA_real_),
        range = max(concentration) - min(concentration),
        height_of_max = height_m[which.max(concentration)],
        # Normalized gradient (dimensionless shape)
        norm_gradient = if (range > 0) gradient / range else 0,
        # Full concentration vector for correlation
        conc_vec = list(concentration),
        height_vec = list(height_m),
        .groups = "drop"
      ) %>%
      rename_with(~ paste0(.x, "_", gas_name), -timeMid)
  }

  ch4 <- prep_gas(profile_ts_ch4, "CH4")
  co2 <- prep_gas(profile_ts_co2, "CO2")

  # Join by closest timestamps (within 5 minutes)
  # Since CH4 and CO2 are measured at different times in the Picarro cycle,
  # we match by nearest timeMid within a tolerance
  combined <- ch4 %>%
    inner_join(co2, by = "timeMid", relationship = "many-to-many") # exact match first

  if (nrow(combined) == 0) {
    # Try fuzzy join by rounding to nearest 30 min
    ch4_r <- ch4 %>% mutate(time_round = round_date(timeMid, "30 minutes"))
    co2_r <- co2 %>% mutate(time_round = round_date(timeMid, "30 minutes"))
    combined <- ch4_r %>%
      inner_join(co2_r, by = "time_round", suffix = c("", ".co2")) %>%
      select(-time_round)
  }

  if (nrow(combined) == 0) return(tibble())

  combined %>%
    mutate(
      # Gradient correlation: do CH4 and CO2 have similar vertical structure?
      gradient_correlation = map2_dbl(
        conc_vec_CH4, conc_vec_CO2,
        function(ch4_c, co2_c) {
          if (length(ch4_c) == length(co2_c) && length(ch4_c) >= 3) {
            cor(ch4_c, co2_c, use = "complete.obs")
          } else NA_real_
        }
      ),
      # Do both gases have their max at the same height?
      max_height_agreement = height_of_max_CH4 == height_of_max_CO2,
      # Divergence: how different are the normalized gradients?
      gradient_divergence = abs(norm_gradient_CH4 - norm_gradient_CO2)
    ) %>%
    select(
      timeMid,
      gradient_CH4 = gradient_CH4,
      gradient_CO2 = gradient_CO2,
      range_CH4, range_CO2,
      height_of_max_CH4, height_of_max_CO2,
      gradient_correlation,
      max_height_agreement,
      gradient_divergence,
      norm_gradient_CH4, norm_gradient_CO2
    ) %>%
    add_season("timeMid") %>%
    add_tod("timeMid") %>%
    mutate(hour = hour(timeMid))
}


#' Summarize tracer comparison by conditions
#'
#' @param tracer_comp Tibble. Output of compare_tracer_profiles()
#' @return Tibble with summary stats by season x TOD
summarize_tracer_comparison <- function(tracer_comp) {

  if (nrow(tracer_comp) == 0) return(tibble())

  tracer_comp %>%
    group_by(season, tod) %>%
    summarise(
      n = n(),
      mean_gradient_corr = mean(gradient_correlation, na.rm = TRUE),
      sd_gradient_corr = sd(gradient_correlation, na.rm = TRUE),
      frac_max_agreement = mean(max_height_agreement, na.rm = TRUE),
      mean_divergence = mean(gradient_divergence, na.rm = TRUE),
      .groups = "drop"
    )
}
