# classify_sites.R
# Classify sites by CH4 profile shape in JJA Night:
#   - Soil sink vs soil source vs neutral (based on ground-level anomaly)
#   - Monotonic vs mid-canopy feature (non-monotonic structure within canopy)

#' Compute site classification based on JJA Night anomaly profile
#'
#' @param all_profiles Named list (by site) of profile_ts tibbles
#' @param all_met Named list (by site) of met time series
#' @param all_data Named list (by site) of loaded site data (needs $attr)
#' @param gas Character. Gas to classify. Default: "CH4"
#' @param ground_threshold Numeric. Abs anomaly below which to call "neutral".
#'   Default: 0.002 (ppb for CH4; for CO2 you'd want larger).
#' @param midcanopy_threshold Numeric. Minimum deviation at interior level
#'   (relative to linear interpolation between endpoints) to call "mid-canopy feature".
#'   Default: 0.003 (ppb for CH4).
#' @return Tibble with columns: site, ground_anomaly, soil_class,
#'   midcanopy_deviation, structure_class, combined_class
classify_sites <- function(all_profiles, all_met, all_data,
                            gas = "CH4",
                            ground_threshold = 0.002,
                            midcanopy_threshold = 0.003) {

  rows <- list()

  for (site in names(all_profiles)) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]

    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

    df <- prof_ts %>%
      inner_join(height_lookup, by = "TowerPosition") %>%
      add_season("timeMid")

    # PAR-based or fixed TOD
    if (!is.null(met) && "PAR" %in% names(met)) {
      df <- df %>%
        mutate(time_round = round_date(timeMid, "30 minutes")) %>%
        left_join(met %>% select(timeMid, PAR) %>% rename(time_round = timeMid),
                  by = "time_round") %>%
        add_tod_par("PAR") %>%
        select(-time_round, -PAR)
    } else {
      df <- df %>% add_tod("timeMid")
    }

    # JJA Night only (strongest structure)
    jja_night <- df %>% filter(season == "JJA", tod == "Night")
    if (nrow(jja_night) == 0) next

    # Mean conc per height
    height_means <- jja_night %>%
      group_by(height_m) %>%
      summarise(mean_conc = mean(concentration, na.rm = TRUE), .groups = "drop") %>%
      arrange(height_m)

    if (nrow(height_means) < 3) next

    # Background = top-of-tower mean
    background <- height_means$mean_conc[height_means$height_m == top_height]
    height_means <- height_means %>% mutate(anomaly = mean_conc - background)

    # Ground anomaly (lowest level)
    ground_anom <- height_means$anomaly[1]

    # Soil classification
    soil_class <- case_when(
      ground_anom >  ground_threshold ~ "soil source",
      ground_anom < -ground_threshold ~ "soil sink",
      TRUE ~ "neutral"
    )

    # Mid-canopy feature: for each interior level, compute deviation from
    # linear interpolation between ground and top. If any interior level
    # deviates by more than the threshold, flag as mid-canopy feature.
    h <- height_means$height_m
    a <- height_means$anomaly
    n <- length(h)

    # Linear interpolation between ground (i=1) and top (i=n)
    expected <- a[1] + (a[n] - a[1]) * (h - h[1]) / (h[n] - h[1])
    deviations <- a - expected

    # Consider only interior levels (exclude endpoints)
    interior_dev <- deviations[2:(n - 1)]
    max_abs_dev <- max(abs(interior_dev))

    structure_class <- if (max_abs_dev > midcanopy_threshold) {
      "mid-canopy feature"
    } else {
      "monotonic"
    }

    combined_class <- paste(soil_class, structure_class, sep = " / ")

    rows[[site]] <- tibble(
      site = site,
      n_levels = n,
      ground_anomaly = ground_anom,
      soil_class = soil_class,
      midcanopy_deviation = max_abs_dev,
      structure_class = structure_class,
      combined_class = combined_class
    )
  }

  classifications <- bind_rows(rows)

  # Order the combined class factor
  class_levels <- c(
    "soil source / mid-canopy feature",
    "soil source / monotonic",
    "soil sink / mid-canopy feature",
    "soil sink / monotonic",
    "neutral / mid-canopy feature",
    "neutral / monotonic"
  )
  classifications %>%
    mutate(combined_class = factor(combined_class,
                                     levels = intersect(class_levels,
                                                         unique(combined_class))))
}
