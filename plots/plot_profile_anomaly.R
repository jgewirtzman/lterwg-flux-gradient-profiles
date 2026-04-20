# plot_profile_anomaly.R
# Concentration anomaly relative to top-of-tower (atmospheric background)
# Red = enriched above background (source signal)
# Blue = depleted below background (sink signal)
# Day and Night stacked on the same panel, distinguished by shape/alpha

#' Anomaly profile relative to top-of-tower background
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param met_ts Tibble or NULL. If provided, uses PAR-based day/night.
#' @param facet_var Character. Variable for faceting panels. Default: "season"
#' @return ggplot object
plot_profile_anomaly <- function(profile_ts, attr_df, gas = "CH4",
                                  site_name = NULL, met_ts = NULL,
                                  facet_var = "season") {

  if (nrow(profile_ts) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  height_lookup <- get_height_lookup(attr_df)
  canopy_h <- get_canopy_height(attr_df)
  top_height <- max(height_lookup$height_m)

  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    add_season("timeMid")

  # PAR-based or fixed day/night
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

  # Mean concentration per height per condition
  height_means <- df %>%
    group_by(across(all_of(c(facet_var, "tod"))), height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  # Get top-of-tower mean per condition (the background reference)
  top_means <- height_means %>%
    filter(height_m == top_height) %>%
    select(all_of(c(facet_var, "tod")), background = mean_conc)

  # Compute anomaly
  anomaly <- height_means %>%
    left_join(top_means, by = c(facet_var, "tod")) %>%
    mutate(anomaly = mean_conc - background)

  # Slight vertical dodge so day/night bars don't overlap
  anomaly <- anomaly %>%
    mutate(
      height_dodge = height_m + ifelse(tod == "Day", 0.5, -0.5)
    )

  title_text <- paste0(gas, " Anomaly Relative to Tower Top")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  colors <- gas_colors(gas)

  p <- ggplot(anomaly, aes(x = anomaly, y = height_dodge)) +
    # Vertical reference line at zero
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
    # Bars
    geom_segment(aes(x = 0, xend = anomaly, y = height_dodge, yend = height_dodge,
                     color = anomaly, linetype = tod),
                 linewidth = 2) +
    # Points
    geom_point(aes(color = anomaly, shape = tod), size = 3) +
    # Error bars
    geom_errorbarh(aes(xmin = anomaly - se_conc, xmax = anomaly + se_conc),
                   height = 0.6, linewidth = 0.3, color = "grey40") +
    # Scales
    scale_color_gradient2(
      low = colors$low, mid = "grey80", high = colors$high,
      midpoint = 0,
      name = paste0("\u0394", gas)
    ) +
    scale_shape_manual(values = c("Day" = 16, "Night" = 17), name = NULL) +
    scale_linetype_manual(values = c("Day" = "solid", "Night" = "dashed"), name = NULL) +
    scale_y_continuous(
      breaks = sort(unique(height_lookup$height_m)),
      labels = paste(round(sort(unique(height_lookup$height_m)), 1), "m")
    ) +
    labs(
      x = paste0(gas, " anomaly from tower top (", gas_units(gas), ")"),
      y = "Height (m)",
      title = title_text,
      subtitle = paste0("Circles = Day (PAR > 10) | Triangles = Night")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      legend.position = "right"
    )

  # Canopy height reference
  if (!is.na(canopy_h)) {
    p <- p + geom_hline(yintercept = canopy_h, linetype = "dashed",
                         color = "forestgreen", linewidth = 0.5) +
      annotate("text", x = Inf, y = canopy_h, label = "canopy ",
               hjust = 1, vjust = -0.5, size = 3, color = "forestgreen")
  }

  # Facet by season
  p <- p + facet_wrap(reformulate(facet_var), nrow = 1)

  p
}
