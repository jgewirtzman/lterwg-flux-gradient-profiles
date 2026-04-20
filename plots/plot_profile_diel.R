# plot_profile_diel.R
# Diel (hour-of-day) concentration profiles colored by tower height

#' Plot diel concentration pattern at each tower height
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param months Numeric vector or NULL. Filter to specific months.
#' @return ggplot object
plot_profile_diel <- function(profile_ts, attr_df, gas = "CH4",
                               site_name = NULL, months = NULL) {

  if (nrow(profile_ts) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  height_lookup <- get_height_lookup(attr_df)

  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    add_season("timeMid") %>%
    mutate(hour = hour(timeMid))

  if (!is.null(months)) {
    df <- df %>% filter(month(timeMid) %in% months)
  }

  # Aggregate by hour x height x season
  diel_summary <- df %>%
    group_by(hour, height_m, season) %>%
    summarise(
      conc_mean = mean(concentration, na.rm = TRUE),
      conc_se = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    ) %>%
    mutate(height_label = paste0(round(height_m, 1), " m"))

  title_text <- paste0(gas, " Diel Concentration by Height")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  ggplot(diel_summary, aes(x = hour, y = conc_mean,
                            color = reorder(height_label, height_m),
                            group = height_label)) +
    geom_ribbon(aes(ymin = conc_mean - conc_se, ymax = conc_mean + conc_se,
                    fill = reorder(height_label, height_m)),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    scale_color_viridis_d(option = "D", end = 0.9, name = "Height") +
    scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
    facet_wrap(~ season, scales = "free_y") +
    labs(
      x = "Hour of day",
      y = paste0(gas, " (", gas_units(gas), ")"),
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}
