# plot_profile_heatmap.R
# Hour-of-day x height heatmap of mean concentration

#' Heatmap of concentration by hour and height
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param months Numeric vector or NULL. Filter to specific months.
#' @return ggplot object
plot_profile_heatmap <- function(profile_ts, attr_df, gas = "CH4",
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

  # Aggregate
  heatmap_data <- df %>%
    group_by(hour, height_m, season) %>%
    summarise(
      conc_mean = mean(concentration, na.rm = TRUE),
      .groups = "drop"
    )

  title_text <- paste0(gas, " Concentration Profile (Hour x Height)")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  ggplot(heatmap_data, aes(x = hour, y = height_m, fill = conc_mean)) +
    geom_tile() +
    scale_x_continuous(breaks = seq(0, 23, 3), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_viridis_c(option = "inferno", name = paste0(gas, "\n(", gas_units(gas), ")")) +
    facet_wrap(~ season) +
    labs(
      x = "Hour of day",
      y = "Height (m)",
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
}
