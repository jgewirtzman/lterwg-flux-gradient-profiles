# plot_multi_site_summary.R
# Cross-site summary visualizations

#' Tile plot summarizing profile characteristics across all sites
#'
#' @param all_metrics Named list (by site) of profile metrics tibbles
#'   (output of calc_profile_metrics)
#' @param gas Character. Gas name
#' @return ggplot object
plot_multi_site_profile_summary <- function(all_metrics, gas = "CH4") {

  # Combine metrics across sites
  combined <- bind_rows(
    lapply(names(all_metrics), function(site) {
      m <- all_metrics[[site]][[gas]]
      if (!is.null(m) && nrow(m) > 0) {
        m %>% mutate(site = site)
      } else NULL
    })
  )

  if (nrow(combined) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  # Summary per site x season x TOD
  summary <- combined %>%
    group_by(site, season, tod) %>%
    summarise(
      mean_gradient = mean(gradient_strength, na.rm = TRUE),
      mean_range = mean(conc_range, na.rm = TRUE),
      mean_below_above_diff = mean(below_above_diff, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )

  ggplot(summary, aes(x = site, y = season, fill = mean_below_above_diff)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = gas_colors(gas)$low, mid = "white", high = gas_colors(gas)$high,
      midpoint = 0,
      name = paste0("Below - Above\nCanopy ", gas, "\n(", gas_units(gas), ")")
    ) +
    facet_wrap(~ tod) +
    labs(
      x = NULL,
      y = NULL,
      title = paste0("Cross-Site ", gas, " Below vs Above Canopy Concentration Difference")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    )
}


#' Dot plot of profile shape frequencies across sites
#'
#' @param all_shapes Named list (by site) of profile shape classification tibbles
#' @param gas Character. Gas name
#' @return ggplot object
plot_multi_site_shape_frequency <- function(all_shapes, gas = "CH4") {

  combined <- bind_rows(
    lapply(names(all_shapes), function(site) {
      s <- all_shapes[[site]][[gas]]
      if (!is.null(s) && nrow(s) > 0) {
        s %>% mutate(site = site)
      } else NULL
    })
  )

  if (nrow(combined) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  # Compute shape frequencies per site x season
  freq <- combined %>%
    add_season("timeMid") %>%
    add_tod("timeMid") %>%
    group_by(site, season, tod, profile_shape) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(site, season, tod) %>%
    mutate(fraction = n / sum(n)) %>%
    ungroup()

  ggplot(freq, aes(x = site, y = fraction, fill = profile_shape)) +
    geom_col(position = "stack") +
    facet_grid(tod ~ season) +
    scale_fill_brewer(palette = "Set2", name = "Profile Shape") +
    labs(
      x = NULL,
      y = "Fraction of profiles",
      title = paste0("Cross-Site ", gas, " Profile Shape Distribution")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      strip.text = element_text(face = "bold")
    )
}
