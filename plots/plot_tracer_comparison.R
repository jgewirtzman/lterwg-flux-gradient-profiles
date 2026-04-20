# plot_tracer_comparison.R
# Visualize CH4 vs CO2 profile shape comparison

#' Scatter plot of CH4 vs CO2 gradient strength
#'
#' Points on the 1:1 line suggest shared forcing (e.g., both driven by
#' soil source + stability). Deviations suggest independent processes.
#'
#' @param tracer_comp Tibble. Output of compare_tracer_profiles()
#' @param site_name Character. Site code for title
#' @return ggplot object
plot_tracer_gradient_scatter <- function(tracer_comp, site_name = NULL) {

  if (nrow(tracer_comp) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  # Use normalized gradients for comparison
  df <- tracer_comp %>% filter(!is.na(norm_gradient_CH4), !is.na(norm_gradient_CO2))

  title_text <- "CH4 vs CO2 Normalized Gradient"
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  ggplot(df, aes(x = norm_gradient_CO2, y = norm_gradient_CH4, color = tod)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(alpha = 0.3, size = 0.8) +
    scale_color_manual(values = c("Day" = "#E6A817", "Night" = "#2C3E50"),
                       name = NULL) +
    facet_wrap(~ season) +
    labs(
      x = expression(CO[2] ~ "normalized gradient"),
      y = expression(CH[4] ~ "normalized gradient"),
      title = title_text,
      subtitle = "Dashed line = 1:1 (identical profile shapes)"
    ) +
    theme_minimal(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))
}


#' Heatmap of gradient correlation between CH4 and CO2 by hour x season
#'
#' @param tracer_comp Tibble. Output of compare_tracer_profiles()
#' @param site_name Character. Site code for title
#' @return ggplot object
plot_tracer_correlation_heatmap <- function(tracer_comp, site_name = NULL) {

  if (nrow(tracer_comp) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  summary <- tracer_comp %>%
    group_by(hour, season) %>%
    summarise(
      mean_corr = mean(gradient_correlation, na.rm = TRUE),
      n = sum(!is.na(gradient_correlation)),
      .groups = "drop"
    )

  title_text <- "CH4-CO2 Profile Correlation by Hour and Season"
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  ggplot(summary, aes(x = hour, y = season, fill = mean_corr)) +
    geom_tile() +
    scale_x_continuous(breaks = seq(0, 23, 3), expand = c(0, 0)) +
    scale_fill_gradient2(
      low = "#B22222", mid = "white", high = "#4682B4",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(
      x = "Hour of day",
      y = NULL,
      title = title_text,
      subtitle = "High correlation = shared forcing; low = independent processes"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank())
}
