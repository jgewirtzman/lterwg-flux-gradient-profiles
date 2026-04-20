# plot_canopy_structure.R
# Visualize canopy height and tower sampling structure across sites

#' Plot canopy structure: canopy height vs tower levels per site
#'
#' Each site shown as a vertical "tower" with measurement levels as dots,
#' canopy height as a dashed green line, colored by canopy_class.
#'
#' @param canopy_metrics Tibble. Output of compute_canopy_metrics_all()
#' @param all_data Named list (by site) of site data with $attr
#' @return ggplot object
plot_canopy_structure <- function(canopy_metrics, all_data) {

  # Build long-format: one row per site x level
  level_rows <- list()
  for (site in canopy_metrics$site) {
    attr_df <- all_data[[site]]$attr
    if (is.null(attr_df)) next
    heights <- as.numeric(attr_df$DistZaxsLvlMeasTow)
    level_rows[[site]] <- tibble(
      site   = site,
      height = heights,
      level  = seq_along(heights)
    )
  }
  levels_df <- bind_rows(level_rows) %>%
    left_join(canopy_metrics %>% select(site, canopy_height, canopy_class),
              by = "site")

  # Order sites by canopy class then by canopy height
  site_order <- canopy_metrics %>%
    arrange(canopy_class, canopy_height) %>%
    pull(site)
  levels_df$site <- factor(levels_df$site, levels = site_order)
  canopy_metrics$site <- factor(canopy_metrics$site, levels = site_order)

  class_colors <- c(
    "No canopy (<1 m)"        = "#E8C547",
    "Short canopy (1-5 m)"    = "#A8D95B",
    "Medium canopy (5-15 m)"  = "#5BBA6F",
    "Tall canopy (15-30 m)"   = "#2D8B3A",
    "Very tall canopy (>30 m)"= "#0F5A1E"
  )

  ggplot() +
    # Canopy height line per site
    geom_segment(data = canopy_metrics,
                 aes(x = site, xend = site,
                     y = 0, yend = canopy_height,
                     color = canopy_class),
                 linewidth = 3, alpha = 0.3) +
    # Measurement levels as dots
    geom_point(data = levels_df,
               aes(x = site, y = height, color = canopy_class),
               size = 2) +
    # Canopy height as horizontal tick
    geom_segment(data = canopy_metrics,
                 aes(x = as.numeric(site) - 0.3,
                     xend = as.numeric(site) + 0.3,
                     y = canopy_height, yend = canopy_height),
                 color = "forestgreen", linewidth = 0.6) +
    scale_color_manual(values = class_colors, name = "Canopy class") +
    labs(
      x = NULL, y = "Height above ground (m)",
      title = "NEON Tower Canopy Structure",
      subtitle = "Dots = measurement levels; green tick = canopy top; bar = canopy depth"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "top"
    )
}
