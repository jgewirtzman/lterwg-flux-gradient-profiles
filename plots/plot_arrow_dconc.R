# plot_arrow_dconc.R
# Arrow visualization of concentration differences between tower levels
#
# Arrows point from high concentration to low concentration.
# Arrow size and color map to dConc magnitude and sign.
# Faceted by season and time of day (or stability).
#
# Uses pairwise differences computed from raw 9-min profile MEANS
# (not the aligned paired-measurement dConc, which attenuates the signal).

#' Arrow plot of concentration differences between tower levels
#'
#' @param profile_diffs Tibble. Output of calc_profile_diffs() or
#'   calc_profile_diffs_stability(). Must have columns: dConc, height_upper,
#'   height_lower, position_upper, position_lower, is_adjacent, plus
#'   faceting variables (e.g., season, tod or stability).
#' @param gas Character. "CH4", "CO2", or "H2O"
#' @param site_name Character. Site code for plot title
#' @param facet_by Character vector. Variables to facet by.
#'   Default: c("season", "tod") -> facet_grid(tod ~ season)
#' @param show_adjacent_panel Logical. Show highlighted adjacent-level panel on right?
#' @return ggplot object
plot_arrow_dconc <- function(profile_diffs,
                              gas = "CH4",
                              site_name = NULL,
                              facet_by = c("season", "tod"),
                              show_adjacent_panel = TRUE) {

  conc_data <- profile_diffs
  if (is.null(conc_data) || nrow(conc_data) == 0) {
    warning("No data to plot")
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  # Arrow direction: from high concentration to low
  conc_data <- conc_data %>%
    mutate(
      arrow_y    = ifelse(dConc > 0, height_upper, height_lower),
      arrow_yend = ifelse(dConc > 0, height_lower, height_upper)
    )

  # X-axis positioning: group by lower position, offset within group
  conc_data <- conc_data %>%
    group_by(position_lower, across(any_of(facet_by))) %>%
    mutate(
      x = position_lower + (row_number() - (n() + 1) / 2) * 0.15
    ) %>%
    ungroup()

  # Axis reference values
  all_heights <- sort(unique(c(conc_data$height_upper, conc_data$height_lower)))
  all_positions <- sort(unique(c(conc_data$position_upper, conc_data$position_lower)))

  # Gas-specific colors
  colors <- gas_colors(gas)

  # Build the plot
  p <- ggplot()

  # Main arrows: all level pairs
  p <- p +
    geom_segment(
      data = conc_data,
      aes(x = x, y = arrow_y,
          xend = x, yend = arrow_yend,
          color = dConc, linewidth = abs(dConc)),
      arrow = arrow(type = "closed", length = unit(0.25, "cm"), angle = 25),
      lineend = "butt"
    )

  # Adjacent-level highlighted panel on right
  if (show_adjacent_panel) {
    adj_data <- conc_data %>% filter(is_adjacent)
    if (nrow(adj_data) > 0) {
      x_adj <- max(conc_data$x, na.rm = TRUE) + 1

      # Black border arrows
      p <- p +
        geom_segment(
          data = adj_data,
          aes(x = x_adj, y = arrow_y,
              xend = x_adj, yend = arrow_yend,
              linewidth = abs(dConc) * 1.3),
          color = "black",
          arrow = arrow(type = "closed", length = unit(0.25, "cm"), angle = 25),
          lineend = "butt"
        ) +
        # Colored interior
        geom_segment(
          data = adj_data,
          aes(x = x_adj, y = arrow_y,
              xend = x_adj, yend = arrow_yend,
              color = dConc, linewidth = abs(dConc)),
          arrow = arrow(type = "closed", length = unit(0.25, "cm"), angle = 25),
          lineend = "butt"
        )
    }
  }

  # Scales
  p <- p +
    scale_color_gradient2(
      low = colors$low, mid = colors$mid, high = colors$high,
      midpoint = 0,
      name = paste0("\u0394", gas, " (", gas_units(gas), ")")
    ) +
    scale_linewidth_continuous(range = c(0.3, 4), guide = "none") +
    scale_y_continuous(
      breaks = all_heights,
      labels = paste(all_heights, "m")
    )

  # X-axis labels
  if (show_adjacent_panel) {
    x_adj <- max(conc_data$x, na.rm = TRUE) + 1
    p <- p + scale_x_continuous(
      breaks = c(all_positions, x_adj),
      labels = c(paste("Level", all_positions), "Adjacent\nOnly")
    )
  } else {
    p <- p + scale_x_continuous(
      breaks = all_positions,
      labels = paste("Level", all_positions)
    )
  }

  # Labels and theme
  title_text <- paste0(gas, " Concentration Differences Between Tower Levels")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  p <- p +
    labs(
      x = NULL,
      y = "Height (m)",
      title = title_text,
      subtitle = "Arrows point from high to low concentration"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
      strip.text = element_text(face = "bold")
    )

  # Faceting
  if (length(facet_by) == 2) {
    p <- p + facet_grid(reformulate(facet_by[1], facet_by[2]))
  } else if (length(facet_by) == 1) {
    p <- p + facet_wrap(reformulate(facet_by))
  }

  p
}
