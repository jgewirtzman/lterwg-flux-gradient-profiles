# plot_profile_dotline.R
# Option C: Dot-and-line vertical profile with error bars
# Clean, publication-ready format with canopy annotation

#' Dot-and-line vertical concentration profile
#'
#' Overlays multiple conditions (e.g., Day/Night) on the same panel,
#' distinguished by color, with one panel per season.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param color_var Character. Variable mapped to color. Default: "tod"
#' @param facet_var Character. Variable for faceting. Default: "season"
#' @return ggplot object
plot_profile_dotline <- function(profile_ts, attr_df, gas = "CH4",
                                  site_name = NULL,
                                  color_var = "tod",
                                  facet_var = "season") {

  if (nrow(profile_ts) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  height_lookup <- get_height_lookup(attr_df)
  canopy_h <- get_canopy_height(attr_df)

  df <- profile_ts %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    add_season("timeMid") %>%
    add_tod("timeMid")

  # Mean + SE per height per condition
  summary <- df %>%
    group_by(across(all_of(c(color_var, facet_var))), height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  title_text <- paste0(gas, " Vertical Profile")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  # Slight vertical offset to avoid overlap
  n_colors <- length(unique(summary[[color_var]]))
  summary <- summary %>%
    group_by(across(all_of(facet_var)), height_m) %>%
    mutate(
      height_dodge = height_m + (as.numeric(factor(.data[[color_var]])) - (n_colors + 1) / 2) * 0.4
    ) %>%
    ungroup()

  p <- ggplot(summary, aes(x = mean_conc, y = height_dodge,
                            color = .data[[color_var]])) +
    # Lines connecting heights
    geom_path(aes(group = .data[[color_var]]), linewidth = 0.6) +
    # Points with error bars
    geom_pointrange(aes(xmin = mean_conc - se_conc,
                        xmax = mean_conc + se_conc),
                    size = 0.5, linewidth = 0.4, fatten = 3) +
    # Color
    scale_color_manual(
      values = c("Day" = "#E6A817", "Night" = "#2C3E50",
                 "unstable" = "#D35400", "neutral" = "#7F8C8D", "stable" = "#2C3E50"),
      name = NULL
    ) +
    # Y-axis: use original heights for labels
    scale_y_continuous(
      breaks = sort(unique(height_lookup$height_m)),
      labels = paste(sort(unique(round(height_lookup$height_m, 1))), "m")
    ) +
    labs(
      x = paste0(gas, " (", gas_units(gas), ")"),
      y = "Height (m)",
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      legend.position = "top"
    )

  # Canopy height
  if (!is.na(canopy_h)) {
    p <- p + geom_hline(yintercept = canopy_h, linetype = "dashed",
                         color = "forestgreen", linewidth = 0.5) +
      annotate("text", x = Inf, y = canopy_h, label = "canopy ",
               hjust = 1, vjust = -0.5, size = 3, color = "forestgreen")
  }

  # Facet
  p <- p + facet_wrap(reformulate(facet_var), nrow = 1)

  p
}
