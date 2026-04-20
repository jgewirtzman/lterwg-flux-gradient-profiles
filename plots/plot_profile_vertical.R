# plot_profile_vertical.R
# Option A: Vertical concentration profile plots
# y = height, x = concentration, one curve per condition bin

#' Vertical concentration profile by season x TOD
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param group_vars Character vector. Faceting variables. Default: c("season", "tod")
#' @return ggplot object
plot_profile_vertical <- function(profile_ts, attr_df, gas = "CH4",
                                   site_name = NULL,
                                   group_vars = c("season", "tod")) {

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
    group_by(across(all_of(group_vars)), height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  title_text <- paste0(gas, " Vertical Concentration Profile")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  p <- ggplot(summary, aes(x = mean_conc, y = height_m)) +
    geom_path(linewidth = 0.8, color = "grey30") +
    geom_pointrange(aes(xmin = mean_conc - se_conc,
                        xmax = mean_conc + se_conc),
                    size = 0.4, linewidth = 0.5, color = "grey20") +
    labs(
      x = paste0(gas, " (", gas_units(gas), ")"),
      y = "Height (m)",
      title = title_text
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
    )

  # Add canopy height reference line
  if (!is.na(canopy_h)) {
    p <- p + geom_hline(yintercept = canopy_h, linetype = "dashed",
                         color = "forestgreen", linewidth = 0.5) +
      annotate("text", x = -Inf, y = canopy_h, label = " canopy",
               hjust = 0, vjust = -0.5, size = 3, color = "forestgreen")
  }

  # Facet
  if (length(group_vars) == 2) {
    p <- p + facet_grid(reformulate(group_vars[1], group_vars[2]),
                        scales = "free_x")
  } else if (length(group_vars) == 1) {
    p <- p + facet_wrap(reformulate(group_vars), scales = "free_x")
  }

  p
}
