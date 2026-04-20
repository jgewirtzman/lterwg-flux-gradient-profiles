# plot_anomaly_dotline.R
# Hybrid: anomaly from tower top (above-canopy atmosphere) shown as
# connected points with SE error bars, day/night overlaid by color,
# faceted by season. Combines the anomaly stat with the clean dotline style.

#' Anomaly-from-atmosphere profile as connected points with day/night overlay
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @param met_ts Tibble or NULL. If provided, uses PAR-based day/night.
#' @param facet_var Character. Variable for faceting. Default: "season"
#' @return ggplot object
plot_anomaly_dotline <- function(profile_ts, attr_df, gas = "CH4",
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

  # Mean + SE per height per condition
  height_stats <- df %>%
    group_by(across(all_of(c(facet_var, "tod"))), height_m) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
      n = sum(!is.na(concentration)),
      .groups = "drop"
    )

  # Top-of-tower reference per condition
  top_ref <- height_stats %>%
    filter(height_m == top_height) %>%
    select(all_of(c(facet_var, "tod")), background = mean_conc, background_se = se_conc)

  # Anomaly: SE combines profile and reference uncertainty in quadrature
  anomaly <- height_stats %>%
    left_join(top_ref, by = c(facet_var, "tod")) %>%
    mutate(
      anomaly = mean_conc - background,
      anomaly_se = sqrt(se_conc^2 + background_se^2)
    )

  title_text <- paste0(gas, " Anomaly from Above-Canopy Atmosphere")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  p <- ggplot(anomaly, aes(x = anomaly, y = height_m, color = tod, group = tod)) +
    # Zero reference line
    geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.4) +
    # Canopy reference (drawn before points so points sit on top)
    {if (!is.na(canopy_h))
      geom_hline(yintercept = canopy_h, linetype = "dashed",
                 color = "forestgreen", linewidth = 0.5)
    } +
    # Lines connecting heights
    geom_path(linewidth = 0.7) +
    # Points with horizontal error bars (SE)
    geom_pointrange(aes(xmin = anomaly - anomaly_se,
                        xmax = anomaly + anomaly_se),
                    size = 0.5, linewidth = 0.5) +
    scale_color_manual(
      values = c("Day" = "#E6A817", "Night" = "#2C3E50"),
      name = NULL
    ) +
    scale_y_continuous(
      breaks = sort(unique(height_lookup$height_m)),
      labels = paste(round(sort(unique(height_lookup$height_m)), 1), "m")
    ) +
    labs(
      x = paste0(gas, " anomaly from tower top (", gas_units(gas), ")"),
      y = "Height (m)",
      title = title_text,
      subtitle = "Positive = enriched vs atmosphere | Negative = depleted"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      legend.position = "top"
    )

  # Canopy annotation on last panel only is tricky; skip for simplicity
  # Facet by season
  p <- p + facet_wrap(reformulate(facet_var), nrow = 1)

  p
}
