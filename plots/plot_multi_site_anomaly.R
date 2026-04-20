# plot_multi_site_anomaly.R
# Cross-site comparison of anomaly-from-atmosphere profiles
#
# Grid: rows = sites, cols = seasons; day/night overlaid

#' Multi-site anomaly-from-atmosphere profile grid
#'
#' @param all_profiles Named list (by site) of profile_ts tibbles
#' @param all_met Named list (by site) of met time series tibbles
#' @param all_data Named list (by site) of loaded site data (needs $attr)
#' @param gas Character. Gas name
#' @param jja_only Logical. If TRUE, show only JJA (summer) panels in a single row
#' @return ggplot object
plot_multi_site_anomaly <- function(all_profiles, all_met, all_data,
                                      gas = "CH4", jja_only = FALSE) {

  rows <- list()

  for (site in names(all_profiles)) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]

    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

    df <- prof_ts %>%
      inner_join(height_lookup, by = "TowerPosition") %>%
      add_season("timeMid")

    # PAR-based or fixed day/night
    if (!is.null(met) && "PAR" %in% names(met)) {
      df <- df %>%
        mutate(time_round = round_date(timeMid, "30 minutes")) %>%
        left_join(met %>% select(timeMid, PAR) %>% rename(time_round = timeMid),
                  by = "time_round") %>%
        add_tod_par("PAR") %>%
        select(-time_round, -PAR)
    } else {
      df <- df %>% add_tod("timeMid")
    }

    # Mean + SE per height per condition
    height_stats <- df %>%
      group_by(season, tod, height_m) %>%
      summarise(
        mean_conc = mean(concentration, na.rm = TRUE),
        se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
        .groups = "drop"
      )

    top_ref <- height_stats %>%
      filter(height_m == top_height) %>%
      select(season, tod, background = mean_conc, background_se = se_conc)

    anomaly <- height_stats %>%
      left_join(top_ref, by = c("season", "tod")) %>%
      mutate(
        anomaly = mean_conc - background,
        anomaly_se = sqrt(se_conc^2 + background_se^2),
        site = site,
        canopy_height = get_canopy_height(attr),
        height_rel = height_m / canopy_height  # height relative to canopy
      )

    rows[[site]] <- anomaly
  }

  if (length(rows) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  all_anomaly <- bind_rows(rows)

  # Filter to JJA if requested
  if (jja_only) {
    all_anomaly <- all_anomaly %>% filter(season == "JJA")
  }

  # Order sites by max canopy height (tallest at top)
  site_order <- all_anomaly %>%
    group_by(site) %>%
    summarise(max_h = max(height_m, na.rm = TRUE)) %>%
    arrange(desc(max_h)) %>%
    pull(site)

  all_anomaly <- all_anomaly %>%
    mutate(site = factor(site, levels = site_order))

  p <- ggplot(all_anomaly, aes(x = anomaly, y = height_m, color = tod, group = tod)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_path(linewidth = 0.5) +
    geom_pointrange(aes(xmin = anomaly - anomaly_se,
                        xmax = anomaly + anomaly_se),
                    size = 0.3, linewidth = 0.3) +
    scale_color_manual(values = c("Day" = "#E6A817", "Night" = "#2C3E50"),
                        name = NULL) +
    labs(
      x = paste0(gas, " anomaly from tower top (", gas_units(gas), ")"),
      y = "Height (m)",
      title = paste0("Cross-site ", gas, " Anomaly from Above-Canopy Atmosphere"),
      subtitle = "Rows ordered by tower height (tallest on top)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "top",
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )

  # Faceting
  if (jja_only) {
    p <- p + facet_wrap(~ site, scales = "free", ncol = 4)
  } else {
    p <- p + facet_grid(site ~ season, scales = "free")
  }

  p
}
