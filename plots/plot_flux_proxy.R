# plot_flux_proxy.R
# "Flux proxy" — concentration anomaly divided by wind speed at that height
# Normalizes out height-varying mixing, giving a rough source/sink intensity

#' Cross-site flux-proxy profiles (JJA Night)
#'
#' flux_proxy = anomaly / ubar at each height
#' Units: ppb·s/m (CH4) or ppm·s/m (CO2) — relative/proportional to flux
#'
#' @param all_profiles Named list (by site) of profile_ts tibbles
#' @param all_met Named list (by site) of met time series
#' @param all_winds Named list (by site) of wind profile tibbles
#' @param all_data Named list (by site) of loaded site data
#' @param classifications Tibble. Output of classify_sites()
#' @param gas Character. Gas name
#' @param min_wind Numeric. Floor on wind speed (m/s) to avoid division-by-near-zero.
#'   Default: 0.05 (5 cm/s)
#' @return ggplot object
plot_multi_site_flux_proxy <- function(all_profiles, all_met, all_winds,
                                         all_data, classifications,
                                         gas = "CH4", min_wind = 0.05) {

  rows <- list()

  for (site in names(all_profiles)) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]
    wind <- all_winds[[site]]

    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next
    if (is.null(wind) || nrow(wind) == 0) next
    if (!site %in% classifications$site) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

    # JJA Night anomaly per height
    df <- prof_ts %>%
      inner_join(height_lookup, by = "TowerPosition") %>%
      add_season("timeMid")

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

    anomaly_df <- df %>%
      filter(season == "JJA", tod == "Night") %>%
      group_by(height_m) %>%
      summarise(mean_conc = mean(concentration, na.rm = TRUE), .groups = "drop")

    if (nrow(anomaly_df) == 0) next
    bg <- anomaly_df$mean_conc[anomaly_df$height_m == top_height]
    anomaly_df <- anomaly_df %>% mutate(anomaly = mean_conc - bg)

    # Wind speed per height under same condition
    wind_tmp <- wind %>% add_season("timeMid")
    if (!is.null(met) && "PAR" %in% names(met)) {
      wind_tmp <- wind_tmp %>%
        mutate(time_round = round_date(timeMid, "30 minutes")) %>%
        left_join(met %>% select(timeMid, PAR) %>% rename(time_round = timeMid),
                  by = "time_round") %>%
        add_tod_par("PAR") %>%
        select(-time_round, -PAR)
    } else {
      wind_tmp <- wind_tmp %>% add_tod("timeMid")
    }
    wind_df <- wind_tmp %>%
      filter(season == "JJA", tod == "Night") %>%
      group_by(height_m) %>%
      summarise(mean_u = mean(ubar, na.rm = TRUE), .groups = "drop")

    # Join anomaly and wind by height
    combined <- anomaly_df %>%
      inner_join(wind_df, by = "height_m") %>%
      mutate(
        flux_proxy = anomaly / pmax(mean_u, min_wind),
        site = site
      )

    rows[[site]] <- combined
  }

  if (length(rows) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  all_fp <- bind_rows(rows) %>%
    left_join(classifications %>% select(site, combined_class, ground_anomaly),
              by = "site") %>%
    mutate(facet_label = paste0(site, "\n", combined_class))

  label_order <- classifications %>%
    arrange(combined_class, desc(abs(ground_anomaly))) %>%
    mutate(facet_label = paste0(site, "\n", combined_class)) %>%
    pull(facet_label)
  all_fp$facet_label <- factor(all_fp$facet_label,
                                 levels = intersect(label_order, unique(all_fp$facet_label)))

  ggplot(all_fp, aes(x = flux_proxy, y = height_m, color = combined_class)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_path(linewidth = 0.6) +
    geom_point(size = 1.5) +
    scale_color_manual(
      values = c(
        "soil source / mid-canopy feature" = "#8B0000",
        "soil source / monotonic"          = "#E74C3C",
        "soil sink / mid-canopy feature"   = "#1F4E79",
        "soil sink / monotonic"            = "#3498DB",
        "neutral / mid-canopy feature"     = "#7F8C8D",
        "neutral / monotonic"              = "#BDC3C7"
      ),
      name = "Classification"
    ) +
    labs(
      x = paste0(gas, " / wind speed   [",
                  gas_units(gas), " \u00b7 s / m]"),
      y = "Height (m)",
      title = paste0("Cross-site ", gas, " Flux Proxy (JJA Night)"),
      subtitle = "(Anomaly from tower top) / (wind speed at that height). Proportional to flux; accounts for weaker mixing near ground."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 7.5),
      legend.position = "top",
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
      panel.spacing = unit(0.3, "lines")
    ) +
    facet_wrap(~ facet_label, scales = "free_x", ncol = 6)
}
