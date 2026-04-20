# plot_anomaly_with_wind.R
# Side-by-side: concentration anomaly profile + wind speed profile
# Lets you see whether big anomalies correspond to stagnant air

library(patchwork)

#' Multi-site grid: anomaly profile + wind speed profile side-by-side per site
#'
#' For each site, builds a two-panel plot (anomaly + wind) and assembles them
#' into a grid. JJA Night only (strongest signal).
#'
#' @param all_profiles Named list (by site) of profile_ts tibbles
#' @param all_met Named list (by site) of met time series
#' @param all_winds Named list (by site) of wind profile tibbles
#' @param all_data Named list (by site) of loaded site data
#' @param classifications Tibble. Output of classify_sites() for ordering
#' @param gas Character. Gas name
#' @return patchwork object
plot_multi_site_anomaly_with_wind <- function(all_profiles, all_met, all_winds,
                                                all_data, classifications,
                                                gas = "CH4") {

  plots <- list()

  # Order sites by classification
  site_order <- classifications %>%
    arrange(combined_class, desc(abs(ground_anomaly))) %>%
    pull(site)

  site_order <- intersect(site_order, names(all_profiles))

  for (site in site_order) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]
    wind <- all_winds[[site]]

    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next
    if (!site %in% classifications$site) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

    # Compute JJA Night anomaly
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
      summarise(
        mean_conc = mean(concentration, na.rm = TRUE),
        se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
        .groups = "drop"
      )

    if (nrow(anomaly_df) == 0) next

    bg <- anomaly_df$mean_conc[anomaly_df$height_m == top_height]
    bg_se <- anomaly_df$se_conc[anomaly_df$height_m == top_height]
    anomaly_df <- anomaly_df %>%
      mutate(
        anomaly = mean_conc - bg,
        anomaly_se = sqrt(se_conc^2 + bg_se^2)
      )

    # Wind profile for JJA Night
    wind_df <- if (!is.null(wind) && nrow(wind) > 0) {
      tmp <- wind %>% add_season("timeMid")
      if (!is.null(met) && "PAR" %in% names(met)) {
        tmp <- tmp %>%
          mutate(time_round = round_date(timeMid, "30 minutes")) %>%
          left_join(met %>% select(timeMid, PAR) %>% rename(time_round = timeMid),
                    by = "time_round") %>%
          add_tod_par("PAR") %>%
          select(-time_round, -PAR)
      } else {
        tmp <- tmp %>% add_tod("timeMid")
      }
      tmp %>%
        filter(season == "JJA", tod == "Night") %>%
        group_by(height_m) %>%
        summarise(
          mean_u = mean(ubar, na.rm = TRUE),
          se_u = sd(ubar, na.rm = TRUE) / sqrt(sum(!is.na(ubar))),
          .groups = "drop"
        )
    } else tibble()

    cls <- classifications$combined_class[classifications$site == site]

    # Classification color
    cls_color <- c(
      "soil source / mid-canopy feature" = "#8B0000",
      "soil source / monotonic"          = "#E74C3C",
      "soil sink / mid-canopy feature"   = "#1F4E79",
      "soil sink / monotonic"            = "#3498DB",
      "neutral / mid-canopy feature"     = "#7F8C8D",
      "neutral / monotonic"              = "#BDC3C7"
    )[as.character(cls)]

    p_anom <- ggplot(anomaly_df, aes(x = anomaly, y = height_m)) +
      geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
      geom_path(linewidth = 0.5, color = cls_color) +
      geom_pointrange(aes(xmin = anomaly - anomaly_se, xmax = anomaly + anomaly_se),
                      size = 0.25, linewidth = 0.3, color = cls_color) +
      labs(x = paste0("\u0394", gas), y = NULL,
           title = site, subtitle = as.character(cls)) +
      theme_minimal(base_size = 8) +
      theme(
        plot.title = element_text(face = "bold", size = 9, color = cls_color),
        plot.subtitle = element_text(size = 6, color = "grey40"),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank()
      )

    p_wind <- if (nrow(wind_df) > 0) {
      ggplot(wind_df, aes(x = mean_u, y = height_m)) +
        geom_path(linewidth = 0.5, color = "#2C3E50") +
        geom_pointrange(aes(xmin = mean_u - se_u, xmax = mean_u + se_u),
                        size = 0.25, linewidth = 0.3, color = "#2C3E50") +
        labs(x = "wind (m/s)", y = NULL, title = NULL, subtitle = NULL) +
        theme_minimal(base_size = 8) +
        theme(
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 6),
          axis.title = element_text(size = 7),
          panel.grid.minor = element_blank()
        )
    } else {
      ggplot() + theme_void()
    }

    plots[[site]] <- p_anom + p_wind + plot_layout(widths = c(1, 0.6))
  }

  if (length(plots) == 0) return(ggplot() + theme_void())

  wrap_plots(plots, ncol = 4) +
    plot_annotation(
      title = paste0("Cross-site ", gas,
                      " JJA Night Anomaly + Wind Speed Profile"),
      subtitle = "Sites grouped by classification (color). Left=\u0394conc from tower top; right=wind speed",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 10))
    )
}
