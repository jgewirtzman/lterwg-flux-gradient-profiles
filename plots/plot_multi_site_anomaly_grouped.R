# plot_multi_site_anomaly_grouped.R
# Cross-site CH4 anomaly profiles grouped by classification
# (soil sink / soil source / neutral x monotonic / mid-canopy feature)

#' Multi-site anomaly grid grouped by soil + structure classification
#'
#' @param all_profiles Named list (by site) of profile_ts tibbles
#' @param all_met Named list (by site) of met time series
#' @param all_data Named list (by site) of loaded site data
#' @param classifications Tibble. Output of classify_sites()
#' @param gas Character. Gas name
#' @param jja_night_only Logical. Show only JJA Night (classification basis).
#'   Default TRUE.
#' @return ggplot object
plot_multi_site_anomaly_grouped <- function(all_profiles, all_met, all_data,
                                              classifications,
                                              gas = "CH4",
                                              jja_night_only = TRUE) {

  rows <- list()

  for (site in names(all_profiles)) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]

    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next
    if (!site %in% classifications$site) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

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

    # Compute anomaly per season x tod
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
        site = site
      )

    if (jja_night_only) {
      anomaly <- anomaly %>% filter(season == "JJA", tod == "Night")
    }

    rows[[site]] <- anomaly
  }

  if (length(rows) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  all_anomaly <- bind_rows(rows) %>%
    left_join(classifications %>% select(site, combined_class, ground_anomaly),
              by = "site")

  # Within each class, order by ground anomaly magnitude (most extreme first)
  all_anomaly <- all_anomaly %>%
    group_by(combined_class) %>%
    mutate(site_ord = reorder(site, -abs(ground_anomaly))) %>%
    ungroup() %>%
    mutate(site = factor(site, levels = unique(site_ord[order(combined_class, -abs(ground_anomaly))])))

  # Simpler: order sites by class first then by ground anomaly
  site_order <- classifications %>%
    arrange(combined_class, desc(abs(ground_anomaly))) %>%
    pull(site)
  all_anomaly$site <- factor(all_anomaly$site, levels = intersect(site_order, unique(all_anomaly$site)))

  subtitle <- if (jja_night_only) {
    "JJA Night only. Grouped by classification based on JJA Night profile."
  } else {
    "All seasons and TOD shown. Grouped by JJA-Night classification."
  }

  p <- ggplot(all_anomaly, aes(x = anomaly, y = height_m)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_path(linewidth = 0.6, color = "#2C3E50") +
    geom_pointrange(aes(xmin = anomaly - anomaly_se,
                        xmax = anomaly + anomaly_se),
                    color = "#2C3E50", size = 0.35, linewidth = 0.4) +
    labs(
      x = paste0(gas, " anomaly from tower top (", gas_units(gas), ")"),
      y = "Height (m)",
      title = paste0("Cross-site ", gas,
                      " Anomaly Grouped by Profile Classification"),
      subtitle = subtitle
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      strip.text.x = element_text(size = 8),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3),
      panel.spacing = unit(0.3, "lines")
    ) +
    facet_wrap(~ combined_class + site, scales = "free", ncol = 6,
               labeller = labeller(.multi_line = FALSE))

  p
}


#' Same as above but with a simpler facet structure: rows = class, cols = site
plot_multi_site_anomaly_grouped_grid <- function(all_profiles, all_met, all_data,
                                                   classifications,
                                                   gas = "CH4") {

  rows <- list()
  for (site in names(all_profiles)) {
    prof_ts <- all_profiles[[site]][[gas]]
    attr <- all_data[[site]]$attr
    met <- all_met[[site]]
    if (is.null(prof_ts) || nrow(prof_ts) == 0 || is.null(attr)) next
    if (!site %in% classifications$site) next

    height_lookup <- get_height_lookup(attr)
    top_height <- max(height_lookup$height_m)

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

    height_stats <- df %>%
      filter(season == "JJA", tod == "Night") %>%
      group_by(height_m) %>%
      summarise(
        mean_conc = mean(concentration, na.rm = TRUE),
        se_conc = sd(concentration, na.rm = TRUE) / sqrt(sum(!is.na(concentration))),
        .groups = "drop"
      )

    if (nrow(height_stats) == 0) next

    background <- height_stats$mean_conc[height_stats$height_m == top_height]
    background_se <- height_stats$se_conc[height_stats$height_m == top_height]

    rows[[site]] <- height_stats %>%
      mutate(
        anomaly = mean_conc - background,
        anomaly_se = sqrt(se_conc^2 + background_se^2),
        site = site
      )
  }

  if (length(rows) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  all_anomaly <- bind_rows(rows) %>%
    left_join(classifications %>% select(site, combined_class, ground_anomaly),
              by = "site")

  # Add site labels with classification info
  all_anomaly <- all_anomaly %>%
    mutate(
      facet_label = paste0(site, "\n", combined_class)
    )

  # Order by class then by ground anomaly magnitude
  label_order <- classifications %>%
    arrange(combined_class, desc(abs(ground_anomaly))) %>%
    mutate(facet_label = paste0(site, "\n", combined_class)) %>%
    pull(facet_label)

  all_anomaly$facet_label <- factor(all_anomaly$facet_label,
                                      levels = intersect(label_order, unique(all_anomaly$facet_label)))

  ggplot(all_anomaly, aes(x = anomaly, y = height_m, color = combined_class)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_path(linewidth = 0.6) +
    geom_pointrange(aes(xmin = anomaly - anomaly_se,
                        xmax = anomaly + anomaly_se),
                    size = 0.35, linewidth = 0.4) +
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
      x = paste0(gas, " anomaly from tower top (", gas_units(gas), ")"),
      y = "Height (m)",
      title = paste0("Cross-site ", gas, " JJA Night Profiles"),
      subtitle = "Grouped by soil class (sink/source/neutral) x structure (monotonic/mid-canopy feature)"
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
