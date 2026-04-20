# plot_flux_profile.R
# Layer-by-layer flux and source/sink profile plots

#' Vertical profile of mean flux at each pair midpoint, faceted by site
#'
#' @param layer_fluxes_all Named list (by site) of layer flux tibbles
#' @param canopy_metrics Tibble
#' @param gas Character
#' @param condition Character or NULL. "JJA_Night" to filter.
#' @return ggplot object
plot_flux_vertical_profile <- function(layer_fluxes_all, canopy_metrics,
                                         gas = "CH4",
                                         condition = "JJA_Night") {

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next

    # Add season/tod
    lf <- lf %>%
      rename(timeMid = time_round) %>%
      add_season("timeMid") %>%
      add_tod("timeMid")

    if (!is.null(condition) && condition == "JJA_Night") {
      lf <- lf %>% filter(season == "JJA", tod == "Night")
    }
    if (nrow(lf) == 0) next

    # Summarize per sensor (flux labeled at the lower sensor of each pair)
    summary <- lf %>%
      group_by(pair_label, z) %>%
      summarise(
        F_mean = mean(F_FG, na.rm = TRUE),
        F_se = sd(F_FG, na.rm = TRUE) / sqrt(sum(!is.na(F_FG))),
        n = sum(!is.na(F_FG)),
        .groups = "drop"
      ) %>%
      mutate(site = site)

    rows[[site]] <- summary
  }
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  # Order sites by canopy height
  site_order <- canopy_metrics %>% arrange(canopy_height) %>% pull(site)
  df$site <- factor(df$site, levels = intersect(site_order, unique(df$site)))

  units <- flux_units_for_gas(gas)
  scale <- flux_scale_for_gas(gas)
  df    <- df %>% mutate(F_mean = F_mean * scale, F_se = F_se * scale)

  # Arrange by z within each site so the connecting line follows height
  df <- df %>% arrange(site, z)

  ggplot(df, aes(x = F_mean, y = z, color = canopy_class, group = site)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_line(orientation = "y", linewidth = 0.6) +
    geom_pointrange(aes(xmin = F_mean - F_se, xmax = F_mean + F_se),
                    size = 0.3, linewidth = 0.4) +
    # Canopy height marker
    geom_hline(data = df %>% distinct(site, canopy_height),
               aes(yintercept = canopy_height),
               linetype = "dashed", color = "forestgreen", linewidth = 0.4) +
    scale_color_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class", guide = "none") +
    facet_wrap(~ site, scales = "free", ncol = 6) +
    labs(
      x = paste0("F_FG (", units, ")"),
      y = "Pair midpoint height (m)",
      title = paste0(gas, " Vertical Flux Profile (", condition, ")"),
      subtitle = "Flux at each adjacent-pair midpoint. Dashed green = canopy height."
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )
}


#' Layer source/sink strength (difference between adjacent pair fluxes)
plot_layer_source_sink <- function(layer_fluxes_all, canopy_metrics,
                                     gas = "CH4",
                                     condition = "JJA_Night") {

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next

    lf <- lf %>%
      rename(timeMid = time_round) %>%
      add_season("timeMid") %>%
      add_tod("timeMid")

    if (!is.null(condition) && condition == "JJA_Night") {
      lf <- lf %>% filter(season == "JJA", tod == "Night")
    }
    if (nrow(lf) == 0) next

    # Summarize per sensor (flux labeled at the lower sensor of each pair)
    level_summary <- lf %>%
      group_by(pair_label, z) %>%
      summarise(F_mean = mean(F_FG, na.rm = TRUE), .groups = "drop") %>%
      arrange(z)

    if (nrow(level_summary) < 2) next

    # Source/sink in each layer between adjacent sensors:
    #   S(layer between sensor i and sensor i+1) = F(z_{i+1}) − F(z_i)
    # No ground boundary is imposed — F at the lowest sensor IS our estimate
    # of the flux at that height (and, by proxy, the near-surface / soil flux).
    z_lev  <- level_summary$z
    F_lev  <- level_summary$F_mean

    layers <- tibble(
      z_lower_layer = z_lev[-length(z_lev)],
      z_upper_layer = z_lev[-1],
      z_mid_layer   = (z_lev[-length(z_lev)] + z_lev[-1]) / 2,
      source_sink   = diff(F_lev),
      site = site
    )

    rows[[site]] <- layers
  }
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  site_order <- canopy_metrics %>% arrange(canopy_height) %>% pull(site)
  df$site <- factor(df$site, levels = intersect(site_order, unique(df$site)))

  units <- flux_units_for_gas(gas)
  scale <- flux_scale_for_gas(gas)
  df    <- df %>% mutate(source_sink = source_sink * scale)

  df <- df %>% arrange(site, z_mid_layer)

  ggplot(df, aes(x = source_sink, y = z_mid_layer, color = canopy_class, group = site)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = source_sink,
                     y = z_mid_layer, yend = z_mid_layer),
                 linewidth = 1.5, alpha = 0.7) +
    geom_point(size = 2) +
    geom_hline(data = df %>% distinct(site, canopy_height),
               aes(yintercept = canopy_height),
               linetype = "dashed", color = "forestgreen", linewidth = 0.4) +
    scale_color_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class", guide = "none") +
    facet_wrap(~ site, scales = "free", ncol = 6) +
    labs(
      x = paste0("Layer source/sink (", units, ")"),
      y = "Layer midpoint height (m)",
      title = paste0(gas, " Layer Source/Sink Profile (", condition, ")"),
      subtitle = "Positive = net source in that layer. Negative = net sink."
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )
}
