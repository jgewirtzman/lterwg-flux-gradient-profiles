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


#' Layer source/sink strength per layer, with optional storage correction.
#'
#' @param layer_fluxes_all Named list (by site) of layer flux tibbles from
#'   compute_layer_fluxes(); used as fallback when source_sink_all is NULL.
#' @param source_sink_all Optional named list (by site) of tibbles from
#'   compute_source_sink_layer(). When provided, uses S_net (= dF + storage)
#'   as the primary source/sink estimate. A secondary marker shows dF alone
#'   so the storage contribution is visible.
plot_layer_source_sink <- function(layer_fluxes_all, canopy_metrics,
                                     gas = "CH4",
                                     condition = "JJA_Night",
                                     source_sink_all = NULL) {

  use_storage <- !is.null(source_sink_all)
  rows <- list()
  site_pool <- if (use_storage) names(source_sink_all) else names(layer_fluxes_all)

  for (site in site_pool) {
    if (use_storage) {
      ss <- source_sink_all[[site]]
      if (is.null(ss) || nrow(ss) == 0) next
      ss <- ss %>%
        rename(timeMid = time_round) %>%
        add_season("timeMid") %>%
        add_tod("timeMid")
      if (!is.null(condition) && condition == "JJA_Night") {
        ss <- ss %>% filter(season == "JJA", tod == "Night")
      }
      if (nrow(ss) == 0) next

      layers <- ss %>%
        group_by(layer_label, z_lower, z_upper) %>%
        summarise(
          source_sink     = mean(S_net, na.rm = TRUE),
          source_sink_dF  = mean(S_net_nostorage, na.rm = TRUE),
          has_storage     = any(!is.na(storage_flux) & storage_flux != 0),
          .groups = "drop"
        ) %>%
        mutate(z_mid_layer = (z_lower + z_upper) / 2, site = site)
    } else {
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

      level_summary <- lf %>%
        group_by(pair_label, z) %>%
        summarise(F_mean = mean(F_FG, na.rm = TRUE), .groups = "drop") %>%
        arrange(z)
      if (nrow(level_summary) < 2) next

      z_lev <- level_summary$z
      F_lev <- level_summary$F_mean
      layers <- tibble(
        z_lower      = z_lev[-length(z_lev)],
        z_upper      = z_lev[-1],
        z_mid_layer  = (z_lev[-length(z_lev)] + z_lev[-1]) / 2,
        source_sink  = diff(F_lev),
        source_sink_dF = diff(F_lev),
        has_storage  = FALSE,
        layer_label  = paste0(round(z_lev[-length(z_lev)], 1), "\u2013",
                               round(z_lev[-1], 1), " m"),
        site = site
      )
    }
    rows[[site]] <- layers
  }

  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  site_order <- canopy_metrics %>% arrange(canopy_height) %>% pull(site)
  df$site <- factor(df$site, levels = intersect(site_order, unique(df$site)))

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  df <- df %>%
    mutate(source_sink    = source_sink    * scale,
           source_sink_dF = source_sink_dF * scale) %>%
    arrange(site, z_mid_layer)

  subtitle_txt <- if (use_storage)
    "Thick bar = S_net (dF + storage). Hollow circle = dF alone (no storage)."
  else
    "Positive = net source in that layer. Negative = net sink. (storage not included)"

  p <- ggplot(df, aes(y = z_mid_layer, color = canopy_class, group = site)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    # Primary bar: S_net (storage-corrected when source_sink_all provided)
    geom_segment(aes(x = 0, xend = source_sink,
                     y = z_mid_layer, yend = z_mid_layer),
                 linewidth = 1.5, alpha = 0.75) +
    geom_point(aes(x = source_sink), size = 2)

  # When using storage: overlay hollow circle at dF-only position
  if (use_storage) {
    p <- p +
      geom_point(aes(x = source_sink_dF), size = 2.5, shape = 1,
                 color = "grey40", alpha = 0.7)
  }

  p +
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
      subtitle = subtitle_txt
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text    = element_text(face = "bold", size = 8),
      panel.border  = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )
}
