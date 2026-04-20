# plot_cross_site_synthesis.R
# Higher-level cross-site summary plots designed to be readable at a glance,
# rather than faceting over all 47 sites.
#
# - plot_flux_profile_by_class():   1-panel-per-canopy-class median profiles
#                                   (5 panels) with IQR ribbon across sites
# - plot_flux_ranked_bars():        one bar per site, top-pair flux magnitude,
#                                   colored by canopy class
# - plot_flux_profile_exemplar():   detailed profile plot for a curated subset
#                                   of ecosystem archetype sites


CANOPY_COLORS <- c(
  "No canopy (<1 m)"        = "#E8C547",
  "Short canopy (1-5 m)"    = "#A8D95B",
  "Medium canopy (5-15 m)"  = "#5BBA6F",
  "Tall canopy (15-30 m)"   = "#2D8B3A",
  "Very tall canopy (>30 m)"= "#0F5A1E"
)


#' Internal helper: aggregate per-site per-z flux stats from layer_fluxes_all,
#' with optional JJA Night filter. Returns one row per (site, z).
#'
#' @param layer_fluxes_all Named list (by site)
#' @param filter_jja_night Logical. If TRUE, keep only JJA Night timestamps.
#' @return Tibble with columns site, z, F_med, F_q25, F_q75, n
aggregate_flux_profiles <- function(layer_fluxes_all, filter_jja_night = TRUE) {
  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    if (filter_jja_night) {
      lf <- lf %>%
        mutate(month = month(time_round), hour = hour(time_round),
               season = case_when(month %in% 6:8 ~ "JJA",
                                  month %in% 3:5 ~ "MAM",
                                  month %in% 9:11 ~ "SON",
                                  TRUE ~ "DJF"),
               tod = ifelse(hour >= 22 | hour < 6, "Night", "Day")) %>%
        filter(season == "JJA", tod == "Night")
      if (nrow(lf) == 0) next
    }
    per_z <- lf %>%
      group_by(z) %>%
      summarise(F_med = median(F_FG, na.rm = TRUE),
                F_q25 = quantile(F_FG, 0.25, na.rm = TRUE),
                F_q75 = quantile(F_FG, 0.75, na.rm = TRUE),
                n = sum(!is.na(F_FG)),
                .groups = "drop") %>%
      arrange(z) %>% mutate(site = site)
    rows[[site]] <- per_z
  }
  bind_rows(rows)
}


#' Median flux profile per canopy class, with IQR ribbon across sites.
#'
#' For each sensor (within a site, JJA Night median), normalizes z by the
#' site's canopy height so profiles from different-height towers can be
#' combined. Then, within each canopy class, plots the median F across sites
#' vs relative height, with an inter-site IQR ribbon.
#'
#' @param layer_fluxes_all Named list (by site) of layer flux tibbles
#' @param canopy_metrics Tibble
#' @param gas Character
#' @return ggplot object
plot_flux_profile_by_class <- function(layer_fluxes_all, canopy_metrics,
                                         gas = "CH4") {

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next

    ch <- canopy_metrics$canopy_height[canopy_metrics$site == site]
    cc <- canopy_metrics$canopy_class[canopy_metrics$site == site]
    if (length(ch) == 0 || is.na(ch) || ch < 0.1) ch <- 1  # fallback for open sites

    lf <- lf %>%
      mutate(month = month(time_round), hour = hour(time_round),
             season = case_when(month %in% 6:8 ~ "JJA",
                                month %in% 3:5 ~ "MAM",
                                month %in% 9:11 ~ "SON",
                                TRUE ~ "DJF"),
             tod = ifelse(hour >= 22 | hour < 6, "Night", "Day")) %>%
      filter(season == "JJA", tod == "Night")
    if (nrow(lf) == 0) next

    per_z <- lf %>%
      group_by(z) %>%
      summarise(F_med = median(F_FG, na.rm = TRUE), .groups = "drop") %>%
      mutate(site = site, canopy_class = cc, canopy_h = ch,
             z_rel = z / ch)
    rows[[site]] <- per_z
  }

  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  # Bin by relative height within each canopy class and compute ensemble stats.
  # Using smooth continuous z_rel directly with IQR over sites at each site's
  # sampled z_rel — approximate via quantile over a 5% sliding window.
  df <- df %>% arrange(canopy_class, z_rel)

  # Simple approach: aggregate to fixed z_rel bins per class
  df <- df %>%
    mutate(z_rel_bin = cut(z_rel, breaks = seq(0, 3, by = 0.1),
                             include.lowest = TRUE))
  class_summary <- df %>%
    group_by(canopy_class, z_rel_bin) %>%
    summarise(
      z_rel_mid = mean(as.numeric(z_rel), na.rm = TRUE),
      F_med = median(F_med, na.rm = TRUE),
      F_q25 = quantile(F_med, 0.25, na.rm = TRUE),
      F_q75 = quantile(F_med, 0.75, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 2, !is.na(z_rel_mid)) %>%
    arrange(canopy_class, z_rel_mid)

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  class_summary <- class_summary %>%
    mutate(F_med = F_med * scale,
           F_q25 = F_q25 * scale, F_q75 = F_q75 * scale)

  ggplot(class_summary,
         aes(x = F_med, y = z_rel_mid, color = canopy_class, fill = canopy_class)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "forestgreen",
               linewidth = 0.4) +
    geom_ribbon(aes(xmin = F_q25, xmax = F_q75, group = canopy_class),
                alpha = 0.25, color = NA) +
    geom_line(orientation = "y", linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_color_manual(values = CANOPY_COLORS, guide = "none") +
    scale_fill_manual(values = CANOPY_COLORS, guide = "none") +
    facet_wrap(~ canopy_class, nrow = 1, scales = "free_x") +
    labs(
      x = paste0(gas, " flux F (", units, ")"),
      y = "Relative height (z / canopy height)",
      title = paste0(gas, " Flux Profile by Canopy Class (JJA Night)"),
      subtitle = "Median across sites in each class; ribbon = inter-site IQR. Green dashed = canopy top."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
    )
}


#' One-bar-per-site ranking of top-pair flux (JJA Night median).
#'
#' Colored by canopy class, ordered by flux magnitude within class.
#'
#' @param layer_fluxes_all Named list (by site)
#' @param canopy_metrics Tibble
#' @param gas Character
#' @return ggplot object
plot_flux_ranked_bars <- function(layer_fluxes_all, canopy_metrics,
                                    gas = "CH4") {
  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    lf <- lf %>%
      mutate(month = month(time_round), hour = hour(time_round),
             season = case_when(month %in% 6:8 ~ "JJA",
                                month %in% 3:5 ~ "MAM",
                                month %in% 9:11 ~ "SON",
                                TRUE ~ "DJF"),
             tod = ifelse(hour >= 22 | hour < 6, "Night", "Day")) %>%
      filter(season == "JJA", tod == "Night")
    if (nrow(lf) == 0) next

    z_vals <- sort(unique(lf$z))
    if (length(z_vals) < 2) next
    z_top <- z_vals[length(z_vals)]
    ftop <- lf$F_FG[lf$z == z_top]
    ftop <- ftop[!is.na(ftop)]
    if (length(ftop) < 20) next

    rows[[site]] <- tibble(
      site = site,
      F_top_med = median(ftop),
      F_top_q25 = quantile(ftop, 0.25),
      F_top_q75 = quantile(ftop, 0.75),
      n = length(ftop)
    )
  }
  stats <- bind_rows(rows) %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site") %>%
    filter(!is.na(canopy_class))

  if (nrow(stats) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  # Order sites by canopy class then by flux magnitude
  cc_order <- c("No canopy (<1 m)", "Short canopy (1-5 m)",
                 "Medium canopy (5-15 m)", "Tall canopy (15-30 m)",
                 "Very tall canopy (>30 m)")
  stats <- stats %>%
    mutate(canopy_class = factor(canopy_class, levels = cc_order)) %>%
    arrange(canopy_class, F_top_med) %>%
    mutate(site = factor(site, levels = site))

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  stats <- stats %>% mutate(F_top_med = F_top_med * scale,
                             F_top_q25 = F_top_q25 * scale,
                             F_top_q75 = F_top_q75 * scale)

  ggplot(stats, aes(y = site, x = F_top_med, fill = canopy_class)) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_col() +
    geom_errorbarh(aes(xmin = F_top_q25, xmax = F_top_q75),
                   height = 0.3, linewidth = 0.3, color = "grey30") +
    scale_fill_manual(values = CANOPY_COLORS, name = "Canopy class") +
    facet_grid(canopy_class ~ ., scales = "free_y", space = "free_y",
               switch = "y") +
    labs(
      x = paste0(gas, " top-pair flux median (", units, ")"),
      y = NULL,
      title = paste0(gas, " Ecosystem Flux by Site, Ranked (JJA Night)"),
      subtitle = "Error bars = IQR across 30-min timestamps. Sites grouped by canopy class."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
      strip.placement = "outside",
      panel.spacing = unit(0.3, "lines"),
      axis.text.y = element_text(size = 8)
    )
}


#' Detailed vertical flux profiles for a curated subset of exemplar sites,
#' covering each ecosystem archetype with a large enough panel to actually read.
#'
#' Exemplar set (default): WREF, SERC (very tall); HARV, BART (tall);
#' ABBY (medium); LAJA (short tropical); CPER (grassland); SRER (desert);
#' BARR (wetland tundra).
#'
#' @param layer_fluxes_all Named list (by site)
#' @param canopy_metrics Tibble
#' @param gas Character
#' @param sites_exemplar Character vector of site codes. Defaults to a
#'   canonical set.
#' @return ggplot object
plot_flux_profile_exemplar <- function(layer_fluxes_all, canopy_metrics,
                                         gas = "CH4",
                                         sites_exemplar = c("WREF", "SERC",
                                                             "HARV", "BART",
                                                             "ABBY", "LAJA",
                                                             "CPER", "SRER",
                                                             "BARR")) {
  rows <- list()
  for (site in intersect(sites_exemplar, names(layer_fluxes_all))) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    lf <- lf %>%
      mutate(month = month(time_round), hour = hour(time_round),
             season = case_when(month %in% 6:8 ~ "JJA",
                                month %in% 3:5 ~ "MAM",
                                month %in% 9:11 ~ "SON",
                                TRUE ~ "DJF"),
             tod = ifelse(hour >= 22 | hour < 6, "Night", "Day")) %>%
      filter(season == "JJA", tod == "Night")
    if (nrow(lf) == 0) next

    per_z <- lf %>%
      group_by(z) %>%
      summarise(
        F_med = median(F_FG, na.rm = TRUE),
        F_q25 = quantile(F_FG, 0.25, na.rm = TRUE),
        F_q75 = quantile(F_FG, 0.75, na.rm = TRUE),
        n = sum(!is.na(F_FG)),
        .groups = "drop"
      ) %>% arrange(z)
    per_z$site <- site
    rows[[site]] <- per_z
  }

  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  # Order sites by canopy height ascending
  site_order <- canopy_metrics %>%
    filter(site %in% df$site) %>%
    arrange(canopy_height) %>% pull(site)
  df <- df %>% mutate(site = factor(site, levels = site_order)) %>%
    arrange(site, z)

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  df <- df %>% mutate(F_med = F_med * scale,
                       F_q25 = F_q25 * scale,
                       F_q75 = F_q75 * scale)

  ggplot(df, aes(x = F_med, y = z, color = canopy_class, group = site)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = F_q25, xmax = F_q75),
                   height = 0.4, linewidth = 0.4, alpha = 0.6) +
    geom_line(orientation = "y", linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(data = df %>% distinct(site, canopy_height),
               aes(yintercept = canopy_height),
               linetype = "dashed", color = "forestgreen", linewidth = 0.5) +
    scale_color_manual(values = CANOPY_COLORS, name = "Canopy class") +
    facet_wrap(~ site, scales = "free", ncol = 3) +
    labs(
      x = paste0(gas, " flux F (", units, ")"),
      y = "Height (m)",
      title = paste0(gas, " Flux Profiles — Ecosystem Archetypes (JJA Night)"),
      subtitle = "Error bars = IQR across 30-min timestamps. Dashed green = canopy top."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold")
    )
}


#' Per-ecosystem-type flux profile grid — all sites faceted, with a shaded
#' canopy band (from 0 to canopy_height) instead of a dashed canopy line.
#'
#' @param layer_fluxes_all Named list (by site)
#' @param canopy_metrics Tibble
#' @param ecosystem_types Named character vector: site -> ecosystem type
#' @param eco_type Character. Which ecosystem to plot.
#' @param gas Character
#' @param filter_jja_night Logical. TRUE = JJA Night only, FALSE = all data.
#' @return ggplot object
plot_flux_profile_by_ecotype <- function(layer_fluxes_all, canopy_metrics,
                                           ecosystem_types,
                                           eco_type = "Forest - Eastern Deciduous",
                                           gas = "CH4",
                                           filter_jja_night = TRUE) {

  eco_sites <- names(ecosystem_types)[ecosystem_types == eco_type]
  keep <- intersect(eco_sites, names(layer_fluxes_all))
  if (length(keep) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for", eco_type)))
  }

  df <- aggregate_flux_profiles(layer_fluxes_all[keep],
                                  filter_jja_night = filter_jja_night)
  if (nrow(df) == 0) {
    return(ggplot() + theme_void() + labs(title = paste("No data for", eco_type)))
  }

  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  site_order <- canopy_metrics %>%
    filter(site %in% df$site) %>%
    arrange(canopy_height) %>% pull(site)
  df <- df %>%
    mutate(site = factor(site, levels = site_order)) %>%
    arrange(site, z)

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  df <- df %>% mutate(F_med = F_med * scale,
                       F_q25 = F_q25 * scale,
                       F_q75 = F_q75 * scale)

  # Canopy band: from 0 to canopy_height per site
  canopy_band <- df %>%
    distinct(site, canopy_height) %>%
    mutate(canopy_height = pmax(canopy_height, 0.1))  # avoid zero band

  # Pick facet columns based on count
  n_sites <- length(unique(df$site))
  ncol <- if (n_sites <= 3) n_sites
         else if (n_sites <= 6) 3
         else if (n_sites <= 12) 4
         else 5

  cond_label <- if (filter_jja_night) "JJA Night" else "All data"

  ggplot(df, aes(x = F_med, y = z, color = canopy_class, group = site)) +
    # Shaded canopy band (0 to canopy height)
    geom_rect(data = canopy_band,
              aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = canopy_height),
              inherit.aes = FALSE,
              fill = "forestgreen", alpha = 0.12) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = F_q25, xmax = F_q75),
                   height = 0.4, linewidth = 0.4, alpha = 0.6) +
    geom_line(orientation = "y", linewidth = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = CANOPY_COLORS, name = "Canopy class") +
    facet_wrap(~ site, scales = "free", ncol = ncol) +
    labs(
      x = paste0(gas, " flux (", units, ")"),
      y = "Height (m)",
      title = paste0(gas, " Flux Profiles — ", eco_type, " (", cond_label, ")"),
      subtitle = "Green band = canopy (0 to canopy height). Error bars = IQR across 30-min obs."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "top",
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.3, "lines")
    )
}


#' Ground-level flux vs ecosystem-level flux scatter per gas.
#' One point per site; x = bottom-sensor flux (≈ soil), y = top-pair flux
#' (≈ ecosystem), error bars = IQR across 30-min timestamps. Colored by
#' canopy class, labeled with site code.
#'
#' @param layer_fluxes_all Named list (by site)
#' @param canopy_metrics Tibble
#' @param gas Character
#' @param filter_jja_night Logical. TRUE = JJA Night only.
#' @return ggplot object
plot_flux_ground_vs_ecosystem <- function(layer_fluxes_all, canopy_metrics,
                                            gas = "CH4",
                                            filter_jja_night = TRUE,
                                            stat = c("median", "mean")) {
  stat <- match.arg(stat)

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    if (filter_jja_night) {
      lf <- lf %>%
        mutate(month = month(time_round), hour = hour(time_round),
               season = case_when(month %in% 6:8 ~ "JJA",
                                  month %in% 3:5 ~ "MAM",
                                  month %in% 9:11 ~ "SON",
                                  TRUE ~ "DJF"),
               tod = ifelse(hour >= 22 | hour < 6, "Night", "Day")) %>%
        filter(season == "JJA", tod == "Night")
      if (nrow(lf) == 0) next
    }

    z_vals <- sort(unique(lf$z))
    if (length(z_vals) < 2) next
    z_top <- z_vals[length(z_vals)]
    z_bot <- z_vals[1]

    # Take stats INDEPENDENTLY (don't require shared timestamps).
    # NEON 9-min raw data measures one level per period (rotating), so very
    # few 30-min bins have BOTH the top AND bottom sensor — pairing would
    # drop most sites. For a site-summary plot this isn't needed: we just
    # want each site's typical bottom and top flux.
    bot_vals <- lf$F_FG[lf$z == z_bot]; bot_vals <- bot_vals[!is.na(bot_vals)]
    top_vals <- lf$F_FG[lf$z == z_top]; top_vals <- top_vals[!is.na(top_vals)]
    if (length(bot_vals) < 20 || length(top_vals) < 20) next

    rows[[site]] <- tibble(
      site = site,
      z_bot = z_bot,
      F_bot_med  = median(bot_vals),
      F_bot_mean = mean(bot_vals),
      F_bot_q25  = quantile(bot_vals, 0.25),
      F_bot_q75  = quantile(bot_vals, 0.75),
      F_top_med  = median(top_vals),
      F_top_mean = mean(top_vals),
      F_top_q25  = quantile(top_vals, 0.25),
      F_top_q75  = quantile(top_vals, 0.75),
      n_bot = length(bot_vals),
      n_top = length(top_vals)
    )
  }
  stats <- bind_rows(rows)
  if (nrow(stats) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  stats <- stats %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site") %>%
    # Fill if bottom sensor is at/below canopy top (true canopy profile).
    # Open if bottom sensor is already above canopy (just an air column —
    # top-vs-bottom differences more likely driven by footprint than local
    # sources/sinks).
    mutate(within_canopy = z_bot <= canopy_height + 0.01)

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  stats <- stats %>% mutate(
    F_bot_med  = F_bot_med  * scale,
    F_bot_mean = F_bot_mean * scale,
    F_top_med  = F_top_med  * scale,
    F_top_mean = F_top_mean * scale
  )

  # Pick the requested stat
  one <- if (stat == "median") {
    stats %>% transmute(site, canopy_class, within_canopy,
                         F_bot = F_bot_med, F_top = F_top_med)
  } else {
    stats %>% transmute(site, canopy_class, within_canopy,
                         F_bot = F_bot_mean, F_top = F_top_mean)
  }

  cond_label <- if (filter_jja_night) "JJA Night" else "All data"

  # Symmetric square axis range so 1:1 line + 1:1 aspect look right.
  finite_vals <- c(one$F_bot, one$F_top)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  rng <- range(finite_vals, na.rm = TRUE)
  lim <- max(abs(rng)) * 1.05

  # Pseudo-log scale: linear within ±sigma, log farther out — handles
  # negatives, doesn't crush small values into a cluster at origin.
  # Pick sigma at the typical near-zero noise floor for the gas.
  sigma <- if (gas == "CH4") 0.1 else if (gas == "CO2") 0.1 else 0.05
  pslog <- scales::pseudo_log_trans(sigma = sigma, base = 10)

  # Choose break points that make sense on a pseudolog axis
  brks <- c(-1000, -100, -10, -1, 0, 1, 10, 100, 1000)
  brks <- brks[abs(brks) <= lim]

  annot <- tibble(
    x = c(-lim * 0.5,  lim * 0.5),
    y = c( lim * 0.5, -lim * 0.5),
    label = c("above 1:1\n=> missing source\n(canopy adds)",
              "below 1:1\n=> missing sink\n(canopy removes)")
  )

  # Rename logical levels to descriptive strings (for a clean legend)
  one <- one %>%
    mutate(shape_cat = factor(
      ifelse(within_canopy,
             "bottom sensor within canopy (local source/sink possible)",
             "bottom sensor above canopy (air column only — footprint-driven)"),
      levels = c("bottom sensor within canopy (local source/sink possible)",
                 "bottom sensor above canopy (air column only — footprint-driven)")
    ))

  ggplot(one, aes(x = F_bot, y = F_top,
                   color = canopy_class, shape = shape_cat)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    geom_point(size = 3, alpha = 0.9, stroke = 1.2) +
    ggrepel::geom_text_repel(aes(label = site),
                              size = 2.8, max.overlaps = 60,
                              box.padding = 0.25, segment.alpha = 0.4,
                              show.legend = FALSE) +
    geom_text(data = annot, aes(x = x, y = y, label = label),
              inherit.aes = FALSE,
              color = "grey45", size = 3, lineheight = 0.95, hjust = 0.5) +
    scale_color_manual(values = CANOPY_COLORS, name = "Canopy class") +
    scale_shape_manual(
      values = c(
        "bottom sensor within canopy (local source/sink possible)"    = 16,  # solid
        "bottom sensor above canopy (air column only — footprint-driven)" = 1   # open
      ),
      name = NULL, drop = FALSE) +
    scale_x_continuous(transform = pslog, breaks = brks,
                        limits = c(-lim, lim)) +
    scale_y_continuous(transform = pslog, breaks = brks,
                        limits = c(-lim, lim)) +
    coord_fixed() +
    labs(
      x = paste0("Ground-level F  (bottom sensor, \u2248 soil)  (", units,
                  ")  [pseudo-log]"),
      y = paste0("Ecosystem F  (top pair, \u2248 atmosphere exchange)  (",
                  units, ")  [pseudo-log]"),
      title = paste0(gas, " Ground vs Ecosystem Flux  (",
                     cond_label, ", per-site ", stat, ")"),
      subtitle = paste0("Pseudo-log axes (linear within \u00b1", sigma,
                         "; log outside). Dashed = 1:1.")
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top",
          legend.box = "vertical")
}
