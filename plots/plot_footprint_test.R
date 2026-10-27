# plot_footprint_test.R
# Figures for the pair-resolved footprint mechanism-2 test.
#
# Headline figures:
#   plot_dc_vs_dfwet()        — scatter of all pairs, colored by canopy class,
#                                with continental regression line
#   plot_dc_vs_dArea()        — same for ΔA (no NWI dependence)
#   plot_per_site_slopes()    — forest plot of per-site dfwet slopes ± SE,
#                                ordered by canopy class
#   plot_mechanism_map()      — map of CONUS+AK+HI+PR with classification dots
#   plot_footprint_geometry() — diagnostic: one example timestamp per site
#                                showing footprints overlaid on NWI polygons
#   plot_fwet_vs_height()     — per-site f_wet vs height (climatology) to show
#                                the monotonic footprint-expansion effect

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(sf)
  if (requireNamespace("scales", quietly = TRUE)) library(scales)
  if (requireNamespace("patchwork", quietly = TRUE)) library(patchwork)
})


#' Scatter: dC vs Δf_wet, colored by canopy class, with continental fit.
#'
#' @param pair_df Tibble. Pair metrics with site, dConc, dfwet columns.
#' @param canopy_metrics Tibble (optional). Adds canopy_class for color.
#' @param continental Single-row tibble from fit_continental_regression().
plot_dc_vs_dfwet <- function(pair_df, canopy_metrics = NULL,
                              continental = NULL, gas = "CH4") {

  df <- pair_df
  if (!is.null(canopy_metrics)) {
    df <- df %>% left_join(canopy_metrics %>% select(site, canopy_class),
                            by = "site")
  } else {
    df$canopy_class <- "all"
  }

  p <- ggplot(df, aes(x = dfwet, y = dConc, color = canopy_class)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = 2) +
    geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                aes(group = 1), formula = y ~ x) +
    scale_color_viridis_d(option = "viridis", end = 0.9) +
    labs(
      x = expression(Delta * italic(f)[wet] ~ "(upper - lower, footprint-weighted)"),
      y = bquote(Delta * .(gas) ~ "(upper - lower)"),
      color = "canopy class",
      title = paste0("Pair-resolved test: ", gas,
                     " concentration difference vs footprint wetness difference"),
      subtitle = if (!is.null(continental))
        sprintf("Continental slope = %.4g ± %.2g  (t=%.1f, n_pairs=%d, n_sites=%d)",
                continental$slope_dfwet, continental$se_dfwet,
                continental$t_dfwet, continental$n_total, continental$n_sites)
        else NULL
    ) +
    theme_minimal()
  p
}


#' Scatter: dC vs ΔA. No NWI needed; tests the pure-geometry version of
#' mechanism 2 (larger upper footprint -> more chance of intercepting anything
#' heterogeneous, wet or otherwise).
plot_dc_vs_dArea <- function(pair_df, canopy_metrics = NULL,
                              gas = "CH4") {

  df <- pair_df
  if (!is.null(canopy_metrics)) {
    df <- df %>% left_join(canopy_metrics %>% select(site, canopy_class),
                            by = "site")
  } else df$canopy_class <- "all"

  ggplot(df, aes(x = dArea, y = dConc, color = canopy_class)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = 2) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                aes(group = 1), formula = y ~ x) +
    scale_x_log10(labels = scales::label_number(scale_cut = scales::cut_si("m^2"))) +
    scale_color_viridis_d(option = "viridis", end = 0.9) +
    labs(
      x = expression(Delta * italic(A) ~ "(upper - lower, m"^2*")"),
      y = bquote(Delta * .(gas) ~ "(upper - lower)"),
      color = "canopy class",
      title = paste0(gas, " concentration vs footprint area expansion")
    ) +
    theme_minimal()
}


#' Forest plot of per-site dfwet slopes ± SE.
plot_per_site_slopes <- function(per_site_df, canopy_metrics = NULL,
                                   gas = "CH4") {

  df <- per_site_df %>%
    filter(!is.na(slope_dfwet)) %>%
    arrange(slope_dfwet)

  if (!is.null(canopy_metrics)) {
    df <- df %>% left_join(canopy_metrics %>% select(site, canopy_class),
                            by = "site")
  } else df$canopy_class <- "all"

  df$site <- factor(df$site, levels = df$site)

  ggplot(df, aes(y = site, x = slope_dfwet, color = canopy_class)) +
    geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
    geom_linerange(aes(xmin = slope_dfwet - se_dfwet,
                        xmax = slope_dfwet + se_dfwet)) +
    geom_point() +
    scale_color_viridis_d(option = "viridis", end = 0.9) +
    labs(
      x = bquote("per-site slope of " * Delta * .(gas) * " ~ " * Delta * italic(f)[wet]),
      y = "site (ordered by slope)",
      color = "canopy class",
      title = paste0("Per-site mechanism-2 slope (", gas, ")")
    ) +
    theme_minimal()
}


#' Map of sites colored by mechanism classification.
#'
#' Uses tower lat/lon from canopy_metrics or pair_df attributes.
#'
#' @param classified_df Output of classify_mechanism().
#' @param site_coords Tibble with site, lat, lon (e.g., from nwi inventory).
plot_mechanism_map <- function(classified_df, site_coords, gas = "CH4") {

  df <- classified_df %>% left_join(site_coords, by = "site") %>%
    filter(!is.na(lat), !is.na(lon))

  pal <- c(
    "mechanism 2 active" = "#1f77b4",
    "mechanism 1 / closed / no wetness" = "#2ca02c",
    "negative dfwet (anomalous)" = "#d62728",
    "untested" = "grey60"
  )

  ggplot(df, aes(x = lon, y = lat, color = mech_class)) +
    geom_point(size = 3) +
    scale_color_manual(values = pal) +
    coord_quickmap() +
    labs(
      x = "longitude", y = "latitude",
      color = "classification",
      title = paste0("Mechanism classification by site (", gas, ")")
    ) +
    theme_minimal()
}


#' Per-site f_wet vs height climatology.
#'
#' Shows how the footprint-weighted wetland fraction grows with measurement
#' height — the *signal* that mechanism-2 acts on.
#'
#' @param pair_df Long-form pair metrics with z_A, z_B, f_wet_A, f_wet_B
plot_fwet_vs_height <- function(pair_df, canopy_metrics = NULL,
                                  facet_by_class = TRUE) {

  # Unfold the pair structure into per-(site, z, f_wet) long form. Doing this
  # explicitly avoids pivot_longer's interaction with a site column.
  long <- bind_rows(
    pair_df %>% select(site, z_m = z_A, f_wet = f_wet_A),
    pair_df %>% select(site, z_m = z_B, f_wet = f_wet_B)
  ) %>%
    filter(!is.na(z_m), !is.na(f_wet))

  if (!is.null(canopy_metrics)) {
    long <- long %>% left_join(canopy_metrics %>%
                                  select(site, canopy_class, canopy_height),
                                by = "site")
  }

  group_cols <- c("site", "z_m")
  if ("canopy_class" %in% names(long)) group_cols <- c(group_cols, "canopy_class")
  summ <- long %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(f_wet_mean = mean(f_wet, na.rm = TRUE),
              .groups = "drop")

  p <- ggplot(summ, aes(x = z_m, y = f_wet_mean, color = site, group = site)) +
    geom_line() + geom_point() +
    scale_color_discrete(guide = "none") +
    labs(x = "measurement height (m)",
         y = "footprint-weighted wetland fraction",
         title = "f_wet vs height (climatology) — mechanism-2 signal") +
    theme_minimal()
  if (facet_by_class && "canopy_class" %in% names(summ)) {
    p <- p + facet_wrap(~ canopy_class, scales = "free_y")
  }
  p
}


#' Diagnostic: plot one example timestamp per site showing the 80% footprint
#' contour at each height, overlaid on the NWI polygons. Used to sanity-check
#' that the footprint code is producing sensible shapes and that the NWI
#' overlap is plausible.
plot_footprint_geometry <- function(site, attr_df, nwi_proj, example_row,
                                      heights, tower_positions,
                                      r_pct = 80) {

  lat <- as.numeric(attr_df$LatTow[1])
  lon <- as.numeric(attr_df$LonTow[1])
  target_crs <- pick_local_aea(lon, lat)

  fps <- compute_footprints_all_heights(
    met_row = example_row, heights = heights,
    tower_positions = tower_positions, lat = lat, r = r_pct
  )

  wind_dir <- if ("wind_dir" %in% names(example_row))
    example_row$wind_dir[1] else NA_real_

  contours <- list()
  for (i in seq_len(nrow(fps))) {
    ffp <- fps$footprint[[i]]
    poly_sf <- footprint_contour_to_sf(ffp, wind_dir, lon, lat, target_crs)
    if (!is.null(poly_sf)) {
      poly_sf$z_m <- fps$z_m[i]
      contours[[i]] <- poly_sf
    }
  }
  contours_sf <- if (length(contours) > 0) do.call(rbind, contours) else NULL

  p <- ggplot() +
    {if (nrow(nwi_proj) > 0)
        geom_sf(data = nwi_proj, fill = "steelblue", alpha = 0.4,
                color = NA)} +
    {if (!is.null(contours_sf))
        geom_sf(data = contours_sf, aes(color = z_m), fill = NA, size = 0.6)} +
    geom_point(aes(x = 0, y = 0), color = "red", size = 2) +
    scale_color_viridis_c(option = "viridis") +
    labs(title = paste0(site, ": 80% footprint contours by height (example)"),
         color = "height (m)") +
    theme_minimal()
  p
}
