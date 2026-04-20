# plot_method_workflow.R
# Publication-ready conceptual methods-workflow figure for a single site.
#
# Shows the full flux computation stack in a clean, narrative layout:
#
#   Row 1 (observations):
#     1. Tower / canopy schematic — sensors and adjacent-pair brackets
#     2. Concentration profile     — with ΔC labels between adjacent sensors
#     3. Atmosphere anomaly         — C deviation from above-canopy reference
#     4. Wind speed profile
#
#   Row 2 (derived / decomposition):
#     5. Eddy diffusivity K(z)
#     6. Vertical flux profile F(z) — each pair labeled; top pair highlighted
#     7. Additive cascade           — how the layers sum from soil up to top
#     8. Soil / vegetation / ecosystem bars
#
# Muted natural palette intended for publication-quality figures.


# ---- Shared muted-natural palette ----
MW_PAL <- list(
  ground       = "#8C6A4A",  # soil brown (only for the soil segment)
  canopy       = "#4F7942",  # forest sage (canopy source)
  canopy_lt    = "#7FA876",  # lighter canopy (canopy source)
  canopy_sink  = "#A86C5C",  # muted clay (canopy sink — photosynthetic uptake)
  atmos        = "#8FA5C4",  # dusty slate blue
  atmos_lt     = "#C9D5E3",  # very light atmos
  conc         = "#2F4858",  # deep slate (measurements)
  wind         = "#B59E5F",  # muted gold
  K            = "#6B5B95",  # muted plum
  flux         = "#B55354",  # brick red (default flux line)
  flux_top     = "#3E2C4F",  # near-black aubergine (whole-ecosystem highlight)
  sensor       = "#2F4858",
  neutral      = "#5C5C5C",
  grid_light   = "#E6E2D3"
)


# Minimal unified theme
theme_mw <- function(base_size = 11) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 1,
                                    hjust = 0, color = "grey15",
                                    margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 2, color = "grey40",
                                    hjust = 0, margin = margin(b = 6)),
      panel.grid.major = element_line(color = MW_PAL$grid_light, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line   = element_line(color = "grey60", linewidth = 0.3),
      axis.ticks  = element_line(color = "grey60", linewidth = 0.3),
      axis.title  = element_text(color = "grey25", size = base_size - 1),
      axis.text   = element_text(color = "grey35", size = base_size - 2),
      plot.margin = margin(4, 8, 4, 4)
    )
}


plot_method_workflow <- function(site, gas,
                                   all_profiles, all_met, all_winds,
                                   all_fluxes, all_data,
                                   filter_jja_night = TRUE,
                                   tod = c("Night", "Day"),
                                   source_sink_df = NULL) {
  tod <- match.arg(tod)

  requireNamespace("patchwork", quietly = TRUE)

  attr_df <- all_data[[site]]$attr
  prof    <- all_profiles[[site]][[gas]]
  met     <- all_met[[site]]
  winds   <- all_winds[[site]]
  flux    <- all_fluxes[[site]][[gas]]
  if (is.null(prof) || is.null(flux) || is.null(attr_df)) {
    stop("Missing data for ", site, " / ", gas)
  }

  height_lookup <- get_height_lookup(attr_df)
  canopy_h <- get_canopy_height(attr_df)
  z_max <- max(height_lookup$height_m) * 1.08
  cond_label <- if (filter_jja_night) paste0("JJA ", tod) else "All data"

  # ---- Apply season/diel filters to the observational panels ----
  if (filter_jja_night) {
    sel_tod <- tod   # avoid shadowing below
    prof <- prof %>% add_season("timeMid") %>% add_tod("timeMid") %>%
      filter(season == "JJA", as.character(tod) == sel_tod)
    flux <- flux %>%
      mutate(month  = lubridate::month(time_round),
             hour   = lubridate::hour(time_round),
             season = case_when(month %in% 6:8 ~ "JJA",
                                month %in% 3:5 ~ "MAM",
                                month %in% 9:11 ~ "SON",
                                TRUE ~ "DJF"),
             tod_ch = ifelse(hour >= 6 & hour < 18, "Day", "Night")) %>%
      filter(season == "JJA", tod_ch == sel_tod)
  }

  scale_gas <- flux_scale_for_gas(gas)
  units_gas <- flux_units_for_gas(gas)
  conc_unit <- switch(gas, "CH4" = "ppm", "CO2" = "ppm", "H2O" = "mmol/mol", "")

  sensor_z <- sort(height_lookup$height_m)
  n_sensors <- length(sensor_z)

  # ===========================================================================
  # PANEL 1 — tower & canopy schematic with pair brackets
  # ===========================================================================
  #
  # A stylized single-tower column: brown ground slab, green canopy band (0 to
  # canopy_h), blue atmosphere above. Sensors as dots. Small brackets to the
  # right of the column indicate each adjacent-sensor pair (F_i). The
  # topmost pair bracket is drawn darker — it is the "whole-ecosystem" flux.

  x_col  <- 0
  x_pair <- 0.9       # where pair brackets sit
  pair_df <- tibble(
    pair_idx  = seq_len(n_sensors - 1),
    z_lower   = sensor_z[-n_sensors],
    z_upper   = sensor_z[-1],
    z_mid     = (sensor_z[-n_sensors] + sensor_z[-1]) / 2,
    label     = paste0("F", seq_len(n_sensors - 1)),
    is_top    = c(rep(FALSE, n_sensors - 2), TRUE)
  )

  p1 <- ggplot() +
    # atmosphere
    annotate("rect", xmin = -0.6, xmax = 0.6,
             ymin = max(canopy_h, 0.1), ymax = z_max,
             fill = MW_PAL$atmos_lt, alpha = 0.5) +
    # canopy band
    annotate("rect", xmin = -0.6, xmax = 0.6, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.5) +
    # canopy top
    geom_hline(yintercept = canopy_h, color = MW_PAL$canopy,
               linetype = "dashed", linewidth = 0.5) +
    # tower column
    annotate("segment", x = x_col, xend = x_col,
             y = 0, yend = z_max, color = MW_PAL$neutral, linewidth = 0.6) +
    # sensors
    annotate("point", x = rep(x_col, n_sensors), y = sensor_z,
             color = MW_PAL$sensor, size = 3.2) +
    # sensor labels
    annotate("text", x = x_col - 0.15, y = sensor_z,
             label = paste0(round(sensor_z, 1), " m"),
             hjust = 1, size = 3, color = "grey30") +
    # pair brackets (right side)
    geom_segment(data = pair_df,
                 aes(x = x_pair, xend = x_pair,
                     y = z_lower, yend = z_upper,
                     color = is_top, linewidth = is_top)) +
    geom_segment(data = pair_df,
                 aes(x = x_pair - 0.04, xend = x_pair,
                     y = z_lower, yend = z_lower,
                     color = is_top, linewidth = is_top)) +
    geom_segment(data = pair_df,
                 aes(x = x_pair - 0.04, xend = x_pair,
                     y = z_upper, yend = z_upper,
                     color = is_top, linewidth = is_top)) +
    geom_text(data = pair_df,
              aes(x = x_pair + 0.08, y = z_mid,
                  label = label, color = is_top),
              hjust = 0, size = 3.2, fontface = "bold") +
    # canopy label
    annotate("text", x = -0.55, y = canopy_h / 2, label = "canopy",
             angle = 90, color = MW_PAL$canopy, size = 3, hjust = 0.5) +
    annotate("text", x = -0.55, y = (canopy_h + z_max)/2, label = "atmosphere",
             angle = 90, color = MW_PAL$atmos, size = 3, hjust = 0.5) +
    scale_color_manual(values = c(`FALSE` = MW_PAL$neutral,
                                    `TRUE`  = MW_PAL$flux_top), guide = "none") +
    scale_linewidth_manual(values = c(`FALSE` = 0.4, `TRUE` = 1.2),
                            guide = "none") +
    scale_x_continuous(limits = c(-0.9, 1.4), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = NULL, y = "Height (m)",
         title = "1. Tower & sensor pairs",
         subtitle = paste0(site, "  ·  canopy ", round(canopy_h, 1), " m",
                            "  ·  ", n_sensors, " sensors, ", n_sensors - 1,
                            " adjacent pairs")) +
    theme_mw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank())

  # ===========================================================================
  # PANEL 2 — Concentration profile with ΔC annotations
  # ===========================================================================
  conc_stats <- prof %>%
    inner_join(height_lookup, by = "TowerPosition") %>%
    group_by(height_m) %>%
    summarise(C_med = median(concentration, na.rm = TRUE),
              C_q25 = quantile(concentration, 0.25, na.rm = TRUE),
              C_q75 = quantile(concentration, 0.75, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(height_m)

  dC_ann <- conc_stats %>%
    mutate(C_next = dplyr::lead(C_med),
            z_next = dplyr::lead(height_m),
            dC     = C_next - C_med,
            z_mid  = (height_m + z_next) / 2,
            C_mid  = (C_med + C_next) / 2) %>%
    filter(!is.na(dC))

  p2 <- ggplot(conc_stats, aes(x = C_med, y = height_m)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.35) +
    geom_errorbarh(aes(xmin = C_q25, xmax = C_q75),
                   height = 0, linewidth = 0.4, color = "grey55") +
    geom_line(orientation = "y", linewidth = 0.6, color = MW_PAL$conc) +
    geom_point(size = 2.8, color = MW_PAL$conc) +
    geom_label(data = dC_ann,
               aes(x = C_mid, y = z_mid,
                    label = sprintf("\u0394C = %+.4f", dC)),
               size = 2.5, color = MW_PAL$flux, label.size = 0,
               label.padding = unit(0.08, "lines"),
               fill = alpha("white", 0.85)) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = paste0(gas, " (", conc_unit, ")"), y = NULL,
         title = paste0("2. ", gas, " profile"),
         subtitle = "ΔC between adjacent sensors") +
    theme_mw()

  # ===========================================================================
  # PANEL 3 — Atmosphere anomaly profile
  # ===========================================================================
  top_ref <- conc_stats$C_med[which.max(conc_stats$height_m)]
  anom <- conc_stats %>% mutate(anom = C_med - top_ref,
                                  anom_lo = C_q25 - top_ref,
                                  anom_hi = C_q75 - top_ref)

  p3 <- ggplot(anom, aes(x = anom, y = height_m)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.35) +
    geom_vline(xintercept = 0, color = MW_PAL$atmos, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = anom_lo, xmax = anom_hi),
                   height = 0, linewidth = 0.4, color = "grey55") +
    geom_line(orientation = "y", linewidth = 0.6, color = MW_PAL$atmos) +
    geom_point(size = 2.8, color = MW_PAL$atmos) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = paste0("C − C(top)  (", conc_unit, ")"), y = NULL,
         title = "3. Anomaly from atmosphere",
         subtitle = "Reference = above-canopy sensor") +
    theme_mw()

  # ===========================================================================
  # PANEL 4 — Wind profile
  # ===========================================================================
  wind_stats <- tibble()
  if (!is.null(winds) && nrow(winds) > 0) {
    w <- winds
    if (filter_jja_night && "timeMid" %in% names(w)) {
      sel_tod <- tod
      w <- w %>% add_season("timeMid") %>% add_tod("timeMid") %>%
        filter(season == "JJA", as.character(tod) == sel_tod)
    }
    wind_stats <- w %>%
      group_by(height_m) %>%
      summarise(u_med = median(ubar, na.rm = TRUE),
                u_q25 = quantile(ubar, 0.25, na.rm = TRUE),
                u_q75 = quantile(ubar, 0.75, na.rm = TRUE),
                .groups = "drop") %>%
      filter(!is.na(u_med)) %>% arrange(height_m)
  }

  p4 <- ggplot(wind_stats, aes(x = u_med, y = height_m)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.35) +
    geom_errorbarh(aes(xmin = u_q25, xmax = u_q75),
                   height = 0, linewidth = 0.4, color = "grey55") +
    geom_line(orientation = "y", linewidth = 0.6, color = MW_PAL$wind) +
    geom_point(size = 2.8, color = MW_PAL$wind) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = "u (m/s)", y = NULL,
         title = "4. Wind speed",
         subtitle = "per sensor height") +
    theme_mw()

  # ===========================================================================
  # PANEL 5 — K(z) profile
  # ===========================================================================
  K_stats <- flux %>%
    group_by(z) %>%
    summarise(K_med = median(K, na.rm = TRUE),
              K_q25 = quantile(K, 0.25, na.rm = TRUE),
              K_q75 = quantile(K, 0.75, na.rm = TRUE),
              .groups = "drop") %>% arrange(z)

  p5 <- ggplot(K_stats, aes(x = K_med, y = z)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.35) +
    geom_errorbarh(aes(xmin = K_q25, xmax = K_q75),
                   height = 0, linewidth = 0.4, color = "grey55") +
    geom_line(orientation = "y", linewidth = 0.6, color = MW_PAL$K) +
    geom_point(size = 2.8, color = MW_PAL$K) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = expression(K~(m^2/s)), y = "Height (m)",
         title = "5. Turbulent diffusivity",
         subtitle = "MOST above canopy + Raupach below") +
    theme_mw()

  # ===========================================================================
  # PANEL 6 — Vertical flux profile (pairs labeled, top highlighted)
  # ===========================================================================
  F_stats <- flux %>%
    group_by(z) %>%
    summarise(F_med = median(F_FG, na.rm = TRUE) * scale_gas,
              F_q25 = quantile(F_FG, 0.25, na.rm = TRUE) * scale_gas,
              F_q75 = quantile(F_FG, 0.75, na.rm = TRUE) * scale_gas,
              .groups = "drop") %>% arrange(z) %>%
    mutate(is_top = z == max(z),
            pair  = paste0("F", seq_len(n())))

  # Separate the top-pair for an explicit "ecosystem" annotation
  F_stats <- F_stats %>%
    mutate(ring = is_top)  # use ring aesthetic for top

  p6 <- ggplot(F_stats, aes(x = F_med, y = z)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = max(canopy_h, 0.1),
             fill = MW_PAL$canopy_lt, alpha = 0.35) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = F_q25, xmax = F_q75),
                   height = 0, linewidth = 0.4, color = "grey55") +
    geom_line(orientation = "y", linewidth = 0.6, color = MW_PAL$flux) +
    # ordinary points
    geom_point(data = ~ subset(.x, !is_top),
               size = 2.8, color = MW_PAL$flux) +
    # top (ecosystem) point — ring + filled
    geom_point(data = ~ subset(.x, is_top),
               size = 5.5, shape = 21,
               color = MW_PAL$flux_top, fill = MW_PAL$flux_top,
               stroke = 1.2) +
    geom_point(data = ~ subset(.x, is_top),
               size = 8, shape = 1, color = MW_PAL$flux_top, stroke = 0.6) +
    # single label for the ecosystem point only
    geom_text(data = ~ subset(.x, is_top),
              aes(label = "ecosystem flux"),
              hjust = -0.2, nudge_x = 0, nudge_y = z_max * 0.04,
              size = 3, fontface = "bold", color = MW_PAL$flux_top) +
    scale_y_continuous(limits = c(0, z_max), expand = c(0, 0)) +
    labs(x = paste0("F  (", units_gas, ")"), y = NULL,
         title = "6. Flux profile",
         subtitle = "F(z) = −K(z) · dC/dz at each sensor pair") +
    theme_mw()

  # ===========================================================================
  # PANEL 7 — Lagrangian flux decomposition (horizontal, bottom-up)
  # ===========================================================================
  # Each row is an INDEPENDENT contribution drawn as a vector from zero:
  # soil flux (brown), then each canopy layer's S_net (dF + storage, when
  # source_sink_df is provided) or dF alone. Hollow tick marks show dF-only
  # so the storage contribution is visible. The final aubergine row is the
  # direct top-pair measurement (independent of the layer components).
  F_sorted <- F_stats %>% arrange(z)
  z_levels <- F_sorted$z
  soil_val <- F_sorted$F_med[1]         # F at bottom pair
  eco_val  <- F_sorted$F_med[nrow(F_sorted)]  # F at top pair

  # Layer source/sink values — prefer S_net from source_sink_df when available
  use_ss <- FALSE
  if (!is.null(source_sink_df) && nrow(source_sink_df) > 0) {
    ss_filt <- source_sink_df
    if (filter_jja_night) {
      ss_filt <- ss_filt %>%
        mutate(.mon  = lubridate::month(time_round),
               .hour = lubridate::hour(time_round)) %>%
        filter(.mon %in% 6:8,
               if (tod == "Night") (.hour < 6 | .hour >= 18) else (.hour >= 6 & .hour < 18))
    }
    if (nrow(ss_filt) > 0) {
      ss_summary <- ss_filt %>%
        group_by(layer_label, z_lower, z_upper) %>%
        summarise(S_net  = mean(S_net,          na.rm = TRUE) * scale_gas,
                  dF_val = mean(S_net_nostorage, na.rm = TRUE) * scale_gas,
                  .groups = "drop") %>%
        arrange(z_lower)
      if (nrow(ss_summary) > 0) {
        use_ss     <- TRUE
        layer_zlo  <- ss_summary$z_lower
        layer_zhi  <- ss_summary$z_upper
        layer_vals <- ss_summary$S_net
        layer_dF   <- ss_summary$dF_val
        p7_subtitle <- "soil + layer S_net (dF + storage)  |  hollow tick = dF alone"
      }
    }
  }
  if (!use_ss) {
    layer_vals <- diff(F_sorted$F_med)
    layer_dF   <- layer_vals
    layer_zlo  <- z_levels[-length(z_levels)]
    layer_zhi  <- z_levels[-1]
    p7_subtitle <- "soil + canopy layer \u0394F = ecosystem flux (by construction)"
  }

  component_rows <- bind_rows(
    tibble(step     = 1L,
           label    = paste0("soil (", round(z_levels[1], 1), " m)"),
           value    = soil_val,
           value_dF = soil_val,
           fill     = MW_PAL$ground,
           is_sum   = FALSE),
    if (length(layer_vals) > 0)
      tibble(step     = seq(2L, length.out = length(layer_vals)),
             label    = sprintf("%.1f \u2013 %.1f m", layer_zlo, layer_zhi),
             value    = layer_vals,
             value_dF = layer_dF,
             fill     = MW_PAL$canopy_lt,
             is_sum   = FALSE)
    else
      tibble()
  )
  # When storage is included, sum of components ≠ top pair; label accordingly.
  sum_label <- if (use_ss) "ecosystem\n(top pair)" else "\u03a3 = ecosystem\n(top pair)"
  sum_row <- tibble(
    step     = nrow(component_rows) + 1L,
    label    = sum_label,
    value    = eco_val,
    value_dF = eco_val,
    fill     = MW_PAL$flux_top,
    is_sum   = TRUE
  )
  casc <- bind_rows(component_rows, sum_row) %>% mutate(ypos = seq_len(n()))
  y_divider <- nrow(component_rows) + 0.5

  xr <- range(c(0, casc$value, casc$value_dF))
  xpad <- 0.08 * diff(xr)
  if (xpad == 0) xpad <- 0.1

  p7 <- ggplot(casc, aes(y = ypos)) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    geom_hline(yintercept = y_divider, linetype = "dashed",
               color = "grey55", linewidth = 0.35) +
    # primary arrows: S_net (or dF when no source_sink_df)
    geom_segment(aes(x = 0, xend = value, yend = ypos,
                     color = fill, linewidth = is_sum),
                 lineend = "butt",
                 arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    # hollow tick at dF-only position (visible when storage differs)
    geom_point(aes(x = value_dF), shape = 124, size = 4, color = "grey40") +
    # signed value labels
    geom_text(aes(x = ifelse(value >= 0, value, 0) + xpad * 0.25,
                   label = sprintf("%+.3f", value),
                   fontface = ifelse(is_sum, "bold", "plain")),
              hjust = 0, size = 3, color = "grey15") +
    # left-gutter labels
    geom_text(aes(x = min(xr) - xpad * 1.1, label = label,
                   fontface = ifelse(is_sum, "bold", "plain")),
              hjust = 1, size = 3, color = "grey25", lineheight = 0.9) +
    scale_color_identity() +
    scale_linewidth_manual(values = c(`FALSE` = 5, `TRUE` = 7), guide = "none") +
    scale_y_continuous(limits = c(0.4, nrow(casc) + 0.6)) +
    scale_x_continuous(limits = c(min(xr) - 2 * xpad,
                                    max(xr) + xpad * 1.4)) +
    labs(x = paste0("F contribution  (", units_gas, ")"), y = NULL,
         title = "7. Lagrangian flux decomposition",
         subtitle = p7_subtitle) +
    theme_mw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank())

  # ===========================================================================
  # PANEL 8 — Soil / vegetation / ecosystem; closure: direct vs bottom-up
  #
  # We show ecosystem flux TWO ways using DISJOINT sensor data, which is the
  # non-trivial Lagrangian closure check:
  #   - "top pair (direct)"     = F at sensor N-1, uses top two sensors only
  #   - "bottom-up (Σ below)"   = F at sensor N-2, uses sensors 1..N-1 only
  #                               (= soil + Σ all canopy layers except topmost,
  #                                via telescoping — no shared data with top)
  # Agreement between the two = the K(z) model is self-consistent.
  # ===========================================================================
  F_sorted <- F_stats %>% arrange(z)
  if (nrow(F_sorted) < 2) {
    p8 <- ggplot() + theme_void() + labs(title = "8. Decomposition")
  } else {
    n_pairs  <- nrow(F_sorted)   # = N-1 pair fluxes
    soil_val <- F_sorted$F_med[1]                       # F_1
    bu_val   <- F_sorted$F_med[n_pairs - 1]             # F_{n-1} (bottom-up eco)
    top_val  <- F_sorted$F_med[n_pairs]                 # F_n (top pair, eco direct)
    veg_val  <- bu_val - soil_val                        # F_{n-1} - F_1 = vegetation

    bars <- tibble(
      what = factor(c(
        "soil\n(F\u2081)",
        "vegetation\n(F\u2099\u208B\u2081 \u2212 F\u2081)",
        "ecosystem\nbottom-up\n(F\u2099\u208B\u2081)",
        "ecosystem\ntop pair\n(F\u2099)"),
        levels = c(
          "soil\n(F\u2081)",
          "vegetation\n(F\u2099\u208B\u2081 \u2212 F\u2081)",
          "ecosystem\nbottom-up\n(F\u2099\u208B\u2081)",
          "ecosystem\ntop pair\n(F\u2099)")),
      value   = c(soil_val, veg_val, bu_val, top_val),
      fillcol = c(MW_PAL$ground, MW_PAL$canopy,
                   MW_PAL$flux_top, MW_PAL$flux_top),
      alpha_v = c(0.9, 0.9, 0.5, 0.9)    # bottom-up lighter so direct reads as primary
    )

    closure_delta <- top_val - bu_val
    closure_note  <- paste0(
      "soil + vegetation = ecosystem bottom-up  \u00b7  ",
      sprintf("closure: top pair \u2212 bottom-up = %+.3f %s",
              closure_delta, units_gas))

    p8 <- ggplot(bars, aes(x = what, y = value,
                            fill = fillcol, alpha = alpha_v)) +
      geom_hline(yintercept = 0, color = "grey45", linewidth = 0.4) +
      geom_vline(xintercept = 2.5, color = "grey70",
                 linetype = "dotted", linewidth = 0.35) +
      geom_col(width = 0.6) +
      geom_text(aes(label = sprintf("%+.3f", value),
                    vjust = ifelse(value >= 0, -0.4, 1.4)),
                size = 3.2, color = "grey15",
                alpha = 1, show.legend = FALSE) +
      scale_fill_identity() +
      scale_alpha_identity() +
      labs(x = NULL, y = paste0("F  (", units_gas, ")"),
           title = "8. Flux decomposition & closure",
           subtitle = closure_note) +
      theme_mw() +
      theme(axis.text.x = element_text(size = 8.5, color = "grey20",
                                          lineheight = 0.9))
  }

  # ===========================================================================
  # Assemble
  # ===========================================================================
  layout <- (p1 | p2 | p3 | p4) /
            (p5 | p6 | p7 | p8) +
    patchwork::plot_annotation(
      title    = paste0("Flux-gradient method walkthrough \u2014 ",
                         site, ", ", gas, "  (", cond_label, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = 15,
                                                color = "grey10"))
    )
  layout
}
