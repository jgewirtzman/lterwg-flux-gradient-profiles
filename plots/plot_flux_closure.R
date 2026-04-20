# plot_flux_closure.R
# Diagnostic: FG-computed top flux vs EC-measured top flux (CO2/H2O only)

#' Scatter of FG top flux vs EC top flux per site
#'
#' @param layer_fluxes_all Named list (by site) of layer flux tibbles for a single gas
#' @param canopy_metrics Tibble. Canopy metrics (for site ordering/coloring)
#' @param gas Character. "CO2" or "H2O"
#' @param subset_condition Optional — tibble with filter conditions
#' @return ggplot object
plot_flux_closure <- function(layer_fluxes_all, canopy_metrics, gas = "CO2") {

  # Extract the topmost sensor's flux from each site. Under Option A the
  # topmost flux lives at z = kept_h[N-1] (pair between second-from-top and
  # top sensor). z_upper = top sensor height → largest z_upper identifies it.
  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    top_pair <- lf %>% filter(z_upper == max(z_upper, na.rm = TRUE))
    if (nrow(top_pair) == 0) next
    rows[[site]] <- top_pair %>% mutate(site = site)
  }
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    filter(!is.na(F_FG), !is.na(F_EC)) %>%
    # Remove extreme outliers for plotting
    filter(abs(F_FG) < 1e4, abs(F_EC) < 1e4)

  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No paired data"))

  # Add canopy metrics
  df <- df %>%
    left_join(canopy_metrics %>% select(site, canopy_class), by = "site")

  # Compute fit statistics per site
  fit_stats <- df %>%
    group_by(site) %>%
    summarise(
      n = n(),
      r = cor(F_FG, F_EC, use = "complete.obs"),
      slope = coef(lm(F_FG ~ F_EC))[2],
      rmse = sqrt(mean((F_FG - F_EC)^2, na.rm = TRUE)),
      .groups = "drop"
    )

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  df <- df %>% mutate(F_FG = F_FG * scale, F_EC = F_EC * scale)

  ggplot(df, aes(x = F_EC, y = F_FG, color = canopy_class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black") +
    scale_color_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class") +
    facet_wrap(~ site, scales = "free", ncol = 6) +
    labs(
      x = paste0("F_EC (", units, ")"),
      y = paste0("F_FG top-pair (", units, ")"),
      title = paste0(gas, " Flux Closure: FG top pair vs Eddy Covariance"),
      subtitle = "Dashed line = 1:1. Good closure validates the K(z) model."
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "top"
    )
}


#' Lagrangian-style closure: above-canopy top-pair flux vs bottom-up
#' reconstruction from all lower pairs.
#'
#' On NEON towers with ≥3 sensors, the top two sensors are typically BOTH
#' above the canopy — so the flux computed from that top pair samples purely
#' turbulent transport with no biological source/sink in between. Meanwhile,
#' starting at the lowest sensor (≈ soil flux) and accumulating all the
#' within-canopy layer source/sinks gives an INDEPENDENT bottom-up
#' reconstruction of the ecosystem flux, using entirely different data.
#'
#' Under Option A (F labeled at the lower sensor of each pair):
#'   F_top_pair   = F at sensor (N-1)   — uses sensors (N-1, N), purely above-canopy
#'   F_recon      = F_1 + Σ(F_{i+1}-F_i) for i in 1..N-3
#'                = F at sensor (N-2)   — uses sensors 1..N-1, everything below
#' These two estimates use DISJOINT data (apart from sensor N-1 playing
#' different roles in each) and different K values, so their agreement is a
#' non-trivial closure check that works without any EC measurement —
#' critical for CH4 where no EC flux exists at all.
#'
#' Sites with fewer than 3 sensors are skipped (no closure possible).
#'
#' @param layer_fluxes_all Named list (by site) of layer flux tibbles for a single gas
#' @param canopy_metrics Tibble
#' @param gas Character
#' @return ggplot object
plot_flux_closure_lagrangian <- function(layer_fluxes_all, canopy_metrics,
                                          gas = "CH4") {

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next

    # Determine the two highest pairs by looking at z (lower-sensor of pair).
    # Under Option A, F is labeled at z; the top-pair flux sits at the
    # second-highest sensor, and the bottom-up reconstruction's endpoint
    # sits at the third-highest sensor.
    z_vals <- sort(unique(lf$z))
    if (length(z_vals) < 2) next   # need ≥3 sensors → ≥2 pair fluxes
    z_top    <- z_vals[length(z_vals)]       # F at sensor N-1 (top pair)
    z_recon  <- z_vals[length(z_vals) - 1]   # F at sensor N-2 (bottom-up)

    top_lf <- lf %>% filter(z == z_top) %>%
      select(time_round, F_top = F_FG, F_EC)
    rec_lf <- lf %>% filter(z == z_recon) %>%
      select(time_round, F_recon = F_FG)
    if (nrow(top_lf) == 0 || nrow(rec_lf) == 0) next

    paired <- inner_join(top_lf, rec_lf, by = "time_round") %>%
      mutate(site = site, z_top = z_top, z_recon = z_recon)
    rows[[site]] <- paired
  }
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  df <- df %>%
    filter(!is.na(F_top), !is.na(F_recon),
           abs(F_top) < 1e4, abs(F_recon) < 1e4) %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site")

  site_order <- canopy_metrics %>% arrange(canopy_height) %>% pull(site)
  df$site <- factor(df$site, levels = intersect(site_order, unique(df$site)))

  scale <- flux_scale_for_gas(gas)
  units <- flux_units_for_gas(gas)
  df <- df %>% mutate(F_top   = F_top   * scale,
                       F_recon = F_recon * scale,
                       F_EC    = F_EC    * scale)

  # Median F_EC per site for reference overlay (CO2/H2O only)
  ec_ref <- df %>%
    filter(!is.na(F_EC), abs(F_EC) < 1e4) %>%
    group_by(site) %>%
    summarise(F_EC_med = median(F_EC, na.rm = TRUE), .groups = "drop")

  p <- ggplot(df, aes(x = F_recon, y = F_top, color = canopy_class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "black") +
    scale_color_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class") +
    facet_wrap(~ site, scales = "free", ncol = 6) +
    labs(
      x = paste0("F_FG bottom-up reconstruction [= F at sensor N-2] (", units, ")"),
      y = paste0("F_FG top pair [= F at sensor N-1, above canopy] (", units, ")"),
      title = paste0(gas, " Lagrangian Closure: above-canopy pair vs below reconstruction"),
      subtitle = paste0(
        "Dashed = 1:1. Uses disjoint data — no EC required. ",
        if (gas %in% c("CO2","H2O")) "Red dashed = median F_EC (reference)." else ""
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "top"
    )

  if (nrow(ec_ref) > 0 && gas %in% c("CO2", "H2O")) {
    p <- p + geom_hline(data = ec_ref, aes(yintercept = F_EC_med),
                        color = "red", linetype = "dashed", linewidth = 0.4)
  }
  p
}


#' Lagrangian closure summary: median(F_top) / median(F_recon) and r per site,
#' where F_top = top-pair flux and F_recon = flux reconstructed from all lower
#' pairs (= F at sensor N-2 under Option A). Works for all gases — no EC needed.
plot_flux_closure_lagrangian_summary <- function(layer_fluxes_all, canopy_metrics,
                                                 gas = "CH4") {
  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    z_vals <- sort(unique(lf$z))
    if (length(z_vals) < 2) next
    z_top   <- z_vals[length(z_vals)]
    z_recon <- z_vals[length(z_vals) - 1]
    top_lf <- lf %>% filter(z == z_top)   %>% select(time_round, F_top   = F_FG)
    rec_lf <- lf %>% filter(z == z_recon) %>% select(time_round, F_recon = F_FG)
    paired <- inner_join(top_lf, rec_lf, by = "time_round") %>%
      filter(!is.na(F_top), !is.na(F_recon),
             abs(F_top) < 1e4, abs(F_recon) < 1e4)
    if (nrow(paired) < 20) next
    rows[[site]] <- tibble(site = site, n = nrow(paired),
                            r = cor(paired$F_top, paired$F_recon),
                            med_top = median(paired$F_top),
                            med_rec = median(paired$F_recon))
  }
  stats <- bind_rows(rows)
  if (nrow(stats) == 0) return(ggplot() + theme_void() + labs(title = "No data"))

  stats <- stats %>%
    mutate(ratio = med_top / med_rec,
           sign_agree = sign(med_top) == sign(med_rec)) %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site") %>%
    filter(!is.na(r))

  stats$site <- reorder(stats$site, stats$r)

  ggplot(stats, aes(x = site, y = r, fill = canopy_class)) +
    geom_col() +
    geom_hline(yintercept = 0,   linewidth = 0.3, color = "grey50") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class") +
    coord_flip() +
    labs(x = NULL, y = "Pearson r  (F_top-pair vs F_recon)",
         title = paste0(gas, " Lagrangian Closure r by Site (no EC needed)"),
         subtitle = "Dashed: r = 0.5 benchmark. All bars use purely internal FG closure.") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
}


#' Summary: r and slope per site across all sites
plot_flux_closure_summary <- function(layer_fluxes_all, canopy_metrics, gas = "CO2") {

  rows <- list()
  for (site in names(layer_fluxes_all)) {
    lf <- layer_fluxes_all[[site]]
    if (is.null(lf) || nrow(lf) == 0) next
    top_pair <- lf %>% filter(z_upper == max(z_upper, na.rm = TRUE))
    rows[[site]] <- top_pair %>% mutate(site = site)
  }
  df <- bind_rows(rows)
  if (nrow(df) == 0) return(ggplot() + theme_void())

  df <- df %>% filter(!is.na(F_FG), !is.na(F_EC), abs(F_FG) < 1e4, abs(F_EC) < 1e4)

  stats <- df %>%
    group_by(site) %>%
    summarise(
      n = n(),
      r = cor(F_FG, F_EC, use = "complete.obs"),
      slope = coef(lm(F_FG ~ F_EC))[2],
      .groups = "drop"
    ) %>%
    left_join(canopy_metrics %>% select(site, canopy_class, canopy_height),
              by = "site") %>%
    filter(!is.na(r), n > 10)

  stats$site <- reorder(stats$site, stats$r)

  ggplot(stats, aes(x = site, y = r, fill = canopy_class)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c(
      "No canopy (<1 m)"        = "#E8C547",
      "Short canopy (1-5 m)"    = "#A8D95B",
      "Medium canopy (5-15 m)"  = "#5BBA6F",
      "Tall canopy (15-30 m)"   = "#2D8B3A",
      "Very tall canopy (>30 m)"= "#0F5A1E"
    ), name = "Canopy class") +
    coord_flip() +
    labs(
      x = NULL, y = "Pearson r (F_FG vs F_EC)",
      title = paste0(gas, " FG-vs-EC Closure Correlation by Site")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
}
