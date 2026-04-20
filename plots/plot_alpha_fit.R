# plot_alpha_fit.R
# Diagnostic plot of the within-canopy wind-profile fit used to calibrate
# Raupach α per site. Lets you see the raw data and the quality of the fit.

#' Diagnostic plot for a single site's α fit.
#'
#' Two stacked panels:
#'   (top)    raw within-canopy wind data on the linear plane of the Raupach
#'            model — y = ln(u), x = (1 − z/h). Annual fit as dashed line.
#'   (bottom) monthly α (12 values, as points) overlaid with the SINUSOIDAL
#'            smooth α(DOY) (solid line, the production model), and the
#'            annual α (red dashed) + default α=3 (grey dotted) for reference.
#'
#' @param site        character
#' @param winds_site  tibble (timeMid, TowerPosition, height_m, ubar)
#' @param canopy_h    numeric
#' @param alpha_tbl   output of fit_alpha_by_site_month()  (monthly)
#' @param sin_tbl     output of fit_alpha_sinusoidal()   (optional, for smooth)
#' @return patchwork plot
plot_alpha_fit_diagnostic <- function(site, winds_site, canopy_h,
                                       alpha_tbl, sin_tbl = NULL,
                                       tier_tbl = NULL) {

  if (is.null(winds_site) || nrow(winds_site) == 0) {
    return(ggplot() + theme_void() +
             labs(title = paste(site, "— no wind data")))
  }
  if (is.na(canopy_h) || canopy_h < 1) {
    return(ggplot() + theme_void() +
             labs(title = paste(site, "— no canopy (open site)"),
                  subtitle = "K(z) uses MOST only; α not applicable"))
  }

  # Within-canopy only
  d <- winds_site %>%
    filter(height_m <= canopy_h, ubar > 0.01, !is.na(ubar)) %>%
    mutate(x = 1 - height_m / canopy_h,
           y = log(ubar),
           month = lubridate::month(timeMid))

  if (nrow(d) < 20) {
    return(ggplot() + theme_void() +
             labs(title = paste(site, "— insufficient within-canopy data")))
  }

  # Pull α fit results for this site
  ann_row <- alpha_tbl %>%
    filter(site == !!site, period_type == "annual")
  mon_rows <- alpha_tbl %>%
    filter(site == !!site, period_type == "month")

  alpha_ann <- ann_row$alpha
  r2_ann    <- ann_row$r2
  n_ann     <- ann_row$n_obs
  beta_ann  <- ann_row$beta

  # ---------- TOP PANEL: raw points + fit line ----------
  # Aggregate by height and month for plotting (median y per bin)
  pts <- d %>%
    group_by(month, height_m, x) %>%
    summarise(y = median(y, na.rm = TRUE), n = n(), .groups = "drop")

  # Intercept of annual fit: a = mean(y) + beta*mean(x)
  a_hat <- mean(d$y, na.rm = TRUE) + beta_ann * mean(d$x, na.rm = TRUE)
  xline <- seq(min(pts$x), max(pts$x), length.out = 100)
  yline <- a_hat - beta_ann * xline
  fit_df <- tibble(x = xline, y = yline)

  p_top <- ggplot(pts, aes(x = x, y = y)) +
    geom_point(aes(color = factor(month), size = n), alpha = 0.8) +
    geom_line(data = fit_df, aes(x = x, y = y),
              color = "black", linewidth = 0.8, linetype = "dashed") +
    scale_color_viridis_d(option = "D", name = "month") +
    scale_size_continuous(range = c(1.5, 4), guide = "none") +
    annotate("label", x = -Inf, y = Inf,
             label = sprintf("annual α = %.2f   β = %.2f   r² = %.2f   n = %d",
                              alpha_ann, beta_ann, r2_ann, n_ann),
             hjust = -0.05, vjust = 1.3,
             label.size = 0, fill = "white",
             size = 3.2, color = "grey20") +
    labs(
      x = "1 − z / canopy_h",
      y = "ln(u)",
      title = paste0(site, "  —  α fit diagnostic (canopy ",
                       round(canopy_h, 1), " m)"),
      subtitle = "points = monthly median wind at each sensor height; dashed = annual fit"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())

  # ---------- BOTTOM PANEL: monthly α ----------
  mon <- mon_rows %>%
    mutate(ok = note == "ok" & !is.na(alpha))

  # Which tier is this site using, and what constant / curve applies?
  tier_row <- NULL
  if (!is.null(tier_tbl)) {
    tr <- tier_tbl[tier_tbl$site == site, , drop = FALSE]
    if (nrow(tr) == 1) tier_row <- tr
  }
  tier_label <- if (!is.null(tier_row))
    sprintf("TIER %d — %s  (α used = %.2f)",
             tier_row$tier, tier_row$reason, tier_row$alpha_used)
  else
    "no tier assignment"

  # Build sinusoid curve if coefficients available for this site.
  # Only draw it as a line if this site actually gets tier 1 (uses the
  # sinusoid in production). If tier 2+, draw a flat constant α instead.
  sin_curve <- NULL
  sin_label <- NULL
  prod_line <- NULL
  if (!is.null(sin_tbl)) {
    srow <- sin_tbl[sin_tbl$site == site, , drop = FALSE]
    if (nrow(srow) == 1 && !is.na(srow$alpha_mean)) {
      # Sinusoid now directly on α: α(DOY) = A0 + Ac·cos + As·sin
      doy_grid <- seq(1, 365, by = 2)
      alpha_raw <- srow$alpha_mean +
                   (if (!is.na(srow$Ac)) srow$Ac * cos(2 * pi * doy_grid / 365) else 0) +
                   (if (!is.na(srow$As)) srow$As * sin(2 * pi * doy_grid / 365) else 0)
      # Per-site cap at observed monthly min/max (what production uses)
      lo <- if (!is.null(srow$alpha_obs_min) && !is.na(srow$alpha_obs_min))
              max(srow$alpha_obs_min, 0) else 0
      hi <- if (!is.null(srow$alpha_obs_max) && !is.na(srow$alpha_obs_max))
              srow$alpha_obs_max else Inf
      alpha_t <- pmin(pmax(alpha_raw, lo), hi)
      month_pos <- 1 + 11 * (doy_grid - 1) / 364
      sin_curve <- tibble(x = month_pos, alpha = alpha_t)
      p_fmt <- if (!is.na(srow$f_p)) sprintf("%.3g", srow$f_p) else "NA"
      sig_txt <- if (isTRUE(srow$seasonal_significant))
        "significant" else "not significant"
      sin_label <- sprintf(
        paste0("sinusoid (on monthly α): A0=%.2f  amp=%.2f  peak DOY %.0f  ",
                "F-test p=%s (%s); clamped to observed [%.2f, %.2f]"),
        srow$alpha_mean, srow$alpha_amp, srow$phase_doy, p_fmt, sig_txt,
        lo, if (is.finite(hi)) hi else Inf)
    }
  }
  if (!is.null(tier_row)) {
    # production curve: sinusoid only if tier 1, else flat
    if (tier_row$tier == 1L && !is.null(sin_curve)) {
      prod_line <- sin_curve %>% mutate(kind = "production")
    } else {
      prod_line <- tibble(x = c(0.5, 13.5),
                           alpha = rep(tier_row$alpha_used, 2),
                           kind = "production")
    }
  }

  p_bot <- ggplot(mon, aes(x = period_value, y = alpha)) +
    geom_hline(yintercept = 3,
               color = "grey55", linetype = "dotted", linewidth = 0.4) +
    geom_hline(yintercept = alpha_ann,
               color = "red", linetype = "dashed", linewidth = 0.5) +
    {if (!is.null(sin_curve))
      geom_line(data = sin_curve, aes(x = x, y = alpha),
                color = "grey55", linetype = "dashed", linewidth = 0.5,
                inherit.aes = FALSE)
     else NULL} +
    {if (!is.null(prod_line))
      geom_line(data = prod_line, aes(x = x, y = alpha),
                color = "#1F4E79", linewidth = 1.3, inherit.aes = FALSE)
     else NULL} +
    geom_point(aes(size = r2, color = ok), alpha = 0.85) +
    scale_size_continuous(range = c(1.5, 5), name = "r²") +
    scale_color_manual(values = c(`TRUE` = "steelblue",
                                    `FALSE` = "grey70"),
                        labels = c(`TRUE` = "monthly fit ok",
                                     `FALSE` = "monthly flagged"),
                        name = NULL) +
    annotate("text", x = 12.3, y = alpha_ann, hjust = 0,
             label = sprintf("annual\nα = %.2f", alpha_ann),
             size = 3, color = "red", lineheight = 0.9) +
    scale_x_continuous(breaks = 1:12,
                        labels = month.abb,
                        limits = c(0.5, 13.5)) +
    labs(x = NULL, y = "α",
         title = NULL,
         subtitle = paste0(
           tier_label,
           "\nBlue solid = α used in production (sinusoid if tier 1, else constant).",
           " Grey dashed = sinusoid fit for reference.",
           if (!is.null(sin_label)) paste0("\n", sin_label) else "")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right")

  p_top / p_bot + patchwork::plot_layout(heights = c(2, 1.2))
}
