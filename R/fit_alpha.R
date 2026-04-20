# fit_alpha.R
# Site-specific (and optionally month-specific) calibration of the Raupach
# canopy attenuation coefficient α from 2D wind-speed profiles.
#
# Theory
# ------
# Within-canopy wind speed under Raupach (1989):
#     u(z) = u(h) * exp(-β * (1 - z/h))                  (z <= h)
# Equivalently, for any within-canopy observations i:
#     ln(u_i) = a - β * (1 - z_i/h)       with a = ln(u(h))
#
# Raupach's same framework for scalar eddy diffusivity gives:
#     K(z) = K(h) * exp(-α * (1 - z/h))
# and mixing-length theory (used consistently in Massman 1997, Harman &
# Finnigan 2007) gives α ≈ 2·β when the in-canopy mixing length decays at the
# same rate as the wind — the standard closure assumption for forest canopies.
#
# We fit β by pooled linear regression within a site × period (month, season,
# or year). α is then 2·β. Diagnostics returned per period:
#   n_obs, n_heights, r2, beta, alpha, note.


#' Fit a single β (and α) from pooled within-canopy wind data.
#'
#' @param df Data frame with columns height_m, ubar, timeMid.
#' @param canopy_h Numeric. Canopy height (m).
#' @param min_obs Minimum total within-canopy observations required.
#' @return Named list: alpha, beta, r2, n_obs, n_heights, note.
.fit_beta_from_pool <- function(df, canopy_h,
                                  min_obs = 50, min_heights = 2) {

  out <- list(alpha = NA_real_, beta = NA_real_, r2 = NA_real_,
              n_obs = nrow(df), n_heights = 0, note = NA_character_)

  # Within canopy only, valid ubar
  d <- df[df$height_m <= canopy_h & df$ubar > 0.01 & !is.na(df$ubar), ,
          drop = FALSE]
  out$n_obs     <- nrow(d)
  out$n_heights <- length(unique(d$height_m))

  if (out$n_obs < min_obs) { out$note <- "insufficient obs"; return(out) }
  if (out$n_heights < min_heights) {
    out$note <- "too few within-canopy heights"; return(out)
  }

  x <- 1 - d$height_m / canopy_h     # so the slope is −β
  y <- log(d$ubar)
  fit <- tryCatch(lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) { out$note <- "lm failed"; return(out) }

  cf <- coef(fit)
  beta <- unname(-cf[2])
  r2   <- summary(fit)$r.squared

  out$beta  <- beta
  out$alpha <- 2 * beta
  out$r2    <- r2
  out$note  <- if (is.na(beta) || beta < 0)
    "negative beta" else if (r2 < 0.1)
    "poor fit"      else "ok"
  out
}


#' Fit α per site per month (and annual), returning a long-format table.
#'
#' Grassland / shrubland / desert / tundra sites with canopy_h < 1 m are
#' skipped — no canopy to attenuate; our K(z) code uses MOST-only for those.
#'
#' @param all_winds Named list (by site) of wind tibbles (extract_wind_profile output).
#' @param canopy_metrics Tibble with site, canopy_height.
#' @param default_alpha Fallback α when the fit is unreliable. Default 3.
#' @param min_obs Minimum within-canopy obs per period to fit.
#' @param min_heights Minimum distinct within-canopy heights required.
#' @return Tibble: site, period_type ("annual" or "month"), period_value,
#'   alpha, beta, r2, n_obs, n_heights, note, canopy_h.
fit_alpha_by_site_month <- function(all_winds, canopy_metrics,
                                      default_alpha = 3,
                                      min_obs = 200,
                                      min_heights = 2) {

  rows <- list()
  for (site in names(all_winds)) {
    wind <- all_winds[[site]]
    if (is.null(wind) || nrow(wind) == 0) next
    ch <- canopy_metrics$canopy_height[canopy_metrics$site == site]
    if (length(ch) == 0 || is.na(ch) || ch < 1) next   # open sites: skip

    wind <- wind %>% mutate(mo = lubridate::month(timeMid))

    # Annual fit (all months pooled)
    ann <- .fit_beta_from_pool(wind, ch, min_obs = min_obs,
                                 min_heights = min_heights)
    rows[[paste(site, "ann", sep = "_")]] <- tibble(
      site = site, period_type = "annual", period_value = NA_integer_,
      alpha = ann$alpha, beta = ann$beta, r2 = ann$r2,
      n_obs = ann$n_obs, n_heights = ann$n_heights,
      note = ann$note, canopy_h = ch)

    # Monthly fits
    for (m in 1:12) {
      dm <- wind[wind$mo == m, , drop = FALSE]
      fm <- .fit_beta_from_pool(dm, ch, min_obs = min_obs,
                                  min_heights = min_heights)
      rows[[paste(site, "m", m, sep = "_")]] <- tibble(
        site = site, period_type = "month", period_value = m,
        alpha = fm$alpha, beta = fm$beta, r2 = fm$r2,
        n_obs = fm$n_obs, n_heights = fm$n_heights,
        note = fm$note, canopy_h = ch)
    }
  }

  bind_rows(rows) %>%
    mutate(alpha_use = pmin(pmax(alpha, 0), 8),    # physical range [0, 8]
            alpha_use = ifelse(is.na(alpha_use) | note %in%
                                  c("insufficient obs",
                                    "too few within-canopy heights",
                                    "lm failed", "negative beta"),
                                NA_real_, alpha_use))
}


#' Fit a sinusoidal seasonal β(DOY) model for one site.
#'
#' β(DOY) = β_mean + A_c·cos(2πDOY/365) + A_s·sin(2πDOY/365)
#'
#' ln(u_i) = a − β(DOY_i) · x_i,  where  x_i = 1 − z_i / h
#'         = a − β_mean·x_i − A_c·(cos·x)_i − A_s·(sin·x)_i
#'
#' So `lm(y ~ x + cos_x + sin_x)` gives all parameters. α = 2·β.
#'
#' @param df        winds tibble with height_m, ubar, timeMid
#' @param canopy_h  canopy height (m)
#' @param min_obs   minimum within-canopy obs
#' @param min_heights minimum distinct within-canopy heights
#' @return list: beta_mean, Ac, As, amplitude, phase_doy, alpha_mean,
#'         alpha_amp, r2, n_obs, n_heights, note
.fit_sinusoidal_from_pool <- function(df, canopy_h,
                                        min_obs = 200, min_heights = 2) {

  out <- list(beta_mean = NA_real_, Ac = NA_real_, As = NA_real_,
              amplitude = NA_real_, phase_doy = NA_real_,
              alpha_mean = NA_real_, alpha_amp = NA_real_,
              r2 = NA_real_, n_obs = nrow(df),
              n_heights = 0, note = NA_character_)

  d <- df[df$height_m <= canopy_h & df$ubar > 0.01 & !is.na(df$ubar), ,
          drop = FALSE]
  out$n_obs     <- nrow(d)
  out$n_heights <- length(unique(d$height_m))
  if (out$n_obs < min_obs) { out$note <- "insufficient obs"; return(out) }
  if (out$n_heights < min_heights) {
    out$note <- "too few within-canopy heights"; return(out)
  }

  doy <- lubridate::yday(d$timeMid)
  omega <- 2 * pi * doy / 365
  x     <- 1 - d$height_m / canopy_h
  y     <- log(d$ubar)

  cos_x <- cos(omega) * x
  sin_x <- sin(omega) * x

  fit <- tryCatch(lm(y ~ x + cos_x + sin_x), error = function(e) NULL)
  if (is.null(fit)) { out$note <- "lm failed"; return(out) }

  cf <- coef(fit)
  # Sign: y = a − β·x  ⇒  coefficient on x is −β_mean, on cos_x is −A_c, etc.
  beta_mean <- unname(-cf["x"])
  Ac        <- unname(-cf["cos_x"])
  As        <- unname(-cf["sin_x"])
  amp       <- sqrt(Ac^2 + As^2)
  phase_doy <- (atan2(As, Ac) / (2 * pi) * 365) %% 365
  r2        <- summary(fit)$r.squared

  out$beta_mean  <- beta_mean
  out$Ac         <- Ac
  out$As         <- As
  out$amplitude  <- amp
  out$phase_doy  <- phase_doy
  out$alpha_mean <- 2 * beta_mean
  out$alpha_amp  <- 2 * amp
  out$r2         <- r2
  out$note       <- if (is.na(beta_mean) || beta_mean < 0)
    "negative beta_mean" else if (r2 < 0.1)
    "poor fit"           else "ok"
  out
}

#' Fit a sinusoid y = A0 + Ac·cos(omega·m) + As·sin(omega·m) to the 12 monthly
#' α values, weighted by the per-month fit r². Also test whether it's a
#' significantly better fit than constant via F-test.
#'
#' This avoids the pooled-pooling problem: the pooled-raw sinusoid conflates
#' seasonal intercept shifts (synoptic wind magnitude) with seasonal slope
#' shifts (canopy attenuation), producing overshoot and phase errors.
#' Fitting on monthly α values uses the already-estimated per-month slope,
#' so only genuine α-vs-DOY variation remains to explain.
#'
#' Also returns the observed min/max of the monthly α values so the sinusoid
#' output can be capped per-site at whatever range the actual data supports
#' (no global priors).
#'
#' @param monthly_rows tibble from alpha_tbl filtered to one site, period_type=="month"
#' @return list(A0, Ac, As, amp, phase_doy, r2_sinusoid, r2_constant, f_p,
#'              seasonal_significant, note, alpha_obs_min, alpha_obs_max)
.fit_sinusoid_on_monthly <- function(monthly_rows, min_good_months = 4,
                                        p_sig = 0.05) {
  out <- list(A0 = NA_real_, Ac = NA_real_, As = NA_real_,
              amp = NA_real_, phase_doy = NA_real_,
              r2_sinusoid = NA_real_, r2_constant = NA_real_,
              f_p = NA_real_, seasonal_significant = FALSE,
              note = NA_character_, n_good_months = 0,
              alpha_obs_min = NA_real_, alpha_obs_max = NA_real_)

  good <- monthly_rows %>%
    filter(note == "ok", !is.na(alpha_use))
  out$n_good_months <- nrow(good)
  if (nrow(good) > 0) {
    out$alpha_obs_min <- min(good$alpha_use, na.rm = TRUE)
    out$alpha_obs_max <- max(good$alpha_use, na.rm = TRUE)
  }

  if (nrow(good) < min_good_months) {
    out$note <- "too few good monthly fits"
    if (nrow(good) > 0) out$A0 <- mean(good$alpha_use)
    return(out)
  }

  # Mid-month DOY
  mid_doy <- c(15, 46, 75, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  good <- good %>%
    mutate(doy = mid_doy[period_value],
            om = 2 * pi * doy / 365,
            cos_om = cos(om),
            sin_om = sin(om))

  w <- good$r2                      # weight monthly α values by fit r²
  w[!is.finite(w) | w < 0] <- 0

  if (sum(w) == 0) w <- rep(1, nrow(good))  # fallback equal weights

  # Constant model
  fit0 <- tryCatch(lm(alpha_use ~ 1, data = good, weights = w),
                    error = function(e) NULL)
  # Sinusoid model
  fit1 <- tryCatch(lm(alpha_use ~ cos_om + sin_om, data = good, weights = w),
                    error = function(e) NULL)
  if (is.null(fit0) || is.null(fit1)) { out$note <- "lm failed"; return(out) }

  r2_0 <- summary(fit0)$r.squared
  r2_1 <- summary(fit1)$r.squared

  # F-test comparing constant vs sinusoid
  fa <- tryCatch(anova(fit0, fit1), error = function(e) NULL)
  p_f <- if (!is.null(fa)) fa$`Pr(>F)`[2] else NA_real_

  cf <- coef(fit1)
  A0 <- unname(cf["(Intercept)"])
  Ac <- unname(cf["cos_om"])
  As <- unname(cf["sin_om"])
  amp <- sqrt(Ac^2 + As^2)
  phase_doy <- (atan2(As, Ac) / (2 * pi) * 365) %% 365

  out$A0 <- A0; out$Ac <- Ac; out$As <- As
  out$amp <- amp; out$phase_doy <- phase_doy
  out$r2_sinusoid <- r2_1
  out$r2_constant <- r2_0
  out$f_p <- p_f
  out$seasonal_significant <- !is.na(p_f) && p_f < p_sig
  out$note <- if (out$seasonal_significant)
    "sinusoid significant" else "sinusoid not significant"
  out
}

#' Site-level α(DOY) fits:
#'   - monthly α values from pooled regression per month
#'   - sinusoid fit on those monthly values (not on raw obs)
#'   - F-test: is sinusoid significantly better than constant?
#'
#' @return tibble with columns:
#'   site, canopy_h, alpha_mean (= sinusoid A0 when significant, else mean of
#'   monthly values), alpha_amp (= sinusoid amplitude), phase_doy,
#'   Ac, As, r2_sinusoid, r2_constant, f_p, seasonal_significant, n_good_months.
fit_alpha_sinusoidal <- function(all_winds, canopy_metrics,
                                   monthly_tbl,
                                   min_good_months = 6,
                                   p_sig = 0.05, ...) {
  rows <- list()
  for (site in names(all_winds)) {
    ch <- canopy_metrics$canopy_height[canopy_metrics$site == site]
    if (length(ch) == 0 || is.na(ch) || ch < 1) next

    m_rows <- monthly_tbl %>%
      filter(site == !!site, period_type == "month")
    if (nrow(m_rows) == 0) next

    f <- .fit_sinusoid_on_monthly(m_rows, min_good_months = min_good_months,
                                     p_sig = p_sig)
    mean_mo <- m_rows %>% filter(note == "ok") %>% pull(alpha_use) %>%
      mean(na.rm = TRUE)
    rows[[site]] <- tibble(
      site = site, canopy_h = ch,
      alpha_mean = if (!is.na(f$A0)) f$A0 else mean_mo,
      alpha_amp  = f$amp,
      phase_doy  = f$phase_doy,
      Ac = f$Ac, As = f$As,
      alpha_obs_min = f$alpha_obs_min,
      alpha_obs_max = f$alpha_obs_max,
      r2_sinusoid = f$r2_sinusoid,
      r2_constant = f$r2_constant,
      f_p = f$f_p,
      seasonal_significant = f$seasonal_significant,
      n_good_months = f$n_good_months,
      note = f$note,
      # legacy columns for back-compat
      beta_mean = if (!is.na(f$A0)) f$A0 / 2 else NA_real_,
      r2 = f$r2_sinusoid,
      n_obs = sum(m_rows$n_obs, na.rm = TRUE),
      n_heights = max(m_rows$n_heights, na.rm = TRUE))
  }
  bind_rows(rows)
}

#' Decide which fit tier each site passes.
#'
#' Gating uses FIT-QUALITY criteria only — no physical priors on phase or
#' amplitude direction. α reports what the air actually does, whether the
#' canopy is leafy, bare, or snow-loaded; whether the attenuation is
#' zero (open site) or strong (dense canopy); whether peak drag is in
#' summer (foliage) or winter (snow). Trust the wind profile.
#'
#' Tiers:
#'   1 sinusoid(DOY): r² ≥ r2_min_sin AND β_mean ≥ 0 AND α_mean ≤ alpha_max
#'   2 constant α_mean from sinusoid fit:
#'       β_mean ≥ 0 AND α_mean ≤ alpha_max (r² below sinusoid threshold,
#'       so seasonal term isn't significant)
#'   25 constant from annual pooled regression:
#'       sinusoid itself unusable (β_mean < 0 or α_mean > alpha_max), but
#'       annual-pooled regression produces a reasonable constant
#'   3 canopy-class median (all site fits bad)
#'   4 global default
#'
#' Note: α_mean = 0 and α_mean close to 0 are PHYSICAL (sparse / leaf-off
#' canopies). α_max is a numerical clamp for very dense canopies.
classify_alpha_tiers <- function(sin_tbl, canopy_metrics,
                                  default_alpha = 3,
                                  p_sig = 0.05,
                                  annual_tbl = NULL) {

  d <- sin_tbl %>%
    left_join(canopy_metrics %>% select(site, canopy_class), by = "site") %>%
    mutate(
      pass_seasonal = !is.na(seasonal_significant) & seasonal_significant &
                        !is.na(f_p) & f_p < p_sig,
      has_mean      = !is.na(alpha_mean) & is.finite(alpha_mean),

      # Tier 1: sinusoid significantly better than constant.
      # Its output is clamped per-site at [alpha_obs_min, alpha_obs_max]
      # so extrapolation never exceeds the observed monthly range.
      tier1_pass = pass_seasonal & has_mean,
      # Tier 2: use site annual α_mean (from sinusoid fit on monthlies; when
      # sinusoid is not significant, A0 = mean of monthly α values).
      tier2_pass = has_mean
    ) %>%
    mutate(
      tier = case_when(
        tier1_pass            ~ 1L,
        tier2_pass            ~ 2L,
        TRUE                  ~ NA_integer_
      )
    )

  # Tier 3 fallback: median α_mean within each canopy class (using tier 1 & 2)
  class_med <- d %>%
    filter(!is.na(tier)) %>%
    group_by(canopy_class) %>%
    summarise(class_alpha = median(alpha_mean, na.rm = TRUE),
              n_in_class = n(), .groups = "drop")

  # Tier 2.5: independent annual pooled regression (no cos/sin term).
  # Use when the sinusoid has nonphysical mean (β < 0 or α > cap).
  ann_fallback <- if (!is.null(annual_tbl)) {
    annual_tbl %>% filter(period_type == "annual") %>%
      transmute(site, alpha_annual_reg = alpha_use)
  } else tibble(site = character(), alpha_annual_reg = numeric())

  d <- d %>%
    left_join(class_med, by = "canopy_class") %>%
    left_join(ann_fallback, by = "site") %>%
    mutate(
      tier = case_when(
        !is.na(tier) ~ tier,
        !is.na(alpha_annual_reg) & is.finite(alpha_annual_reg) &
          alpha_annual_reg >= 0 ~ 25L,
        !is.na(class_alpha) ~ 3L,
        TRUE ~ 4L
      ),
      alpha_used = case_when(
        tier == 1L  ~ alpha_mean,              # sinusoid (temporally varying)
        tier == 2L  ~ alpha_mean,              # constant from sinusoid fit
        tier == 25L ~ alpha_annual_reg,        # constant from annual regression
        tier == 3L  ~ class_alpha,             # class median
        tier == 4L  ~ default_alpha,           # global fallback
        TRUE        ~ NA_real_
      ),
      reason = case_when(
        tier == 1L  ~ "sinusoid significantly better than constant (F-test); capped at observed range",
        tier == 2L  ~ "sinusoid not significant → site α_mean (annual)",
        tier == 25L ~ "site α_mean unavailable → annual regression",
        tier == 3L  ~ "site fits unusable → class median",
        tier == 4L  ~ "no class median → default"
      )
    )
  d
}


#' Lookup function for α at any (site, date), with 4-tier gated cascade.
#'
#' @param sin_tbl        output of fit_alpha_sinusoidal()
#' @param canopy_metrics tibble with canopy_class per site
#' @param default_alpha  fallback
#' @return list(fn, tier_tbl) — fn is function(site, date_or_doy),
#'         tier_tbl records which tier each site used.
make_alpha_lookup_sinusoidal <- function(sin_tbl, canopy_metrics,
                                           default_alpha = 3,
                                           annual_tbl = NULL, ...) {

  tier_tbl <- classify_alpha_tiers(sin_tbl, canopy_metrics,
                                     default_alpha = default_alpha,
                                     annual_tbl = annual_tbl, ...)

  fn <- function(site, date_or_doy) {
    doy <- if (inherits(date_or_doy, c("Date", "POSIXt")))
      as.numeric(lubridate::yday(date_or_doy)) else as.numeric(date_or_doy)

    row <- tier_tbl[tier_tbl$site == site, , drop = FALSE]
    if (nrow(row) == 0) return(default_alpha)

    if (isTRUE(row$tier == 1L)) {
      # Sinusoid fit is on α directly (from monthly α values).
      # α(DOY) = A0 + Ac·cos(ω) + As·sin(ω),  A0 = row$alpha_mean.
      # Cap output at the site's OBSERVED monthly min/max — no global priors,
      # just don't let extrapolated months exceed what was actually measured.
      alpha_t <- row$alpha_mean +
                 row$Ac * cos(2 * pi * doy / 365) +
                 row$As * sin(2 * pi * doy / 365)
      lo <- if (!is.na(row$alpha_obs_min)) row$alpha_obs_min else 0
      hi <- if (!is.na(row$alpha_obs_max)) row$alpha_obs_max else Inf
      # keep a hard non-negativity floor (wind can't increase into the canopy)
      lo <- pmax(lo, 0)
      return(pmin(pmax(alpha_t, lo), hi))
    }
    # tiers 2/3/4: constant
    if (is.na(row$alpha_used)) return(default_alpha)
    return(row$alpha_used)
  }

  list(fn = fn, tier_tbl = tier_tbl)
}


#' Build an (site, month) -> α lookup function with graceful fallbacks.
#'
#' (Legacy — kept for comparison with sinusoidal version.)
#' Priority:
#'   1. Monthly α for that (site, month)
#'   2. Annual α for that site
#'   3. Global default α
make_alpha_lookup <- function(alpha_tbl, default_alpha = 3) {

  ann <- alpha_tbl %>% filter(period_type == "annual") %>%
    select(site, alpha_ann = alpha_use)

  mon <- alpha_tbl %>% filter(period_type == "month") %>%
    select(site, mo = period_value, alpha_mo = alpha_use)

  function(site, month) {
    am <- mon$alpha_mo[mon$site == site & mon$mo == month]
    if (length(am) && !is.na(am)) return(am)
    aa <- ann$alpha_ann[ann$site == site]
    if (length(aa) && !is.na(aa)) return(aa)
    default_alpha
  }
}
