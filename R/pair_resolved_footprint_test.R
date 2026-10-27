# pair_resolved_footprint_test.R
# The headline mechanism-2 test for the manuscript.
#
# For each adjacent sensor pair (z_lower, z_upper) at each 30-min:
#   dC(CH4)  ~  ΔA + Δf_wet
# Where:
#   ΔA      = area(footprint_upper) - area(footprint_lower)   (m^2)
#   Δf_wet  = f_wet(footprint_upper) - f_wet(footprint_lower) (unitless)
#   dC      = mean(upper) - mean(lower)                       (gas conc units)
#
# Mechanism-2 (hidden wet microsites) predicts dC scales with Δf_wet (and
# secondarily with ΔA since larger footprints are more likely to intercept
# wet patches). Mechanism-1 (canopy emissions) predicts no relationship at
# upland-only or low-statured sites, with the canopy-interface signal showing
# up as site-level residual structure.
#
# Levels of analysis:
#   per-site regression       -> site-level slope estimates
#   per-canopy-class           -> stratified pooled regressions
#   continental pooled         -> grand-mean slope (mixed model with site RE)
#
# Residual classification:
#   sites whose mean dC at fixed Δf_wet bins is *above* the continental
#   regression line have unexplained CH4 enhancement beyond what footprint
#   wetness explains -> mechanism-1 candidate
#   sites with strong positive Δf_wet slopes -> mechanism-2 dominated
#   sites with neither -> neutral / closed canopy / low-flux

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(broom)
  if (requireNamespace("lme4", quietly = TRUE)) library(lme4)
})


#' Restrict to adjacent pairs only (level difference == 1) for the headline test.
#' Optionally also accept bottom-top pairs (long differences) as a sensitivity check.
#'
#' @param pair_df Tibble. Output of compute_footprint_pair_metrics().
#' @param adjacent_only Logical. Default TRUE.
filter_pair_metrics <- function(pair_df, adjacent_only = TRUE,
                                  min_dArea = 1, drop_na_dfwet = TRUE,
                                  jja_night_only = FALSE) {

  out <- pair_df %>%
    filter(!is.na(dConc), !is.na(dArea), is.finite(dConc), is.finite(dArea))

  if (drop_na_dfwet) {
    out <- out %>% filter(!is.na(dfwet))
  }

  if (adjacent_only) {
    out <- out %>%
      mutate(level_a = as.numeric(sub("_.*", "", dLevelsAminusB)),
             level_b = as.numeric(sub(".*_", "", dLevelsAminusB))) %>%
      filter(abs(level_a - level_b) == 1) %>%
      select(-level_a, -level_b)
  }

  # Require non-zero footprint expansion (some pairs have ΔA ≈ 0)
  out <- out %>% filter(abs(dArea) > min_dArea)

  if (jja_night_only) {
    out <- out %>%
      mutate(month = lubridate::month(timeMid),
             hour  = lubridate::hour(timeMid)) %>%
      filter(month %in% c(6, 7, 8), hour >= 22 | hour < 5)
  }

  out
}


#' Per-site regressions: dConc ~ dfwet and dConc ~ dArea (separately).
#'
#' Returns a tidy per-site tibble of slopes, SE, p, R^2, n.
#'
#' @param pair_df Tibble. Filtered pair metrics for ONE gas across sites
#'   (must contain a `site` column, added by the caller).
#' @return Tibble (one row per site) with slope, SE, R^2, n.
fit_per_site_regression <- function(pair_df, gas = "CH4") {

  fit_one <- function(df) {
    safe_coef <- function(formula_str) {
      # Return all-NA (typed as numeric) if predictor has zero variance,
      # if lm fails, or if the predictor isn't estimable.
      na_out <- list(slope = NA_real_, se = NA_real_,
                     p = NA_real_, r2 = NA_real_)
      pred <- strsplit(formula_str, " ~ ")[[1]][2]
      if (!pred %in% names(df) || sum(!is.na(df[[pred]])) < 3 ||
          length(unique(df[[pred]][!is.na(df[[pred]])])) < 2) {
        return(na_out)
      }
      f <- tryCatch(lm(as.formula(formula_str), data = df),
                    error = function(e) NULL)
      if (is.null(f) || !(pred %in% rownames(summary(f)$coefficients))) {
        return(na_out)
      }
      list(
        slope = as.numeric(summary(f)$coefficients[pred, "Estimate"]),
        se    = as.numeric(summary(f)$coefficients[pred, "Std. Error"]),
        p     = as.numeric(summary(f)$coefficients[pred, "Pr(>|t|)"]),
        r2    = as.numeric(summary(f)$r.squared)
      )
    }
    fw <- safe_coef("dConc ~ dfwet")
    ar <- safe_coef("dConc ~ dArea")
    tibble(
      n = nrow(df),
      slope_dfwet = fw$slope, se_dfwet = fw$se,
      p_dfwet = fw$p, r2_dfwet = fw$r2,
      slope_dArea = ar$slope, se_dArea = ar$se,
      p_dArea = ar$p, r2_dArea = ar$r2,
      gas = gas
    )
  }

  pair_df %>%
    split(.$site) %>%
    purrr::imap_dfr(function(df, s) {
      cbind(site = s, fit_one(df))
    }) %>%
    as_tibble()
}


#' Continental pooled regression with site as a random effect.
#' Falls back to fixed-effect site if lme4 not available.
#'
#' @return Single-row tibble with grand slope, SE, t-value, n_total, n_sites.
fit_continental_regression <- function(pair_df, gas = "CH4") {

  if (requireNamespace("lme4", quietly = TRUE)) {
    fit <- tryCatch(
      lme4::lmer(dConc ~ dfwet + dArea + (1 + dfwet | site),
                  data = pair_df, REML = TRUE),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      fe <- summary(fit)$coefficients
      return(tibble(
        gas = gas,
        slope_dfwet = fe["dfwet", "Estimate"],
        se_dfwet    = fe["dfwet", "Std. Error"],
        t_dfwet     = fe["dfwet", "t value"],
        slope_dArea = fe["dArea", "Estimate"],
        se_dArea    = fe["dArea", "Std. Error"],
        t_dArea     = fe["dArea", "t value"],
        n_total     = nrow(pair_df),
        n_sites     = length(unique(pair_df$site)),
        method      = "lmer(random slope by site)"
      ))
    }
  }

  # Fallback: pooled OLS with site fixed effect
  fit <- lm(dConc ~ dfwet + dArea + factor(site), data = pair_df)
  coefs <- summary(fit)$coefficients
  tibble(
    gas = gas,
    slope_dfwet = coefs["dfwet", "Estimate"],
    se_dfwet    = coefs["dfwet", "Std. Error"],
    t_dfwet     = coefs["dfwet", "t value"],
    slope_dArea = coefs["dArea", "Estimate"],
    se_dArea    = coefs["dArea", "Std. Error"],
    t_dArea     = coefs["dArea", "t value"],
    n_total = nrow(pair_df),
    n_sites = length(unique(pair_df$site)),
    method = "lm(fixed effect site)"
  )
}


#' Classify each site as mechanism-1-dominant, mechanism-2-dominant, or mixed,
#' based on its dfwet slope and unexplained intercept relative to the
#' continental regression.
#'
#' Rule of thumb (refine later with data):
#'   slope_dfwet positive AND p < 0.05 AND r2 > 0.05      -> mechanism 2 active
#'   intercept_residual positive AND > 1 SE above 0       -> mechanism 1 candidate
#'   neither                                              -> neutral / closed canopy
#'
#' Sites with no NWI wetlands within the buffer (Δf_wet ≡ 0) get classified by
#' ΔArea slope alone — mechanism 2 cannot be tested at these.
classify_mechanism <- function(per_site_df, canopy_metrics = NULL,
                                 slope_threshold_p = 0.05, r2_threshold = 0.05) {
  per_site_df %>%
    mutate(
      mech_class = case_when(
        is.na(slope_dfwet) ~ "untested",
        slope_dfwet > 0 & p_dfwet < slope_threshold_p & r2_dfwet > r2_threshold ~ "mechanism 2 active",
        slope_dfwet < 0 & p_dfwet < slope_threshold_p ~ "negative dfwet (anomalous)",
        TRUE ~ "mechanism 1 / closed / no wetness"
      )
    ) %>%
    {
      if (!is.null(canopy_metrics)) {
        left_join(., canopy_metrics %>% select(site, canopy_class, canopy_height),
                  by = "site")
      } else .
    }
}
