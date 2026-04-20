# compute_fluxes.R
# Timestamp-level (30-min) flux computation using canopy K(z) model
#
# For each 30-min window: join raw per-level concentrations with met (u*, L,
# wind profile), compute K(z) at each measurement height, then compute fluxes
# at each adjacent level pair. Also attach EC reference fluxes where available.
#
# ---------------- QC / filter scope (IMPORTANT) ----------------
# The filters applied in this file (ustar threshold, gradient SNR, Bad_Eddy
# outlier flag, Stability_500 flag) are FLUX-SPECIFIC. They are needed because
# flux = K * dC/dz depends on (a) sufficient turbulence for K(z) to be valid
# and (b) a detectable gradient above measurement noise.
#
# Concentration-only analyses — anomaly plots, profile heatmaps, storage flux,
# tracer-pair comparisons — use `prepare_profile_timeseries()` output directly.
# That function applies only NEON's basic sensor quality flag (qfFinl == 0)
# and keeps all valid concentration observations regardless of u* / SNR /
# stability. Do NOT gate those analyses behind the flux filters.
#
# Filter actions (matching parent lterwg-flux-gradient workflow):
#   qfFinl == 0     -> REMOVE (applied in prepare_profile_timeseries, upstream)
#   ustar >= 0.1    -> REMOVE
#   dC_SNR >= 3     -> REMOVE
#   Bad_Eddy (IQR outlier on F_FG)  -> FLAG ONLY (column F_outlier)
#   Stability_500 (stable/unstable/neutral based on L) -> FLAG ONLY

#' Aggregate raw 9-min profile data to 30-min per-level means AND propagate
#' measurement uncertainty so we can compute a pair SNR downstream.
#'
#' Each 9-min record has a within-9-min variance (`vari`) and sample count
#' (`numSamp`). SE of that 9-min mean = sqrt(vari / numSamp). The SE of the
#' 30-min mean of N 9-min means is sqrt(sum(vari_i / numSamp_i)) / N.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @return Tibble: time_round, TowerPosition, mean_conc, se_conc, n_obs
aggregate_profile_30min <- function(profile_ts) {
  if (nrow(profile_ts) == 0) return(tibble())
  profile_ts %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    group_by(time_round, TowerPosition) %>%
    summarise(
      mean_conc = mean(concentration, na.rm = TRUE),
      # Variance propagation: var of 30-min mean = (1/N^2) * sum(vari_i / numSamp_i)
      se_conc   = {
        contrib <- vari / numSamp
        contrib <- contrib[is.finite(contrib)]
        if (length(contrib) == 0) NA_real_ else sqrt(sum(contrib)) / length(contrib)
      },
      n_obs = sum(!is.na(concentration)),
      .groups = "drop"
    ) %>%
    filter(!is.na(mean_conc))
}


#' Extract 30-min EC reference flux for CO2 or H2O from aligned data
#'
#' @param aligned_data Data frame. One gas from min9Diff.list
#' @param gas Character. Gas name
#' @return Tibble: time_round, F_EC
extract_EC_flux <- function(aligned_data, gas) {

  if (is.null(aligned_data) || nrow(aligned_data) == 0) return(tibble())

  flux_col <- switch(gas,
    "CO2" = "FC_turb_interp",
    "H2O" = "LE_turb_interp",   # W/m² — will convert below
    NULL
  )
  if (is.null(flux_col) || !flux_col %in% names(aligned_data)) return(tibble())

  aligned_data %>%
    mutate(timeMid = if ("timeMid" %in% names(.)) {
      ymd_hms(timeMid, quiet = TRUE)
    } else if ("timeEnd_A" %in% names(.)) {
      ymd_hms(timeEnd_A, quiet = TRUE)
    } else NA_real_) %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    group_by(time_round) %>%
    summarise(F_EC = median(.data[[flux_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(time_round), !is.na(F_EC)) %>%
    {
      if (gas == "H2O") {
        # Convert LE (W/m²) to H2O flux (mmol/m²/s)
        # F_H2O = LE / (lambda * M_H2O) with lambda = 2.45e6 J/kg, M_H2O = 0.018 kg/mol
        # Result in mol/m²/s; convert to mmol/m²/s (×1000)
        mutate(., F_EC = F_EC / (2.45e6 * 0.018) * 1000)
      } else {
        .
      }
    }
}


#' Compute timestamp-level layer fluxes using K(z) model
#'
#' Joins 30-min per-level concentrations with met data, applies canopy K(z),
#' computes flux at each adjacent pair midpoint.
#'
#' @param profile_ts Tibble. Output of prepare_profile_timeseries()
#' @param met_ts Tibble. Output of extract_met_timeseries() — needs ustar, zoL
#' @param attr_df Data frame. Site attributes (canopy height, level heights)
#' @param aligned_data Data frame. Aligned data for EC reference flux (CO2/H2O)
#' @param gas Character. Gas name
#' @param alpha Numeric. Raupach decay parameter. Default 3.
#' @return Tibble: time_round, pair_label, z_lower, z_upper, z_mid, K, dCdz, F_FG,
#'   F_EC (where available), u_star, zoL, canopy_h
compute_layer_fluxes <- function(profile_ts, met_ts, attr_df,
                                  aligned_data = NULL, gas = "CH4",
                                  alpha = 3,
                                  ustar_threshold = 0.1,
                                  snr_threshold = 3,
                                  alpha_fn = NULL, site = NULL) {
  # alpha_fn: optional function(site, month) -> alpha. If supplied, used
  # per-timestamp in place of the scalar `alpha`. Requires `site` arg.

  if (nrow(profile_ts) == 0 || is.null(met_ts) || nrow(met_ts) == 0) {
    return(tibble())
  }

  height_lookup <- get_height_lookup(attr_df)
  canopy_h <- get_canopy_height(attr_df)

  # Aggregate concentrations to 30-min per level
  conc_30 <- aggregate_profile_30min(profile_ts)

  # Pivot wide: one row per time_round, columns per TowerPosition
  conc_wide <- conc_30 %>%
    left_join(height_lookup, by = "TowerPosition") %>%
    select(time_round, TowerPosition, height_m, mean_conc, se_conc)

  # Join with met
  met_30 <- met_ts %>%
    rename(time_round = timeMid) %>%
    select(time_round, ustar, zoL, any_of(c("L_obukhov", "rhoa_kgm3")))

  # Merge in EC reference flux if available
  ec <- if (!is.null(aligned_data)) extract_EC_flux(aligned_data, gas) else tibble()

  # Compute flux AT each tower level (not at pair midpoints).
  # dC/dz uses central differences where possible, one-sided at ends.
  # K is evaluated at the level height z_i directly.
  positions <- sort(unique(height_lookup$TowerPosition))
  heights   <- sapply(positions, function(p)
    height_lookup$height_m[height_lookup$TowerPosition == p])
  n_lev <- length(positions)
  if (n_lev < 2) return(tibble())

  # Pivot wide: rows = time_round, cols = C at each level (ordered low->high).
  # Also carry SE (se_conc) in a parallel wide matrix so we can compute
  # pair-level SNR downstream.
  conc_mat <- conc_wide %>%
    select(time_round, TowerPosition, mean_conc) %>%
    tidyr::pivot_wider(names_from = TowerPosition,
                       values_from = mean_conc,
                       names_prefix = "C_lev")
  se_mat <- conc_wide %>%
    select(time_round, TowerPosition, se_conc) %>%
    tidyr::pivot_wider(names_from = TowerPosition,
                       values_from = se_conc,
                       names_prefix = "SE_lev")
  lev_cols <- paste0("C_lev", positions)
  lev_cols <- intersect(lev_cols, names(conc_mat))
  if (length(lev_cols) < 2) return(tibble())

  mat <- as.matrix(conc_mat[, lev_cols, drop = FALSE])
  kept <- rowSums(!is.na(mat)) >= 2
  conc_mat <- conc_mat[kept, , drop = FALSE]
  se_mat   <- se_mat[kept, , drop = FALSE]
  mat      <- mat[kept, , drop = FALSE]

  kept_pos <- as.integer(sub("C_lev", "", lev_cols))
  kept_h   <- sapply(kept_pos, function(p)
    height_lookup$height_m[height_lookup$TowerPosition == p])
  ord <- order(kept_h)
  kept_pos <- kept_pos[ord]
  kept_h   <- kept_h[ord]
  lev_cols     <- paste0("C_lev",  kept_pos)
  lev_cols_se  <- paste0("SE_lev", kept_pos)
  mat      <- mat[, lev_cols, drop = FALSE]

  # Join met (this filters rows to those with met data)
  conc_mat <- conc_mat %>% inner_join(met_30, by = "time_round")
  # Re-align se_mat to the remaining time_rounds
  se_mat   <- se_mat[match(conc_mat$time_round, se_mat$time_round), , drop = FALSE]
  se_cols_present <- intersect(lev_cols_se, names(se_mat))
  semat <- as.matrix(se_mat[, se_cols_present, drop = FALSE])
  mat   <- as.matrix(conc_mat[, lev_cols, drop = FALSE])

  nl <- length(kept_h)

  # --- u* threshold filter ---
  ustar_ok <- !is.na(conc_mat$ustar) & conc_mat$ustar >= ustar_threshold
  conc_mat <- conc_mat[ustar_ok, , drop = FALSE]
  mat      <- mat[ustar_ok, , drop = FALSE]
  semat    <- semat[ustar_ok, , drop = FALSE]
  if (nrow(conc_mat) == 0) return(tibble())

  # Build long data frame: one row per (time_round, level).
  # Option A — pair flux labeled at the LOWER sensor of each pair:
  #   F(z_i) = -K(z_i) * (C_{i+1} - C_i) / (z_{i+1} - z_i) * rho_mol
  # The top sensor has no flux (no pair above it).
  per_level <- list()
  disp_h <- if (is.na(canopy_h) || canopy_h < 1) 0 else 0.66 * canopy_h
  rho_mol_vec <- if ("rhoa_kgm3" %in% names(conc_mat) && any(is.finite(conc_mat$rhoa_kgm3)))
    conc_mat$rhoa_kgm3 / 0.028964   # dry air molar density from stage-4 (mol/m³)
  else
    rep(42.3, nrow(conc_mat))

  for (i in seq_len(nl - 1)) {     # i from 1..N-1 — no flux at top sensor
    z_i     <- kept_h[i]
    z_above <- kept_h[i + 1]
    C_i     <- mat[, i]
    C_above <- mat[, i + 1]
    SE_i      <- if (i       <= ncol(semat)) semat[, i]     else rep(NA_real_, length(C_i))
    SE_above  <- if ((i + 1) <= ncol(semat)) semat[, i + 1] else rep(NA_real_, length(C_i))

    dC   <- C_above - C_i
    dCdz <- dC / (z_above - z_i)

    # Pair SNR = |dC| / sqrt(SE_upper^2 + SE_lower^2)
    pair_sd <- sqrt(SE_above^2 + SE_i^2)
    dC_SNR  <- abs(dC) / pair_sd

    # K at z_i: use L_obukhov directly (z_i/zoL would give z_i·L/z_tower, not L)
    L_vec <- if ("L_obukhov" %in% names(conc_mat))
      conc_mat$L_obukhov
    else
      rep(NA_real_, nrow(conc_mat))

    # Per-row α: look up by (site, date) if alpha_fn supplied, else scalar.
    # alpha_fn should accept (site, date) — the sinusoidal version uses DOY
    # extracted from the date; the legacy monthly version ignores the date
    # beyond pulling the month out.
    if (!is.null(alpha_fn) && !is.null(site)) {
      dates <- conc_mat$time_round
      alpha_vec <- vapply(dates, function(d) alpha_fn(site, d), numeric(1))
    } else {
      alpha_vec <- rep(alpha, nrow(conc_mat))
    }

    K_vec <- mapply(function(us, L, al) {
      canopy_K_profile(z_i, canopy_h, us, L, disp_h, al)
    }, conc_mat$ustar, L_vec, alpha_vec)

    per_level[[i]] <- tibble(
      time_round   = conc_mat$time_round,
      TowerPosition = kept_pos[i],
      level_rank   = i,
      z            = z_i,
      z_above      = z_above,
      C            = C_i,
      C_above      = C_above,
      dCdz         = dCdz,
      dC_SNR       = dC_SNR,
      K            = as.numeric(K_vec),
      ustar        = conc_mat$ustar,
      zoL          = conc_mat$zoL,
      L_obukhov    = L_vec,
      F_FG         = -as.numeric(K_vec) * dCdz * rho_mol_vec
    )
  }

  result <- bind_rows(per_level)

  # Back-compat columns so existing plot code keeps working.
  # Under Option A, "flux at sensor z_i" uses the pair (z_i, z_{i+1}).
  # We fill the legacy pair-based columns so that z_upper==max picks out the
  # TOPMOST flux (the highest sensor with a pair above it, i.e. kept_h[N-1])
  # and z_lower==min picks out the BOTTOM flux (kept_h[1]).
  result <- result %>%
    mutate(
      z_lower    = z,
      z_upper    = z_above,
      z_mid      = z,                          # we label F at the lower sensor
      pair_label = paste0(round(z, 1), " m")
    )

  # Join EC reference flux (only relevant for topmost pair)
  if (nrow(ec) > 0) {
    result <- result %>% left_join(ec, by = "time_round")
  } else {
    result$F_EC <- NA_real_
  }

  # --- SNR filter on the concentration gradient ---
  # Drop pairs where the concentration difference is within measurement noise
  # (|dC| < snr_threshold * sigma_dC). Matches parent repo's dConcSNR filter.
  result <- result %>%
    filter(!is.na(dC_SNR), dC_SNR >= snr_threshold)

  # --- Parent-repo "Bad_Eddy" flag (IQR outlier on F_FG, per level) ---
  # Parent repo FLAGS outliers but does NOT remove them from the main dataset
  # (downstream analyses decide whether to use the flag). We match that
  # behavior: add a F_outlier column, keep the row.
  result <- result %>%
    group_by(z) %>%
    mutate(
      .q1  = quantile(F_FG, 0.25, na.rm = TRUE),
      .q3  = quantile(F_FG, 0.75, na.rm = TRUE),
      .iqr = .q3 - .q1,
      F_outlier = F_FG < (.q1 - 1.5 * .iqr) | F_FG > (.q3 + 1.5 * .iqr)
    ) %>%
    ungroup() %>%
    select(-.q1, -.q3, -.iqr)

  # --- Stability flag (parent-repo reporting column, not a filter) ---
  # Stability_500: stable if 0 < L < 500; unstable if -500 < L < 0;
  # neutral if |L| > 500. L_obukhov carried from per_level (stage-4 value).
  result <- result %>%
    mutate(
      Stability_500 = case_when(
        is.na(L_obukhov)                                  ~ NA_character_,
        abs(L_obukhov) > 500                              ~ "neutral",
        L_obukhov > 0 & L_obukhov <= 500                  ~ "stable",
        L_obukhov < 0 & L_obukhov >= -500                 ~ "unstable",
        TRUE                                              ~ NA_character_
      )
    )

  result %>%
    mutate(gas = gas) %>%
    filter(!is.na(F_FG))
}


#' Compute per-layer source/sink as flux divergence plus storage correction
#'
#' S(layer) = [F(z_upper) - F(z_lower)] + storage_flux
#'
#' F(z) comes from compute_layer_fluxes(), labeled at the LOWER sensor of each
#' pair. The flux divergence between two adjacent pairs gives the net transport
#' budget for the enclosed air layer. Adding the layer's storage rate converts
#' this from an apparent transport divergence to a true source/sink estimate.
#'
#' @param flux_df Tibble. Output of compute_layer_fluxes() for one site/gas.
#' @param storage_df Tibble or NULL. Output of compute_layer_storage_flux().
#'   If NULL, S_net = dF (no storage correction; valid for seasonal means).
#' @return Tibble: time_round, pos_lower, pos_upper, z_lower, z_upper,
#'   layer_label, F_lower, F_upper, dF, storage_flux, S_net, S_net_nostorage,
#'   plus ustar / L_obukhov / Stability_500 from the lower-pair row.
compute_source_sink_layer <- function(flux_df, storage_df = NULL) {

  if (is.null(flux_df) || nrow(flux_df) == 0) return(tibble())

  positions <- sort(unique(flux_df$TowerPosition))
  if (length(positions) < 2) return(tibble())

  # Map each position to the next one in the sorted sequence
  next_pos <- c(positions[-1], NA_integer_)
  names(next_pos) <- as.character(positions)

  # Lower boundary flux for each layer (from the pair whose lower sensor = pos_i)
  fl <- flux_df %>%
    select(time_round, pos_lower = TowerPosition,
           z_lower = z, z_upper = z_above,
           F_lower = F_FG, K_lower = K,
           ustar, zoL, L_obukhov, Stability_500) %>%
    mutate(pos_upper = next_pos[as.character(pos_lower)]) %>%
    filter(!is.na(pos_upper))

  # Upper boundary flux (from the pair whose lower sensor = pos_{i+1})
  fu <- flux_df %>%
    select(time_round, pos_upper = TowerPosition, F_upper = F_FG)

  result <- fl %>%
    left_join(fu, by = c("time_round", "pos_upper")) %>%
    mutate(dF = F_upper - F_lower)

  # Join storage correction (optional)
  if (!is.null(storage_df) && nrow(storage_df) > 0) {
    result <- result %>%
      left_join(
        storage_df %>% select(time_round, pos_lower, pos_upper, storage_flux),
        by = c("time_round", "pos_lower", "pos_upper")
      )
  } else {
    result$storage_flux <- NA_real_
  }

  result %>%
    mutate(
      S_net            = dF + coalesce(storage_flux, 0),
      S_net_nostorage  = dF,
      layer_label      = paste0(round(z_lower, 1), "\u2013", round(z_upper, 1), " m")
    ) %>%
    filter(!is.na(F_lower), !is.na(F_upper))
}
