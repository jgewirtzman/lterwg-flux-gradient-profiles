# run_all.R
# Master orchestration script for height-resolved concentration analysis
#
# This script runs the full analysis pipeline:
#   1. Download data from Google Drive (toggle off after first run)
#   2. Load data for selected sites
#   3. Process profiles, storage flux, shape analysis, footprint flags
#   4. Multi-tracer comparison (CH4 vs CO2)
#   5. Generate all visualizations
#
# QC NOTE: two analysis branches with DIFFERENT QC requirements:
#
#   Concentration analyses (anomaly, heatmap, vertical profile, storage flux,
#   tracer comparison): use raw profile data — only NEON sensor QC
#   (qfFinl == 0) is applied in prepare_profile_timeseries(). No u*, SNR, or
#   stability filters — these would inappropriately exclude valid low-wind
#   concentration data.
#
#   Flux analyses (K(z) fluxes, source/sink, closure): compute_layer_fluxes()
#   applies ustar >= 0.1, gradient SNR >= 3, and adds Bad_Eddy + Stability
#   flag columns (matching the parent lterwg-flux-gradient workflow).
#
# Usage: source("config.R") first, then source this file or run interactively.

source("config.R")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Sites to process
sites_to_process <- all_sites
# sites_to_process <- forest_sites
# sites_to_process <- c("HARV")         # single site for debugging

# Gases to analyze
gases <- c("CH4", "CO2", "H2O")

# Toggle data download (set FALSE after first successful download)
do_download <- FALSE

# Footprint filter for arrow plots: "strict", "moderate", "relaxed", or NULL
arrow_footprint_filter <- NULL  # show all data; change to "moderate" for filtered

# ============================================================================
# STEP 1: DOWNLOAD DATA
# ============================================================================

if (do_download) {
  message("=== Downloading data from Google Drive ===")
  download_all_sites(
    sites = sites_to_process,
    gdrive_folder_id = gdrive_folder_id,
    local_dir = data_dir,
    file_types = c("9min", "aligned_conc_flux_9min", "attr"),
    overwrite = FALSE
  )
}

# ============================================================================
# STEP 2: LOAD DATA
# ============================================================================

message("=== Loading data ===")
all_data <- load_all_sites(sites_to_process, data_dir)

# ============================================================================
# STEP 3: PROCESS EACH SITE
# ============================================================================

# Storage for results
all_profiles    <- list()  # long-format profile time series
all_diffs       <- list()  # pairwise conc diffs from raw profiles (for arrows)
all_storage     <- list()  # storage flux summaries
all_metrics     <- list()  # profile shape metrics
all_shapes      <- list()  # profile shape classifications
all_met         <- list()  # met/stability time series from aligned data
all_winds       <- list()  # per-level wind profile time series
all_fluxes      <- list()  # layer-by-layer K(z)-based fluxes per site+gas
all_source_sink <- list()  # per-layer source/sink (flux divergence + storage)
all_tracer_comp <- list()  # tracer comparison results

# Checkpoint directory — one RDS per processed site
checkpoint_dir <- file.path(proj_dir, "output", "checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# STEP 3a: FIT SITE-SPECIFIC α FROM WIND PROFILES (before fluxes)
# ============================================================================
# α (Raupach in-canopy attenuation) is fit per site × month from the 2D
# wind profile. Pre-compute it so fluxes use site-specific α instead of
# the textbook default of 3. Evergreen sites should be stable seasonally;
# deciduous sites will show big leaf-on vs leaf-off differences.
message("=== Fitting canopy α from wind profiles ===")
winds_for_fit <- list()
for (site in sites_to_process) {
  ckpt_path <- file.path(checkpoint_dir, paste0(site, ".rds"))
  if (file.exists(ckpt_path)) {
    ck <- readRDS(ckpt_path); winds_for_fit[[site]] <- ck$winds
  } else if (!is.null(all_data[[site]]$aligned) &&
              !is.null(all_data[[site]]$attr)) {
    ga <- intersect(names(all_data[[site]]$aligned),
                     c("CH4", "CO2", "H2O"))[1]
    if (!is.na(ga)) {
      winds_for_fit[[site]] <- extract_wind_profile(
        all_data[[site]]$aligned[[ga]], all_data[[site]]$attr)
    }
  }
}
canopy_metrics_pre <- compute_canopy_metrics_all(all_data)

# Monthly + annual fits (for diagnostics / plotting)
alpha_tbl <- fit_alpha_by_site_month(winds_for_fit, canopy_metrics_pre)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(alpha_tbl, file.path(tab_dir, "alpha_by_site_month.csv"),
          row.names = FALSE)

# SINUSOIDAL α(DOY) — smooth, physically motivated seasonal model, with a
# 4-tier quality-gated fallback cascade:
#   tier 1: sinusoid (r²>=0.25, β>0, α_mean∈[1,5], amp/α_mean<=0.6, phase ok)
#   tier 2: constant α_mean from sinusoid fit (if valid)
#   tier 3: canopy-class median α from tier 1/2 sites
#   tier 4: default α = 3
alpha_sin_tbl <- fit_alpha_sinusoidal(winds_for_fit, canopy_metrics_pre,
                                        monthly_tbl = alpha_tbl)
write.csv(alpha_sin_tbl, file.path(tab_dir, "alpha_sinusoidal.csv"),
          row.names = FALSE)
alpha_lookup <- make_alpha_lookup_sinusoidal(alpha_sin_tbl,
                                               canopy_metrics_pre,
                                               default_alpha = 3,
                                               annual_tbl = alpha_tbl)
alpha_fn       <- alpha_lookup$fn
alpha_tier_tbl <- alpha_lookup$tier_tbl
write.csv(alpha_tier_tbl %>% select(site, canopy_class, tier, alpha_used,
                                       alpha_mean, alpha_amp, phase_doy, r2,
                                       reason),
          file.path(tab_dir, "alpha_tier_assignments.csv"),
          row.names = FALSE)

message("  α tier counts:")
print(table(tier = alpha_tier_tbl$tier, useNA = "ifany"))

ann_summary <- alpha_tbl %>%
  filter(period_type == "annual") %>%
  select(site, alpha, r2, n_obs, note)
message(sprintf("  Fit α for %d forested sites (canopy > 1 m)",
                 nrow(ann_summary)))

# Per-site α fit diagnostic figures (monthly points + sinusoidal smooth)
alpha_fig_dir <- file.path(fig_dir, "alpha_fits")
dir.create(alpha_fig_dir, showWarnings = FALSE, recursive = TRUE)
for (site in names(winds_for_fit)) {
  ch <- canopy_metrics_pre$canopy_height[canopy_metrics_pre$site == site]
  if (length(ch) == 0 || is.na(ch) || ch < 1) next
  tryCatch({
    p <- plot_alpha_fit_diagnostic(site, winds_for_fit[[site]], ch,
                                      alpha_tbl, alpha_sin_tbl,
                                      alpha_tier_tbl)
    ggsave(file.path(alpha_fig_dir, paste0("alpha_fit_", site, ".png")),
           p, width = 11, height = 7.5, dpi = 140)
  }, error = function(e) message("  α fit plot failed for ", site, ": ", e$message))
}
message(sprintf("  Wrote α-fit diagnostic figures to %s", alpha_fig_dir))

for (site in sites_to_process) {
  message("=== Processing: ", site, " ===")

  # Check for existing checkpoint — if present, load and skip
  ckpt_path <- file.path(checkpoint_dir, paste0(site, ".rds"))
  if (file.exists(ckpt_path)) {
    message("  Loading from checkpoint: ", ckpt_path)
    ckpt <- readRDS(ckpt_path)
    all_profiles[[site]]    <- ckpt$profiles
    all_diffs[[site]]       <- ckpt$diffs
    all_storage[[site]]     <- ckpt$storage
    all_metrics[[site]]     <- ckpt$metrics
    all_shapes[[site]]      <- ckpt$shapes
    all_met[[site]]         <- ckpt$met
    all_winds[[site]]       <- ckpt$winds
    all_fluxes[[site]]      <- ckpt$fluxes
    all_source_sink[[site]] <- ckpt$source_sink
    all_tracer_comp[[site]] <- ckpt$tracer_comp
    next
  }

  tryCatch({
    site_data <- all_data[[site]]

    # Check data availability
    has_raw     <- !is.null(site_data$raw_9min)
    has_aligned <- !is.null(site_data$aligned)
    has_attr    <- !is.null(site_data$attr)

    if (!has_attr) {
      message("  No attributes for ", site, " -- skipping")
      next
    }

    # --- Extract met/stability time series from aligned data (once per site) ---
    if (has_aligned) {
      met_gas <- intersect(names(site_data$aligned), c("CH4", "CO2", "H2O"))[1]
      if (!is.na(met_gas)) {
        message("  Extracting met time series (PAR, stability, ustar, wind profile)")
        all_met[[site]]   <- extract_met_timeseries(site_data$aligned[[met_gas]])
        all_winds[[site]] <- extract_wind_profile(site_data$aligned[[met_gas]], site_data$attr)
      }
    }

    for (gas in gases) {
      message("  Gas: ", gas)

      # --- Profile preparation (from raw 9min data) ---
      if (has_raw && gas %in% names(site_data$raw_9min)) {
        prof_ts <- prepare_profile_timeseries(site_data$raw_9min, gas)

        if (nrow(prof_ts) > 0) {
          all_profiles[[site]][[gas]] <- prof_ts

          # Pairwise concentration diffs from raw profile means (for arrows)
          all_diffs[[site]][[gas]] <- calc_profile_diffs(
            prof_ts, site_data$attr, gas,
            met_ts = all_met[[site]],
            group_vars = c("season", "tod")
          )

          # Storage flux
          storage <- calc_storage_flux_level(prof_ts, gas)
          all_storage[[site]][[gas]] <- summarize_storage_flux_diel(storage, site_data$attr)

          # Profile shape classification
          all_shapes[[site]][[gas]] <- classify_profile_shape(prof_ts, site_data$attr, gas)

          # Profile metrics
          all_metrics[[site]][[gas]] <- calc_profile_metrics(prof_ts, site_data$attr, gas)

          # Layer-by-layer fluxes using canopy K(z) model with site-
          # specific α from wind-profile fit (monthly, with annual fallback)
          aligned_for_gas <- if (!is.null(site_data$aligned) &&
                                  gas %in% names(site_data$aligned))
            site_data$aligned[[gas]] else NULL
          all_fluxes[[site]][[gas]] <- compute_layer_fluxes(
            prof_ts, all_met[[site]], site_data$attr,
            aligned_data = aligned_for_gas, gas = gas,
            alpha_fn = alpha_fn, site = site
          )

          # Head-to-head: add F_AE and F_WP columns (parent-repo methods)
          all_fluxes[[site]][[gas]] <- add_ae_wp_fluxes(
            all_fluxes[[site]][[gas]], all_met[[site]],
            all_winds[[site]], site_data$attr
          )

          # Layer storage flux (30-min, per sensor pair, same units as F_FG)
          layer_storage <- compute_layer_storage_flux(
            prof_ts, site_data$attr, all_met[[site]], gas
          )

          # Per-layer source/sink: flux divergence + storage correction
          all_source_sink[[site]][[gas]] <- compute_source_sink_layer(
            all_fluxes[[site]][[gas]], layer_storage
          )
        }
      } else {
        message("    No raw 9min data for ", gas)
      }
    }

    # --- Multi-tracer comparison ---
    if (!is.null(all_profiles[[site]][["CH4"]]) &&
        !is.null(all_profiles[[site]][["CO2"]])) {
      message("  Tracer comparison: CH4 vs CO2")
      all_tracer_comp[[site]] <- compare_tracer_profiles(
        all_profiles[[site]][["CH4"]],
        all_profiles[[site]][["CO2"]],
        site_data$attr
      )
    }

    # --- Save checkpoint ---
    saveRDS(list(
      profiles    = all_profiles[[site]],
      diffs       = all_diffs[[site]],
      storage     = all_storage[[site]],
      metrics     = all_metrics[[site]],
      shapes      = all_shapes[[site]],
      met         = all_met[[site]],
      winds       = all_winds[[site]],
      fluxes      = all_fluxes[[site]],
      source_sink = all_source_sink[[site]],
      tracer_comp = all_tracer_comp[[site]]
    ), ckpt_path)

  }, error = function(e) {
    message("  ERROR processing ", site, ": ", e$message)
  })
}

# ============================================================================
# STEP 4: GENERATE VISUALIZATIONS
# ============================================================================

message("=== Generating visualizations ===")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Per-site plotting loop. Set `do_per_site_plots <- FALSE` in config / CLI
# to skip when iterating fast on cross-site / flux-closure diagnostics.
# ---------------------------------------------------------------------------
do_per_site_plots <- if (exists("do_per_site_plots")) do_per_site_plots else TRUE
if (do_per_site_plots) {
for (site in sites_to_process) {
  tryCatch({
    site_fig_dir <- file.path(fig_dir, site)
    dir.create(site_fig_dir, showWarnings = FALSE, recursive = TRUE)

    for (gas in gases) {
      message("  Plotting: ", site, " ", gas)
      tryCatch({
        # Arrow plot (from raw profile means, not aligned paired dConc)
        if (!is.null(all_diffs[[site]][[gas]])) {
          p <- plot_arrow_dconc(all_diffs[[site]][[gas]], gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_arrows.png")),
                 p, width = 12, height = 8, dpi = 150)
        }

        # Profile plots (from raw 9min)
        if (!is.null(all_profiles[[site]][[gas]])) {
          attr <- all_data[[site]]$attr

          p <- plot_profile_diel(all_profiles[[site]][[gas]], attr, gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_diel_profile.png")),
                 p, width = 14, height = 10, dpi = 150)

          p <- plot_profile_heatmap(all_profiles[[site]][[gas]], attr, gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_heatmap.png")),
                 p, width = 14, height = 10, dpi = 150)

          p <- plot_profile_vertical(all_profiles[[site]][[gas]], attr, gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_vertical_profile.png")),
                 p, width = 12, height = 8, dpi = 150)

          p <- plot_profile_anomaly(all_profiles[[site]][[gas]], attr, gas, site,
                                     met_ts = all_met[[site]])
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_anomaly.png")),
                 p, width = 14, height = 6, dpi = 150)

          p <- plot_anomaly_dotline(all_profiles[[site]][[gas]], attr, gas, site,
                                     met_ts = all_met[[site]])
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_anomaly_dotline.png")),
                 p, width = 14, height = 6, dpi = 150)

          p <- plot_profile_dotline(all_profiles[[site]][[gas]], attr, gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_dotline.png")),
                 p, width = 14, height = 6, dpi = 150)
        }

        # Storage flux
        if (!is.null(all_storage[[site]][[gas]])) {
          p <- plot_storage_flux_diel(all_storage[[site]][[gas]], gas, site)
          ggsave(file.path(site_fig_dir, paste0(site, "_", gas, "_storage_flux.png")),
                 p, width = 12, height = 8, dpi = 150)
        }
      }, error = function(e) {
        message("    ERROR plotting ", site, " ", gas, ": ", e$message)
      })
    }

    # Tracer comparison
    tryCatch({
      if (!is.null(all_tracer_comp[[site]])) {
        p <- plot_tracer_gradient_scatter(all_tracer_comp[[site]], site)
        ggsave(file.path(site_fig_dir, paste0(site, "_tracer_gradient_scatter.png")),
               p, width = 12, height = 10, dpi = 150)

        p <- plot_tracer_correlation_heatmap(all_tracer_comp[[site]], site)
        ggsave(file.path(site_fig_dir, paste0(site, "_tracer_correlation_heatmap.png")),
               p, width = 10, height = 6, dpi = 150)
      }
    }, error = function(e) {
      message("    ERROR tracer comparison ", site, ": ", e$message)
    })
  }, error = function(e) {
    message("  ERROR plotting site ", site, ": ", e$message)
  })
}
}  # end if (do_per_site_plots)

# ============================================================================
# STEP 5: CROSS-SITE SUMMARY (only if multiple sites)
# ============================================================================

if (length(sites_to_process) > 1) {
  message("=== Cross-site summaries ===")

  # Compute canopy metrics once (used by multiple downstream plots)
  canopy_metrics <- compute_canopy_metrics_all(all_data)
  write.csv(canopy_metrics,
            file.path(tab_dir, "canopy_metrics.csv"),
            row.names = FALSE)

  tryCatch({
    p <- plot_canopy_structure(canopy_metrics, all_data)
    ggsave(file.path(fig_dir, "canopy_structure_overview.png"),
           p, width = 16, height = 7, dpi = 150)
  }, error = function(e) message("  ERROR canopy structure plot: ", e$message))

  for (gas in gases) {
    tryCatch({
      # ---------- DISABLED: produced "No data" plots ----------
      # calc_profile_metrics / classify_profile_shape group by exact timeMid,
      # and raw 9-min NEON data measures one level per timestamp — so grouping
      # never yields ≥3 levels. Would need 30-min aggregation first.
      # p <- plot_multi_site_profile_summary(all_metrics, gas)
      # ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_summary.png")), p, ...)
      # p <- plot_multi_site_shape_frequency(all_shapes, gas)
      # ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_shapes.png")), p, ...)

      # ---------- DISABLED: unreadable (188 subplots) ----------
      # Full 4-season x all-site anomaly grid — too dense to read
      # p <- plot_multi_site_anomaly(all_profiles, all_met, all_data, gas,
      #                               jja_only = FALSE)
      # ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_anomaly_grid.png")),
      #        p, width = 14, height = 2 + length(sites_to_process), dpi = 130, limitsize = FALSE)

      # JJA-only compact grid (easier to see summer patterns across sites)
      p <- plot_multi_site_anomaly(all_profiles, all_met, all_data, gas,
                                    jja_only = TRUE)
      ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_anomaly_JJA.png")),
             p, width = 14, height = 10, dpi = 150)

      # Classification-grouped grid (soil sink/source x monotonic/mid-canopy)
      classifications <- classify_sites(all_profiles, all_met, all_data, gas = gas)
      write.csv(classifications,
                file.path(tab_dir, paste0("site_classifications_", gas, ".csv")),
                row.names = FALSE)

      p <- plot_multi_site_anomaly_grouped_grid(all_profiles, all_met, all_data,
                                                  classifications, gas)
      ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_anomaly_classified.png")),
             p, width = 16, height = 12, dpi = 150)

      # Side-by-side: anomaly + wind profile per site (Option A)
      p <- plot_multi_site_anomaly_with_wind(all_profiles, all_met, all_winds,
                                              all_data, classifications, gas)
      ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_anomaly_with_wind.png")),
             p, width = 18, height = 20, dpi = 140)

      # Flux proxy: anomaly / wind speed at each height (Option B)
      p <- plot_multi_site_flux_proxy(all_profiles, all_met, all_winds,
                                       all_data, classifications, gas)
      ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_flux_proxy.png")),
             p, width = 16, height = 12, dpi = 150)

      # --- K(z) based layer fluxes ---
      # Build per-site flux list for this gas
      gas_fluxes <- lapply(all_fluxes, function(f) f[[gas]])
      gas_fluxes <- gas_fluxes[!sapply(gas_fluxes, is.null)]
      gas_fluxes <- gas_fluxes[sapply(gas_fluxes, function(f) nrow(f) > 0)]

      if (length(gas_fluxes) > 0) {
        # Closure check vs EC (only meaningful for CO2 / H2O where F_EC exists)
        if (gas %in% c("CO2", "H2O")) {
          p <- plot_flux_closure(gas_fluxes, canopy_metrics, gas)
          ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_closure_scatter.png")),
                 p, width = 18, height = 16, dpi = 140, limitsize = FALSE)

          p <- plot_flux_closure_summary(gas_fluxes, canopy_metrics, gas)
          ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_closure_summary.png")),
                 p, width = 10, height = 10, dpi = 150)
        }

        # Lagrangian internal closure: top-pair flux vs bottom-up reconstruction.
        # Works for all gases including CH4 (no EC needed).
        p <- plot_flux_closure_lagrangian(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_closure_lagrangian.png")),
               p, width = 18, height = 16, dpi = 140, limitsize = FALSE)

        # Lagrangian closure SUMMARY (one bar per site — all gases)
        p <- plot_flux_closure_lagrangian_summary(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_closure_lagrangian_summary.png")),
               p, width = 10, height = 10, dpi = 150)

        # Vertical flux profile (all gases)
        p <- plot_flux_vertical_profile(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_flux_profile.png")),
               p, width = 18, height = 16, dpi = 140, limitsize = FALSE)

        # Layer source/sink profile
        p <- plot_layer_source_sink(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_layer_source_sink.png")),
               p, width = 18, height = 16, dpi = 140, limitsize = FALSE)

        # ---- Synthesis plots (readable at a glance) ----
        # 1. By-canopy-class median profile with IQR ribbon (5 panels, normalized height)
        p <- plot_flux_profile_by_class(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("synth_", gas, "_profile_by_class.png")),
               p, width = 16, height = 6, dpi = 150)

        # 3. Ranked bar chart: one bar per site, sorted within canopy class
        p <- plot_flux_ranked_bars(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("synth_", gas, "_ranked_bars.png")),
               p, width = 10, height = 14, dpi = 150)

        # 4. Exemplar sites — detailed profile (9 archetype sites)
        p <- plot_flux_profile_exemplar(gas_fluxes, canopy_metrics, gas)
        ggsave(file.path(fig_dir, paste0("synth_", gas, "_exemplar_profiles.png")),
               p, width = 14, height = 12, dpi = 150)

        # 5. Per-ecosystem-type profile grids (all sites in each type),
        #    JJA-night version + all-data version. Shaded canopy band.
        ecotype_dir <- file.path(fig_dir, "by_ecotype")
        dir.create(ecotype_dir, showWarnings = FALSE, recursive = TRUE)
        for (eco_type in unique(ecosystem_types)) {
          n_sites <- sum(ecosystem_types == eco_type &
                          names(ecosystem_types) %in% names(gas_fluxes))
          if (n_sites == 0) next
          eco_slug <- gsub("[^A-Za-z0-9]+", "_", eco_type)
          h <- max(6, 2 + ceiling(n_sites / 4) * 3.2)

          for (cond in c("JJAnight", "all")) {
            p <- plot_flux_profile_by_ecotype(
              gas_fluxes, canopy_metrics, ecosystem_types,
              eco_type = eco_type, gas = gas,
              filter_jja_night = (cond == "JJAnight"))
            ggsave(file.path(ecotype_dir,
                             paste0("ecotype_", eco_slug, "_", gas, "_", cond, ".png")),
                   p, width = 14, height = h, dpi = 140, limitsize = FALSE)
          }
        }

        # 6. Ground-level vs ecosystem-level flux scatter
        #    (per gas × condition × statistic — 4 plots per gas)
        for (cond in c("JJAnight", "all")) {
          for (stat_kind in c("median", "mean")) {
            p <- plot_flux_ground_vs_ecosystem(
              gas_fluxes, canopy_metrics, gas,
              filter_jja_night = (cond == "JJAnight"),
              stat = stat_kind)
            ggsave(file.path(fig_dir,
                             paste0("synth_", gas, "_ground_vs_ecosystem_",
                                    cond, "_", stat_kind, ".png")),
                   p, width = 10, height = 10, dpi = 150)
          }
        }
      }
    }, error = function(e) {
      message("  ERROR cross-site summary ", gas, ": ", e$message)
    })
  }

  # --- Per-ecosystem-type cross-site plots (DISABLED for speed) ---
  # These produce canopy-class subset duplicates of the multi-site figures.
  # Re-enable when zooming in on a specific ecosystem. Keeps the pipeline fast.
  if (FALSE) {
  message("=== Per-ecosystem cross-site summaries ===")

  # Use canopy-structure class as the grouping variable
  eco_map <- setNames(as.character(canopy_metrics$canopy_class),
                       canopy_metrics$site)

  # Loop over each ecosystem type
  for (eco_type in unique(eco_map)) {
    eco_sites <- names(eco_map)[eco_map == eco_type]
    if (length(eco_sites) < 2) next  # need at least 2 sites for a grid

    # Filter the big lists to this ecosystem
    eco_profiles <- all_profiles[eco_sites]
    eco_met      <- all_met[eco_sites]
    eco_winds    <- all_winds[eco_sites]
    eco_data     <- all_data[eco_sites]

    # Drop NULL entries (sites that failed to load)
    eco_profiles <- eco_profiles[!sapply(eco_profiles, is.null)]
    if (length(eco_profiles) < 2) next

    eco_slug <- gsub("[^A-Za-z0-9]+", "_", eco_type)
    message("  Ecosystem: ", eco_type, " (", length(eco_profiles), " sites)")

    eco_dir <- file.path(fig_dir, "by_ecosystem", eco_slug)
    dir.create(eco_dir, showWarnings = FALSE, recursive = TRUE)

    for (gas in gases) {
      tryCatch({
        # JJA compact grid
        p <- plot_multi_site_anomaly(eco_profiles, eco_met, eco_data, gas,
                                      jja_only = TRUE)
        ggsave(file.path(eco_dir, paste0(eco_slug, "_", gas, "_anomaly_JJA.png")),
               p, width = 14, height = max(6, 1.5 * length(eco_profiles)),
               dpi = 140, limitsize = FALSE)

        # Classifications within this ecosystem
        eco_cls <- classify_sites(eco_profiles, eco_met, eco_data, gas = gas)
        write.csv(eco_cls,
                  file.path(tab_dir, paste0(eco_slug, "_classifications_", gas, ".csv")),
                  row.names = FALSE)

        # Classification-grouped grid
        p <- plot_multi_site_anomaly_grouped_grid(eco_profiles, eco_met, eco_data,
                                                    eco_cls, gas)
        ggsave(file.path(eco_dir, paste0(eco_slug, "_", gas, "_anomaly_classified.png")),
               p, width = 14, height = max(6, 2 + length(eco_profiles) * 0.8),
               dpi = 140, limitsize = FALSE)

        # Anomaly + wind side-by-side
        p <- plot_multi_site_anomaly_with_wind(eco_profiles, eco_met, eco_winds,
                                                eco_data, eco_cls, gas)
        ggsave(file.path(eco_dir, paste0(eco_slug, "_", gas, "_anomaly_with_wind.png")),
               p, width = 16, height = max(8, 2 + length(eco_profiles) * 1.2),
               dpi = 140, limitsize = FALSE)

        # Flux proxy
        p <- plot_multi_site_flux_proxy(eco_profiles, eco_met, eco_winds,
                                         eco_data, eco_cls, gas)
        ggsave(file.path(eco_dir, paste0(eco_slug, "_", gas, "_flux_proxy.png")),
               p, width = 14, height = max(6, 2 + length(eco_profiles) * 0.8),
               dpi = 140, limitsize = FALSE)
      }, error = function(e) {
        message("    ERROR ", eco_slug, " ", gas, ": ", e$message)
      })
    }
  }
  }  # end if(FALSE) per-ecosystem block
}

# ============================================================================
# STEP 6: PER-SITE METHOD-WORKFLOW FIGURES
# Two per site: CH4 (JJA Night) and CO2 (JJA Day)
# ============================================================================
message("=== Generating per-site method-workflow figures ===")
mw_dir <- file.path(fig_dir, "method_workflow")
dir.create(mw_dir, showWarnings = FALSE, recursive = TRUE)

for (site in sites_to_process) {
  if (is.null(all_fluxes[[site]])) next
  for (gt in list(list("CH4", "Night"), list("CO2", "Day"))) {
    gas <- gt[[1]]; tt <- gt[[2]]
    if (is.null(all_fluxes[[site]][[gas]]) ||
        is.null(all_profiles[[site]][[gas]])) next
    tryCatch({
      p <- plot_method_workflow(site, gas,
                                  all_profiles, all_met, all_winds,
                                  all_fluxes, all_data,
                                  filter_jja_night = TRUE, tod = tt)
      ggsave(file.path(mw_dir,
                       paste0("method_workflow_", site, "_", gas, ".png")),
             p, width = 20, height = 11, dpi = 150)
    }, error = function(e) {
      message("  method workflow ", site, " ", gas, ": ", e$message)
    })
  }
}
message(sprintf("  Wrote %d method-workflow figures to %s",
                 length(list.files(mw_dir, pattern = "\\.png$")), mw_dir))

message("=== Done! Figures saved to: ", fig_dir, " ===")
