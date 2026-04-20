# raw_cache.R
# Fast-iteration cache: save each site's INPUT data (profile_ts, met_ts,
# wind_ts, aligned_data, attr) separately from its computed fluxes.
#
# This lets us re-compute fluxes in minutes (not 30+ min) when only the α
# lookup, filter thresholds, or flux method changes — since raw parsing and
# gap-filling of the RData files is the slow step.
#
# Usage:
#   # first ever run (full parse):
#   build_raw_cache(sites, all_data, gases, cache_dir)
#
#   # subsequent iterations:
#   inputs <- load_raw_cache("HARV", cache_dir)   # instant
#   fluxes <- compute_layer_fluxes(inputs$prof[[gas]], inputs$met, inputs$attr,
#                                    aligned_data = inputs$aligned[[gas]],
#                                    gas = gas, alpha_fn = alpha_fn, site = s)


#' Build or refresh the raw input cache for the given sites.
#'
#' For each site, saves a single RDS: profile_ts per gas, met_ts, wind_ts,
#' attr, and the relevant aligned EC columns (minimal subset to keep RDS
#' small). Existing caches are skipped unless overwrite=TRUE.
build_raw_cache <- function(sites, all_data, gases,
                              cache_dir, overwrite = FALSE) {
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  for (site in sites) {
    fpath <- file.path(cache_dir, paste0(site, "_raw.rds"))
    if (file.exists(fpath) && !overwrite) next

    site_data <- all_data[[site]]
    if (is.null(site_data$raw_9min) || is.null(site_data$attr)) next

    # Met + winds (shared across gases)
    met_ts <- NULL; wind_ts <- NULL
    if (!is.null(site_data$aligned)) {
      met_gas <- intersect(names(site_data$aligned), c("CH4", "CO2", "H2O"))[1]
      if (!is.na(met_gas)) {
        met_ts  <- extract_met_timeseries(site_data$aligned[[met_gas]])
        wind_ts <- extract_wind_profile(site_data$aligned[[met_gas]],
                                           site_data$attr)
      }
    }

    # Profiles per gas
    prof_list <- list()
    aligned_list <- list()
    for (gas in gases) {
      if (gas %in% names(site_data$raw_9min)) {
        prof_list[[gas]] <- prepare_profile_timeseries(site_data$raw_9min, gas)
      }
      # keep only the EC-flux-bearing columns of aligned data, to trim RDS
      if (!is.null(site_data$aligned) && gas %in% names(site_data$aligned)) {
        keep_cols <- intersect(c("timeMid", "timeEnd_A", "FC_turb_interp",
                                  "LE_turb_interp", "ustar_interp"),
                                names(site_data$aligned[[gas]]))
        aligned_list[[gas]] <- site_data$aligned[[gas]][, keep_cols, drop = FALSE]
      }
    }

    saveRDS(list(
      site    = site,
      attr    = site_data$attr,
      prof    = prof_list,
      met     = met_ts,
      winds   = wind_ts,
      aligned = aligned_list
    ), fpath)
    message("  cached raw inputs for ", site)
  }
  invisible(TRUE)
}


#' Load raw input cache for a single site.
load_raw_cache <- function(site, cache_dir) {
  fpath <- file.path(cache_dir, paste0(site, "_raw.rds"))
  if (!file.exists(fpath)) return(NULL)
  readRDS(fpath)
}


#' Fast re-compute path: given an already-built raw cache, compute fluxes
#' for all sites with the current α lookup and filter settings, without
#' re-parsing RData.
#'
#' @param sites Character vector of sites to process.
#' @param gases Character vector of gases.
#' @param cache_dir Directory with <site>_raw.rds files.
#' @param alpha_fn Function(site, date) -> α scalar.
#' @param ustar_threshold / snr_threshold passed to compute_layer_fluxes.
#' @param add_ae_wp If TRUE, also add F_AE / F_WP columns per pair.
#' @return Named list by site, each a named list by gas of flux tibbles.
recompute_fluxes_from_cache <- function(sites, gases, cache_dir,
                                           alpha_fn = NULL,
                                           ustar_threshold = 0.1,
                                           snr_threshold = 3,
                                           add_ae_wp = TRUE) {
  out <- list()
  for (site in sites) {
    inp <- load_raw_cache(site, cache_dir)
    if (is.null(inp)) next
    out[[site]] <- list()
    for (gas in gases) {
      if (is.null(inp$prof[[gas]])) next
      fl <- compute_layer_fluxes(inp$prof[[gas]], inp$met, inp$attr,
                                  aligned_data = inp$aligned[[gas]],
                                  gas = gas,
                                  alpha_fn = alpha_fn, site = site,
                                  ustar_threshold = ustar_threshold,
                                  snr_threshold = snr_threshold)
      if (add_ae_wp) {
        fl <- add_ae_wp_fluxes(fl, inp$met, inp$winds, inp$attr)
      }
      out[[site]][[gas]] <- fl
    }
    message("  recomputed fluxes for ", site)
  }
  out
}
