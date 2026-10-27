# load_fluxfoot_met.R
# Extract per-timestamp wind direction and measured sigma_v from the NEON
# 30-min product ({SITE}_30min.Rdata, element FluxFoot = dp04/data/foot/stat).
#
# These are the two Kljun inputs we were previously guessing at:
#   angZaxsErth   -> wind direction (deg from north). Previously NULL, which
#                    left every footprint pinned to one arbitrary bearing.
#   veloYaxsHorSd -> sigma_v measured directly. Previously parameterised as
#                    ~1.9*u* via estimate_sigma_v().
#
# Note distZaxsAbl (boundary-layer height) is a constant 1000 m in this product,
# so it is no better than the existing default and is not used. distObkv
# contains +/-Inf under neutral conditions; the project's own L_obukhov is
# preferable.
#
# Staging dir is deliberately OUTSIDE the Google Drive mount: writing large
# files into the mounted Drive folder silently produces 0-byte files.

suppressPackageStartupMessages({ library(dplyr); library(lubridate) })

FLUXFOOT_DIR <- "~/NEON/fluxgrad_min30"

#' Per-site wind direction and sigma_v on the 30-min grid.
#'
#' @param site Character. NEON site code.
#' @param dir Character. Directory holding {SITE}_30min.zip / extracted data.
#' @return Tibble (time_round, wind_dir, sigma_v_meas, ustar_neon) or NULL.
load_fluxfoot_met <- function(site, dir = FLUXFOOT_DIR) {
  dir <- path.expand(dir)
  rd  <- file.path(dir, "data", site, paste0(site, "_30min.Rdata"))
  if (!file.exists(rd)) {
    zp <- file.path(dir, paste0(site, "_30min.zip"))
    if (!file.exists(zp)) return(NULL)
    tryCatch(utils::unzip(zp, exdir = dir), error = function(e) NULL)
    if (!file.exists(rd)) return(NULL)
  }
  e <- new.env()
  ok <- tryCatch({ load(rd, envir = e); TRUE }, error = function(e) FALSE)
  if (!ok || !length(ls(e))) return(NULL)
  L <- get(ls(e)[1], envir = e)
  ff <- L$FluxFoot
  if (is.null(ff) || !"angZaxsErth" %in% names(ff)) return(NULL)

  out <- tibble(
    time_round   = round_date(as.POSIXct(ff$timeBgn, tz = "UTC"), "30 minutes"),
    wind_dir     = suppressWarnings(as.numeric(ff$angZaxsErth)),
    sigma_v_meas = suppressWarnings(as.numeric(ff$veloYaxsHorSd)),
    ustar_neon   = if ("veloFric" %in% names(ff))
                     suppressWarnings(as.numeric(ff$veloFric)) else NA_real_
  ) %>%
    filter(!is.na(time_round)) %>%
    # wind direction must be a valid bearing; sigma_v must be positive
    mutate(wind_dir     = ifelse(is.finite(wind_dir) & wind_dir >= 0 &
                                   wind_dir <= 360, wind_dir, NA_real_),
           sigma_v_meas = ifelse(is.finite(sigma_v_meas) & sigma_v_meas > 0,
                                 sigma_v_meas, NA_real_)) %>%
    group_by(time_round) %>%
    summarise(wind_dir = first(na.omit(wind_dir)),
              sigma_v_meas = first(na.omit(sigma_v_meas)),
              ustar_neon = first(na.omit(ustar_neon)),
              .groups = "drop")
  out
}


#' Which tall sites currently have usable FluxFoot met on disk.
available_fluxfoot_sites <- function(dir = FLUXFOOT_DIR) {
  dir <- path.expand(dir)
  z <- sub("_30min\\.zip$", "", basename(list.files(dir, pattern = "_30min\\.zip$")))
  d <- if (dir.exists(file.path(dir, "data")))
         list.dirs(file.path(dir, "data"), full.names = FALSE, recursive = FALSE) else character(0)
  sort(unique(c(z, d)))
}
