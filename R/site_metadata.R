# site_metadata.R
# Per-site metadata, read from the same AmeriFlux/NEON table the methane
# workflow uses, so site classes here match the ones in the methane paper.
#
# Source: lterwg-flux-gradient/metadata/Ameriflux_NEON field-sites.csv
# (IGBP class, Koeppen climate, MAP, MAT, latitude, longitude, elevation).
#
# This replaces hand-coded site labels that cited "NEON site characteristics,
# NLCD" but matched neither: LAJA (1.0 m canopy, IGBP grassland) was labelled
# tropical forest, HEAL (0.3 m, open shrubland) boreal forest, NIWO (0.6 m)
# subalpine forest, BLAN (4.0 m) deciduous forest.
#
# Two classifications are kept on purpose, because they measure different
# things and legitimately disagree at some sites:
#   igbp          - landscape-scale class (what the methane paper uses)
#   canopy_height - NEON DistZaxsCnpy, the canopy where the tower stands
#                   (what vertical-structure analysis needs)
# They conflict at NIWO (IGBP ENF, 0.6 m canopy), CLBJ (IGBP GRA, 16 m) and
# BLAN (IGBP DBF, 4 m). Do not resolve that by picking one; choose the
# classification that matches the question.

suppressPackageStartupMessages({ library(dplyr) })

IGBP_LABELS <- c(
  CRO = "Croplands",
  CVM = "Cropland/Natural Vegetation Mosaic",
  DBF = "Deciduous Broadleaf Forest",
  EBF = "Evergreen Broadleaf Forest",
  ENF = "Evergreen Needleleaf Forest",
  GRA = "Grasslands",
  MF  = "Mixed Forest",
  OSH = "Open Shrublands",
  SAV = "Savannas",
  WET = "Permanent Wetlands"
)

#' Locate the AmeriFlux/NEON site metadata table.
find_site_metadata_file <- function(parent_repo = NULL) {
  env <- Sys.getenv("SITE_METADATA_FILE", unset = "")
  if (nzchar(env)) {
    if (!file.exists(env)) {
      stop("SITE_METADATA_FILE is set to '", env, "', which does not exist.",
           call. = FALSE)
    }
    return(env)
  }
  roots <- unique(c(parent_repo,
                    file.path(dirname(proj_dir), "lterwg-flux-gradient")))
  cand <- file.path(roots, "metadata", "Ameriflux_NEON field-sites.csv")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("Site metadata not found. Looked in:\n  ", paste(cand, collapse = "\n  "),
         "\nSet SITE_METADATA_FILE to the path of 'Ameriflux_NEON field-sites.csv'.",
         call. = FALSE)
  }
  hit[1]
}

#' Load per-site metadata as a tidy tibble, one row per NEON site.
#'
#' @return Tibble: site, ameriflux_id, name, igbp, igbp_label, koeppen,
#'   map_mm, mat_C, lat, lon, elev_m.
load_site_metadata <- function(parent_repo = NULL) {
  f <- find_site_metadata_file(parent_repo)
  m <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  col <- function(nm) {
    if (!nm %in% names(m)) stop("Site metadata is missing column '", nm, "' (", f, ")",
                                call. = FALSE)
    m[[nm]]
  }
  out <- tibble(
    site         = trimws(col("Site_Id NEON")),
    ameriflux_id = trimws(col("Site Id AF")),
    name         = trimws(col("Name")),
    igbp         = trimws(col("Vegetation Abbreviation (IGBP)")),
    koeppen      = trimws(col("Climate Class Abbreviation (Koeppen)")),
    map_mm       = suppressWarnings(as.numeric(col("Mean Average Precipitation (mm)"))),
    mat_C        = suppressWarnings(as.numeric(col("Mean Average Tempurature (degrees C)"))),
    lat          = suppressWarnings(as.numeric(col("Latitude (degrees)"))),
    lon          = suppressWarnings(as.numeric(col("Longitude (degrees)"))),
    elev_m       = suppressWarnings(as.numeric(col("Elevation (m)")))
  ) %>%
    mutate(igbp_label = unname(IGBP_LABELS[igbp])) %>%
    relocate(igbp_label, .after = igbp)

  bad <- out$site[!is.finite(out$lat) | !is.finite(out$lon)]
  if (length(bad)) stop("Site metadata has no coordinates for: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  out
}

#' Longitude (degrees east) for one or more NEON sites.
site_longitude <- function(site, metadata = SITE_METADATA) {
  lon <- metadata$lon[match(site, metadata$site)]
  if (anyNA(lon)) {
    stop("No longitude for site(s): ",
         paste(unique(site[is.na(lon)]), collapse = ", "), call. = FALSE)
  }
  lon
}

#' Approximate local solar hour from a UTC timestamp and a longitude.
#'
#' solar_hour = UTC hour + longitude / 15. This ignores the equation of time
#' (up to about +/-16 minutes), which is negligible for classifying day and
#' night. Timestamps in this project are UTC, so hour() on its own gives the
#' wrong local time by 4-10 hours depending on site.
solar_hour <- function(time_utc, lon) {
  if (!inherits(time_utc, "POSIXct")) stop("solar_hour() needs a POSIXct time.",
                                          call. = FALSE)
  h <- as.numeric(format(time_utc, "%H", tz = "UTC")) +
       as.numeric(format(time_utc, "%M", tz = "UTC")) / 60
  (h + lon / 15) %% 24
}

#' Integer local solar hour (0-23), for binning diel curves.
#'
#' Use this instead of hour(timeMid): timestamps are UTC, so hour() shifts
#' every diel curve by the site's UTC offset (5 h at HARV, 8 h at the western
#' sites), putting the midday peak in the evening.
local_hour <- function(time_utc, lon) {
  as.integer(floor(solar_hour(time_utc, lon))) %% 24L
}

#' Tower longitude from a NEON attribute table (attr.df / attr_df).
#'
#' The tower's own LonTow is the most precise location available for a site,
#' and the per-site functions already receive attr_df.
attr_longitude <- function(attr_df) {
  if (is.null(attr_df) || !"LonTow" %in% names(attr_df)) {
    stop("attr_df has no LonTow column, so the local time of day cannot be ",
         "computed.", call. = FALSE)
  }
  lon <- suppressWarnings(as.numeric(attr_df$LonTow[1]))
  if (!is.finite(lon)) stop("attr_df$LonTow is not a finite number.", call. = FALSE)
  lon
}
