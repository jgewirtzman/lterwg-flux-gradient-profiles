# fetch_nwi.R
# National Wetlands Inventory (NWI) extraction per NEON site
#
# Pulls polygons from USFWS NWI MapServer (REST/FeatureServer) within a buffer
# around each tower, saves per-site GeoPackage, and builds a coverage / wet-
# fraction inventory across all 47 NEON sites.
#
# Single-source design: NWI is the *only* wetness layer used. Covers CONUS,
# Alaska, Puerto Rico, Hawaii. Static (no ephemeral saturation) — we accept
# that limitation and discuss it.
#
# Endpoint: https://fwspublicservices.wim.usgs.gov/wetlandsmapservice/rest/services/Wetlands/MapServer/0
# Layer 0 = NWI wetland polygons. Returns geojson with WETLAND_TYPE attribute.

suppressPackageStartupMessages({
  library(sf)
  library(httr)
  library(jsonlite)
  library(dplyr)
})

# Default NWI ArcGIS REST query endpoint.
#
# The older fwsprimary.wim.usgs.gov/server/... and www.fws.gov/wetlands/arcgis/...
# hosts are dead as of 2026-08: the former answers queries with ArcGIS error
# bodies, the latter 404s. wetlandsmapservice on fwspublicservices is the live
# service. Kept the old hosts as trailing fallbacks in case they come back.
NWI_ENDPOINTS <- c(
  "https://fwspublicservices.wim.usgs.gov/wetlandsmapservice/rest/services/Wetlands/MapServer/0/query",
  "https://fwsprimary.wim.usgs.gov/server/rest/services/Wetlands/MapServer/0/query",
  "https://www.fws.gov/wetlands/arcgis/rest/services/Wetlands/MapServer/0/query"
)

# Fields on the live service are qualified by source table ("Wetlands.ACRES").
# Requesting the bare names returns an error. We ask for the qualified names
# and strip the prefix off the returned columns so downstream code sees
# ATTRIBUTE / WETLAND_TYPE / ACRES regardless of which host answered.
NWI_OUT_FIELDS <- "Wetlands.ATTRIBUTE,Wetlands.WETLAND_TYPE,Wetlands.ACRES"

# The live service caps a page at 1000 features. Requesting more than the
# server's maxRecordCount is actively harmful here: the paging loop stops when
# a page comes back short, so asking for 2000 and receiving the 1000-feature
# cap looks like "last page" and silently truncates the layer.
NWI_MAX_RECORDS <- 1000

# Floor for the page-size backoff below. Below this, a failure is a real
# failure rather than an over-large response.
NWI_MIN_RECORDS <- 125L

# Retries per page before giving up on an endpoint. The service returns 500
# both for over-large responses and for requests sent too quickly, so the
# retry loop alternates "wait and retry" with "wait and halve the page".
NWI_MAX_ATTEMPTS <- 6L


#' ArcGIS REST reports failures as HTTP 200 with an {"error": ...} JSON body,
#' so httr::stop_for_status() never fires and an unparseable response reads as
#' "no more features". Detect the error body explicitly and return its message.
#'
#' @param txt Character. Raw response body.
#' @return Character error message, or NULL if the body is not an error.
nwi_error_message <- function(txt) {
  if (!grepl("\"error\"", substr(txt, 1, 2000), fixed = TRUE)) return(NULL)
  parsed <- tryCatch(fromJSON(txt, simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$error)) return(NULL)
  msg <- parsed$error$message
  code <- parsed$error$code
  if (is.null(msg) || !nzchar(msg)) msg <- "unspecified ArcGIS error"
  sprintf("ArcGIS error %s: %s", if (is.null(code)) "?" else code, msg)
}


#' Strip the "Wetlands." / "NWI_Wetland_Codes." table qualifier that the live
#' service prepends to field names, so callers see stable bare column names.
normalize_nwi_fields <- function(x) {
  names(x) <- sub("^(Wetlands|NWI_Wetland_Codes)\\.", "", names(x))
  x
}


#' Build an envelope around (lon, lat) of half-width buffer_m metres,
#' expressed in EPSG:4326 (degrees). Uses a local equirectangular approximation
#' that is plenty accurate for buffers of a few km at NEON-tower latitudes
#' (worst case Alaska: ~5% N-S, 10% E-W error in km->deg — negligible because
#' NWI itself is at 1:24,000 nominal scale).
#'
#' @return Named list with xmin, ymin, xmax, ymax (in WGS84 degrees)
build_envelope <- function(lon, lat, buffer_m) {
  dlat <- buffer_m / 111320
  dlon <- buffer_m / (111320 * cos(lat * pi / 180))
  list(xmin = lon - dlon, ymin = lat - dlat,
       xmax = lon + dlon, ymax = lat + dlat)
}


#' Fetch NWI polygons within a buffer of a single point from the USFWS REST
#' API. Pages through results if the server caps response size.
#'
#' @param lon Numeric. Tower longitude (WGS84).
#' @param lat Numeric. Tower latitude (WGS84).
#' @param buffer_m Numeric. Half-width of bounding box (m). Default 5000.
#' @param max_records Integer. Per-page max from server. Default 2000.
#' @return sf object with NWI polygons in EPSG:4326, or empty sf if none found.
fetch_nwi_polygons <- function(lon, lat, buffer_m = 5000,
                                max_records = NWI_MAX_RECORDS) {

  env <- build_envelope(lon, lat, buffer_m)
  geom_json <- toJSON(list(
    xmin = env$xmin, ymin = env$ymin,
    xmax = env$xmax, ymax = env$ymax,
    spatialReference = list(wkid = 4326)
  ), auto_unbox = TRUE)

  base_params <- list(
    geometry       = geom_json,
    geometryType   = "esriGeometryEnvelope",
    inSR           = "4326",
    outSR          = "4326",
    spatialRel     = "esriSpatialRelIntersects",
    outFields      = NWI_OUT_FIELDS,
    returnGeometry = "true",
    f              = "geojson"
  )

  empty_result <- st_sf(WETLAND_TYPE = character(0),
                        ATTRIBUTE    = character(0),
                        ACRES        = numeric(0),
                        geometry     = st_sfc(crs = 4326))

  for (ep in NWI_ENDPOINTS) {
    # Per-endpoint state: a partial failure must not leak an offset or a
    # half-filled feature list into the next endpoint's attempt.
    all_features  <- list()
    result_offset <- 0
    cur_records   <- max_records
    attempt       <- 0L
    ok            <- TRUE

    repeat {
      params <- c(base_params,
                  list(resultOffset = result_offset,
                       resultRecordCount = cur_records))
      resp <- tryCatch(GET(ep, query = params, timeout(180)),
                       error = function(e) e)

      failure <- NULL; oversize <- FALSE
      if (inherits(resp, "error")) {
        failure <- conditionMessage(resp)
      } else if (http_error(resp)) {
        failure <- paste0("HTTP ", status_code(resp))
        # 5xx means the server choked on the response; 4xx means it rejected
        # the request itself (wrong host, bad fields) and a smaller page
        # will not help.
        oversize <- status_code(resp) >= 500
      } else {
        failure <- nwi_error_message(content(resp, as = "text",
                                             encoding = "UTF-8"))
        oversize <- !is.null(failure)
      }

      # Where wetland polygons are geometrically dense (long estuarine
      # shorelines — SERC is the pathological case) a full-size page exceeds
      # the server's response budget and it answers 500. The request is fine;
      # it is only too big. Halve the page and retry the same offset before
      # giving up on the endpoint.
      if (!is.null(failure)) {
        if (oversize && attempt < NWI_MAX_ATTEMPTS) {
          attempt <- attempt + 1L
          # The service also 500s when hit rapidly, so a failure is not
          # necessarily about size. Pause first and retry the same page size
          # once; only shrink the page if it fails again.
          if (attempt %% 2L == 0L) {
            Sys.sleep(5)
            message("  NWI request failed (", failure, ") — retrying same page")
          } else {
            cur_records <- max(NWI_MIN_RECORDS, cur_records %/% 2L)
            Sys.sleep(5)
            message("  NWI request failed (", failure, ") — retrying at ",
                    cur_records, " features/page")
          }
          next
        }
        message("  NWI endpoint failed (", ep, "): ", failure)
        ok <- FALSE; break
      }
      attempt <- 0L   # page succeeded; reset the retry budget

      txt <- content(resp, as = "text", encoding = "UTF-8")
      gj <- tryCatch(read_sf(txt, quiet = TRUE), error = function(e) NULL)
      if (is.null(gj)) {
        message("  NWI endpoint returned unparseable body (", ep, ")")
        ok <- FALSE; break
      }
      if (nrow(gj) == 0) break          # genuinely no (more) features

      all_features[[length(all_features) + 1]] <- normalize_nwi_fields(gj)
      if (nrow(gj) < cur_records) break  # short page => last page
      result_offset <- result_offset + nrow(gj)
    }

    if (!ok) {
      message("  — trying next endpoint")
      next
    }

    if (length(all_features) == 0) {
      # Endpoint answered cleanly and there is genuinely no wetland here.
      attr(empty_result, "endpoint") <- ep
      return(empty_result)
    }
    bound <- do.call(rbind, all_features)
    attr(bound, "endpoint") <- ep
    return(bound)
  }

  # Every endpoint errored — distinct from "answered, found nothing". Signal it
  # so a service outage is never silently recorded as a wetland-free site.
  stop("all NWI endpoints failed for lon=", lon, ", lat=", lat)
}


#' Per-site NWI inventory: fetch polygons, project to local UTM, save to disk,
#' and tally wet-area within the buffer.
#'
#' @param site Character. NEON site code.
#' @param lon Numeric. Tower longitude.
#' @param lat Numeric. Tower latitude.
#' @param buffer_m Numeric. Half-width of bounding box (m). Default 5000.
#' @param out_dir Character. Directory for per-site GeoPackages.
#' @return One-row tibble of inventory metrics for this site.
fetch_nwi_for_site <- function(site, lon, lat, buffer_m = 5000,
                                out_dir = NULL, overwrite = FALSE) {

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    gpkg_path <- file.path(out_dir, paste0(site, "_nwi.gpkg"))
  } else {
    gpkg_path <- NULL
  }

  # Skip if cached
  if (!overwrite && !is.null(gpkg_path) && file.exists(gpkg_path)) {
    nwi <- tryCatch(st_read(gpkg_path, quiet = TRUE), error = function(e) NULL)
    if (!is.null(nwi)) {
      message("  ", site, ": cached (", nrow(nwi), " polygons)")
      return(summarize_nwi_inventory(site, nwi, lon, lat, buffer_m))
    }
  }

  message("  ", site, ": fetching NWI (buffer=", buffer_m, " m)")
  nwi <- fetch_nwi_polygons(lon, lat, buffer_m = buffer_m)

  # Save raw (WGS84) and projected versions to a single GeoPackage
  if (!is.null(gpkg_path) && nrow(nwi) > 0) {
    st_write(nwi, gpkg_path, layer = "nwi_wgs84", delete_dsn = TRUE,
             quiet = TRUE)
    # Also save in local equal-area projection for area math
    utm_crs <- pick_local_aea(lon, lat)
    nwi_proj <- st_transform(nwi, utm_crs)
    st_write(nwi_proj, gpkg_path, layer = "nwi_proj",
             append = TRUE, quiet = TRUE)
  }

  summarize_nwi_inventory(site, nwi, lon, lat, buffer_m)
}


#' Choose a local Albers-equal-area-like projection for tower coordinates.
#' We use Lambert Azimuthal Equal Area centred on the tower — works everywhere
#' (Alaska, tropics) and preserves area for f_wet integration.
pick_local_aea <- function(lon, lat) {
  sprintf("+proj=laea +lat_0=%f +lon_0=%f +datum=WGS84 +units=m +no_defs",
          lat, lon)
}


#' The query envelope in the local equal-area projection, as an sf polygon.
#' Built by projecting the WGS84 envelope corners rather than assuming a
#' (2*buffer_m)^2 square, so its area reflects the actual queried region.
nwi_envelope_proj <- function(lon, lat, buffer_m, crs) {
  env <- build_envelope(lon, lat, buffer_m)
  ring <- matrix(c(env$xmin, env$ymin,
                   env$xmax, env$ymin,
                   env$xmax, env$ymax,
                   env$xmin, env$ymax,
                   env$xmin, env$ymin),
                 ncol = 2, byrow = TRUE)
  poly <- st_sfc(st_polygon(list(ring)), crs = 4326)
  st_transform(poly, crs)
}


#' Summarize NWI features within buffer: total wet area, fraction wet, count.
#' Computed in the local equal-area projection.
#'
#' The service returns every polygon that *intersects* the query envelope, not
#' the polygons clipped to it. Summing those raw areas overstates wet area
#' wherever wetlands are large relative to the buffer — at TOOL it produced a
#' "fraction" of 3.7. We clip to the envelope before integrating.
summarize_nwi_inventory <- function(site, nwi_sf, lon, lat, buffer_m) {
  utm_crs <- pick_local_aea(lon, lat)
  envelope <- nwi_envelope_proj(lon, lat, buffer_m, utm_crs)
  buffer_area_m2 <- as.numeric(st_area(envelope))

  if (nrow(nwi_sf) == 0) {
    return(tibble(
      site = site, n_polygons = 0L,
      total_wet_area_m2 = 0, buffer_area_m2 = buffer_area_m2,
      fraction_wet = 0, dominant_wetland_type = NA_character_,
      lat = lat, lon = lon, buffer_m = buffer_m
    ))
  }

  nwi_proj <- st_transform(nwi_sf, utm_crs)
  nwi_proj <- suppressWarnings(st_make_valid(nwi_proj))
  clipped  <- suppressWarnings(st_intersection(nwi_proj, envelope))

  if (nrow(clipped) == 0) {
    return(tibble(
      site = site, n_polygons = 0L,
      total_wet_area_m2 = 0, buffer_area_m2 = buffer_area_m2,
      fraction_wet = 0, dominant_wetland_type = NA_character_,
      lat = lat, lon = lon, buffer_m = buffer_m
    ))
  }

  clipped$area_m2 <- as.numeric(st_area(clipped))

  dom <- if ("WETLAND_TYPE" %in% names(clipped)) {
    clipped %>%
      st_drop_geometry() %>%
      group_by(WETLAND_TYPE) %>%
      summarise(a = sum(area_m2), .groups = "drop") %>%
      arrange(desc(a)) %>%
      slice(1) %>%
      pull(WETLAND_TYPE)
  } else NA_character_

  wet_area <- sum(clipped$area_m2)

  tibble(
    site = site,
    n_polygons = nrow(clipped),
    total_wet_area_m2 = wet_area,
    buffer_area_m2 = buffer_area_m2,
    fraction_wet = wet_area / buffer_area_m2,
    dominant_wetland_type = if (length(dom) == 0) NA_character_ else dom,
    lat = lat, lon = lon, buffer_m = buffer_m
  )
}


#' Inventory NWI for all sites with loaded attributes.
#'
#' @param all_data Named list (by site) — output of load_all_sites().
#' @param buffer_m Numeric. Half-width of bounding box (m). Default 5000.
#' @param out_dir Character. Where to save per-site GeoPackages.
#' @return Tibble (one row per site) with inventory metrics.
inventory_nwi_all_sites <- function(all_data, buffer_m = 5000, out_dir,
                                    overwrite = FALSE) {
  rows <- list()
  for (site in names(all_data)) {
    attr <- all_data[[site]]$attr
    if (is.null(attr) || !all(c("LatTow", "LonTow") %in% names(attr))) {
      message("  ", site, ": no coords — skipping")
      next
    }
    lat <- as.numeric(attr$LatTow[1])
    lon <- as.numeric(attr$LonTow[1])
    if (is.na(lat) || is.na(lon)) next
    rows[[site]] <- tryCatch(
      fetch_nwi_for_site(site, lon, lat, buffer_m = buffer_m,
                          out_dir = out_dir, overwrite = overwrite),
      error = function(e) {
        message("  ", site, ": fetch failed — ", e$message)
        tibble(site = site, n_polygons = NA_integer_,
               total_wet_area_m2 = NA, buffer_area_m2 = (2*buffer_m)^2,
               fraction_wet = NA, dominant_wetland_type = NA_character_,
               lat = lat, lon = lon, buffer_m = buffer_m)
      }
    )
  }
  bind_rows(rows)
}


#' Load a previously-fetched NWI layer for a site (projected version).
#'
#' @return sf object in local LAEA projection, or NULL if not available.
load_nwi_for_site <- function(site, nwi_dir) {
  gpkg <- file.path(nwi_dir, paste0(site, "_nwi.gpkg"))
  if (!file.exists(gpkg)) return(NULL)
  layers <- st_layers(gpkg)$name
  if ("nwi_proj" %in% layers) {
    st_read(gpkg, layer = "nwi_proj", quiet = TRUE)
  } else if ("nwi_wgs84" %in% layers) {
    st_read(gpkg, layer = "nwi_wgs84", quiet = TRUE)
  } else NULL
}
