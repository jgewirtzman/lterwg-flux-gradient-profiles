# fetch_all_nwi.R
# Populate data/nwi/ with NWI wetland polygons for every NEON site in data/.
#
# Standalone on purpose: it needs only the tower coordinates, which live in the
# small per-site *_attr.Rdata files, so it does not load the multi-GB
# concentration data or source the parent repo the way run_all.R does.
#
# Writes one GeoPackage per site (layers nwi_wgs84 + nwi_proj) plus
# output/tables/nwi_inventory.csv. Re-runs skip sites already on disk unless
# overwrite = TRUE.
#
# Usage:  Rscript fetch_all_nwi.R [SITE ...]     (no args = all sites)

suppressPackageStartupMessages({
  library(sf); library(httr); library(jsonlite); library(dplyr)
})

overwrite <- FALSE
# Buffer radius (m). Overridable so a wider buffer can be fetched without
# clobbering the 5 km cache — footprints at tall towers reach well past 5 km,
# and beyond the buffer wetland is silently counted as dry.
buffer_m  <- as.numeric(Sys.getenv("NWI_BUFFER_M", "5000"))

proj_dir <- normalizePath(".")
data_dir <- file.path(proj_dir, "data")
nwi_dir  <- if (buffer_m == 5000) file.path(data_dir, "nwi") else
            file.path(data_dir, sprintf("nwi_%dkm", round(buffer_m / 1000)))
tab_dir  <- file.path(proj_dir, "output", "tables")
dir.create(nwi_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(proj_dir, "R", "fetch_nwi.R"))

# ---- Site coordinates, read straight from the attr files ----
args  <- commandArgs(trailingOnly = TRUE)
sites <- if (length(args) > 0) args else {
  d <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
  sort(d[!grepl("^nwi", d) & nzchar(d)])
}

site_coords <- function(site) {
  f <- file.path(data_dir, site, paste0(site, "_attr.Rdata"))
  if (!file.exists(f)) return(NULL)
  env <- new.env()
  ok <- tryCatch({ load(f, envir = env); TRUE }, error = function(e) FALSE)
  if (!ok || !"attr.df" %in% ls(env)) return(NULL)
  a <- env$attr.df
  if (!all(c("LatTow", "LonTow") %in% names(a))) return(NULL)
  lat <- suppressWarnings(as.numeric(a$LatTow[1]))
  lon <- suppressWarnings(as.numeric(a$LonTow[1]))
  if (is.na(lat) || is.na(lon)) return(NULL)
  list(lat = lat, lon = lon)
}

message("=== NWI fetch: ", length(sites), " sites, buffer = ", buffer_m, " m ===")

rows <- list()
for (i in seq_along(sites)) {
  site <- sites[i]
  cc <- site_coords(site)
  if (is.null(cc)) {
    message(sprintf("[%2d/%d] %s: no usable coordinates — skipping",
                    i, length(sites), site))
    next
  }

  gpkg <- file.path(nwi_dir, paste0(site, "_nwi.gpkg"))
  if (!overwrite && file.exists(gpkg)) {
    nwi <- tryCatch(st_read(gpkg, layer = "nwi_wgs84", quiet = TRUE),
                    error = function(e) NULL)
    if (!is.null(nwi)) {
      rows[[site]] <- summarize_nwi_inventory(site, nwi, cc$lon, cc$lat, buffer_m)
      message(sprintf("[%2d/%d] %s: cached (%d polygons)",
                      i, length(sites), site, nrow(nwi)))
      next
    }
  }

  t0 <- Sys.time()
  res <- tryCatch(
    fetch_nwi_for_site(site, cc$lon, cc$lat, buffer_m = buffer_m,
                       out_dir = nwi_dir, overwrite = overwrite),
    error = function(e) {
      message(sprintf("[%2d/%d] %s: FETCH FAILED — %s",
                      i, length(sites), site, conditionMessage(e)))
      # NA (not 0) so a service failure is never mistaken for "no wetlands".
      tibble(site = site, n_polygons = NA_integer_,
             total_wet_area_m2 = NA_real_, buffer_area_m2 = (2 * buffer_m)^2,
             fraction_wet = NA_real_, dominant_wetland_type = NA_character_,
             lat = cc$lat, lon = cc$lon, buffer_m = buffer_m)
    }
  )
  rows[[site]] <- res
  if (!is.na(res$n_polygons)) {
    message(sprintf("[%2d/%d] %s: %d polygons, f_wet=%.4f, %s (%.0fs)",
                    i, length(sites), site, res$n_polygons, res$fraction_wet,
                    ifelse(is.na(res$dominant_wetland_type), "-",
                           res$dominant_wetland_type),
                    as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  # Be polite to the USFWS service. At a 10 km buffer a sustained run gets
  # rate-limited partway through (every site from SERC onward failed), so the
  # pause is overridable and should be larger for wide buffers.
  Sys.sleep(as.numeric(Sys.getenv("NWI_SLEEP_S", "1")))
}

inv <- bind_rows(rows)
out_csv <- file.path(tab_dir, sprintf("nwi_inventory%s.csv", if (buffer_m == 5000) "" else sprintf("_%dkm", round(buffer_m/1000))))

# Merge into any existing inventory rather than overwriting it: a single-site
# re-run (retrying one failure) must not wipe the other 46 rows.
if (file.exists(out_csv)) {
  prev <- tryCatch(read.csv(out_csv, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  if (!is.null(prev) && nrow(prev) > 0 && "site" %in% names(prev)) {
    inv <- bind_rows(inv, prev[!prev$site %in% inv$site, , drop = FALSE])
  }
}
inv <- inv[order(inv$site), , drop = FALSE]
write.csv(inv, out_csv, row.names = FALSE)

message("\n=== Summary ===")
message(sprintf("  sites attempted : %d", nrow(inv)))
message(sprintf("  fetch failures  : %d", sum(is.na(inv$n_polygons))))
message(sprintf("  with wetlands   : %d",
                sum(!is.na(inv$n_polygons) & inv$n_polygons > 0)))
message(sprintf("  genuinely dry   : %d",
                sum(!is.na(inv$n_polygons) & inv$n_polygons == 0)))
message("  wrote ", out_csv)
