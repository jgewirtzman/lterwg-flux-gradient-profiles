# load_data.R
# Load downloaded RData files into standardized structures

#' Load all data for a single site
#'
#' @param site Character. NEON site code
#' @param local_dir Character. Local data directory containing site subfolders
#' @return Named list with elements:
#'   $raw_9min  - min9.list: list of 3 DFs (CH4, CO2, H2O) with per-level concentrations
#'   $aligned   - min9Diff.list: list of 3 DFs with paired dConc + met/flux
#'   $attr      - attr.df: site attributes (tower heights, canopy info)
load_site_data <- function(site, local_dir) {

  site_dir <- file.path(local_dir, site)
  if (!dir.exists(site_dir)) {
    stop("Site directory not found: ", site_dir)
  }

  result <- list(raw_9min = NULL, aligned = NULL, attr = NULL)

  # Load raw 9-min concentration data
  f9 <- list.files(site_dir, pattern = paste0(site, "_9min\\.(RData|Rdata|rdata)$"),
                   full.names = TRUE)
  if (length(f9) > 0) {
    env <- new.env()
    load(f9[1], envir = env)
    # The file contains min9.list (list of CH4, CO2, H2O dataframes)
    if ("min9.list" %in% ls(env)) {
      result$raw_9min <- env$min9.list
    } else {
      # Try to find any list object
      objs <- ls(env)
      for (obj in objs) {
        val <- get(obj, envir = env)
        if (is.list(val) && any(c("CH4", "CO2", "H2O") %in% names(val))) {
          result$raw_9min <- val
          break
        }
      }
    }
  }

  # Load aligned concentration-flux data
  fa <- list.files(site_dir, pattern = "aligned_conc_flux_9min\\.(RData|Rdata|rdata)$",
                   full.names = TRUE)
  if (length(fa) > 0) {
    env <- new.env()
    load(fa[1], envir = env)
    # The file contains min9Diff.list (list of CH4, CO2, H2O dataframes)
    if ("min9Diff.list" %in% ls(env)) {
      result$aligned <- env$min9Diff.list
    } else {
      objs <- ls(env)
      for (obj in objs) {
        val <- get(obj, envir = env)
        if (is.list(val) && any(c("CH4", "CO2", "H2O") %in% names(val))) {
          result$aligned <- val
          break
        }
      }
    }
  }

  # Load site attributes
  fat <- list.files(site_dir, pattern = paste0(site, "_attr\\.(RData|Rdata|rdata)$"),
                    full.names = TRUE)
  if (length(fat) > 0) {
    env <- new.env()
    load(fat[1], envir = env)
    if ("attr.df" %in% ls(env)) {
      result$attr <- env$attr.df
    } else {
      objs <- ls(env)
      for (obj in objs) {
        val <- get(obj, envir = env)
        if (is.data.frame(val) && "TowerPosition" %in% names(val)) {
          result$attr <- val
          break
        }
      }
    }
  }

  result
}


#' Load data for multiple sites
#'
#' @param sites Character vector of NEON site codes
#' @param local_dir Character. Local data directory
#' @return Named list (by site) of site data lists
load_all_sites <- function(sites, local_dir) {
  results <- list()
  for (site in sites) {
    message("Loading: ", site)
    results[[site]] <- tryCatch(
      load_site_data(site, local_dir),
      error = function(e) {
        warning("Failed to load ", site, ": ", e$message)
        list(raw_9min = NULL, aligned = NULL, attr = NULL)
      }
    )
  }
  results
}
