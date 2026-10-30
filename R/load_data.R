# load_data.R
# Load per-site RData files into standardized structures.
#
# `local_dir` may be one directory or several. Several directories are treated
# as consecutive NEON download periods -- for the 2021-2026 record,
# NEON_Aligned_Concentrations (Aug 2021 - Jun 2024) and
# NEON_Aligned_Concentrations_202407_202608 (Jul 2024 - Aug 2026) -- and are
# stacked into one continuous record. Stacking mirrors stack.periods() in
# lterwg-flux-gradient/workflows/flow.calc.flux.batch.R: bind rows per gas,
# drop duplicate windows at period boundaries on key columns, sort by time.
# Using the same rule keeps this project's record identical to the one the
# version 3 fluxes were computed from.
#
# Each directory holds per-site subfolders: <dir>/<SITE>/<SITE>_<suffix>.

# Keys that identify one record; duplicates at period boundaries are dropped.
# The aligned keys are exactly those used upstream.
LOAD_KEYS <- list(
  raw_9min = c("timeBgn", "TowerPosition"),
  aligned  = c("timeBgn_A", "timeBgn_B", "TowerPosition_A", "TowerPosition_B")
)

#' Load one RData file and return the object holding per-gas data.
#'
#' Tries `obj_name` first, then any list with CH4/CO2/H2O elements.
load_gas_list <- function(path, obj_name) {
  env <- new.env()
  load(path, envir = env)
  if (obj_name %in% ls(env)) return(get(obj_name, envir = env))
  for (obj in ls(env)) {
    val <- get(obj, envir = env)
    if (is.list(val) && any(c("CH4", "CO2", "H2O") %in% names(val))) return(val)
  }
  NULL
}

#' Find a site's file for `suffix` in one directory, or NA.
#' Anchored on the full name: an unanchored "9min" would also match
#' "_aligned_conc_flux_9min".
site_file <- function(dir, site, suffix) {
  f <- list.files(file.path(dir, site),
                  pattern = paste0("^", site, suffix, "\\.(RData|Rdata|rdata)$"),
                  full.names = TRUE)
  if (length(f)) f[1] else NA_character_
}

#' Stack per-gas lists from several periods, dropping boundary duplicates.
stack_gas_lists <- function(parts, keys) {
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  if (length(parts) == 1) return(parts[[1]])
  gases <- unique(unlist(lapply(parts, names)))
  out <- lapply(gases, function(g) {
    d <- dplyr::bind_rows(lapply(parts, function(x) x[[g]]))
    k <- intersect(keys, names(d))
    if (length(k)) d <- d[!duplicated(d[, k, drop = FALSE]), , drop = FALSE]
    if (length(k)) d <- d[order(d[[k[1]]]), , drop = FALSE]
    d
  })
  names(out) <- gases
  out
}

#' Load all data for a single site.
#'
#' @param site Character. NEON site code.
#' @param local_dir Character. One data directory, or several download-period
#'   directories to stack, each containing site subfolders.
#' @return Named list:
#'   $raw_9min  - min9.list: per-gas data frames of per-level concentrations
#'   $aligned   - min9Diff.list: per-gas data frames of paired dConc + met/flux
#'   $attr      - attr.df: site attributes (tower heights, canopy info)
load_site_data <- function(site, local_dir) {
  have <- local_dir[dir.exists(file.path(local_dir, site))]
  if (!length(have)) {
    stop("No folder for ", site, " in: ", paste(local_dir, collapse = ", "),
         call. = FALSE)
  }

  raw <- lapply(have, function(d) {
    f <- site_file(d, site, "_9min"); if (is.na(f)) NULL else load_gas_list(f, "min9.list")
  })
  ali <- lapply(have, function(d) {
    f <- site_file(d, site, "_aligned_conc_flux_9min")
    if (is.na(f)) NULL else load_gas_list(f, "min9Diff.list")
  })

  # Attributes are not stacked. Use the most recent period's table, and warn
  # if measurement heights differ between periods: stacking then mixes
  # records from different inlet configurations under one height lookup.
  attrs <- lapply(have, function(d) {
    f <- site_file(d, site, "_attr")
    if (is.na(f)) return(NULL)
    env <- new.env(); load(f, envir = env)
    if ("attr.df" %in% ls(env)) return(env$attr.df)
    for (obj in ls(env)) {
      val <- get(obj, envir = env)
      if (is.data.frame(val) && "TowerPosition" %in% names(val)) return(val)
    }
    NULL
  })
  attrs <- Filter(Negate(is.null), attrs)
  attr_df <- if (length(attrs)) attrs[[length(attrs)]] else NULL
  if (length(attrs) > 1) {
    hts <- lapply(attrs, function(a) sort(round(as.numeric(a$DistZaxsLvlMeasTow), 2)))
    if (!all(vapply(hts, identical, logical(1), hts[[1]]))) {
      warning(site, ": measurement heights differ between download periods; ",
              "using the most recent attribute table for all of them.", call. = FALSE)
    }
  }

  list(raw_9min = stack_gas_lists(raw, LOAD_KEYS$raw_9min),
       aligned  = stack_gas_lists(ali, LOAD_KEYS$aligned),
       attr     = attr_df)
}


#' Load data for multiple sites.
#'
#' @param sites Character vector of NEON site codes.
#' @param local_dir Character. One data directory, or several to stack.
#' @return Named list (by site) of site data lists.
load_all_sites <- function(sites, local_dir) {
  results <- list()
  for (site in sites) {
    message("Loading: ", site)
    results[[site]] <- tryCatch(
      load_site_data(site, local_dir),
      error = function(e) {
        warning("Failed to load ", site, ": ", conditionMessage(e), call. = FALSE)
        list(raw_9min = NULL, aligned = NULL, attr = NULL)
      }
    )
  }
  results
}
