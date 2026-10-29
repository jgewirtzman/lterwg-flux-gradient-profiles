# load_kljun.R
# Load the Kljun et al. (2015) Flux Footprint Prediction (FFP) model.
#
# Kljun, N., Calanca, P., Rotach, M. W., and Schmid, H. P. (2015): A simple
# two-dimensional parameterisation for Flux Footprint Prediction (FFP),
# Geosci. Model Dev., 8, 3695-3713, doi:10.5194/gmd-8-3695-2015.
#
# The R port (v1.42, Copyright (C) 2015-2023 Natascha Kljun) is not vendored
# here: the file carries a copyright notice but no license granting
# redistribution. It is read from the LTER flux-gradient repos instead.
#
# That file cannot be sourced directly. Its function definitions are followed
# by scratch/demo code that fails to parse (line ~1302), so we parse only the
# prefix before the "PLOT" divider and evaluate nothing but function
# definitions, which also guards against any top-level side effects an upstream
# edit might introduce. Loading fails loudly: a silent skip previously left
# every footprint call to break later with a confusing error.

KLJUN_REQUIRED <- c("calc_footprint_FFP", "checkinput",
                    "calc_footprint_FFP_climatology", "checkinput_climat")

#' Locate the Kljun FFP source file.
#'
#' Honours KLJUN_FFP_FILE if set, otherwise looks in the sibling LTER repos.
find_kljun_file <- function(parent_repo = NULL) {
  env <- Sys.getenv("KLJUN_FFP_FILE", unset = "")
  if (nzchar(env)) {
    if (!file.exists(env)) {
      stop("KLJUN_FFP_FILE is set to '", env, "', which does not exist.",
           call. = FALSE)
    }
    return(env)
  }
  roots <- unique(c(parent_repo,
                    file.path(dirname(proj_dir), "lterwg-flux-gradient"),
                    file.path(dirname(proj_dir), "lterwg-flux-gradient-eval")))
  cand <- file.path(roots, "exploratory", "fluxfootprint_kljun.R")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("Kljun FFP source not found. Looked in:\n  ",
         paste(cand, collapse = "\n  "),
         "\nSet KLJUN_FFP_FILE to the path of fluxfootprint_kljun.R.",
         call. = FALSE)
  }
  hit[1]
}

#' Load the Kljun FFP functions into `envir` and verify they work.
load_kljun <- function(parent_repo = NULL, envir = globalenv()) {
  f <- find_kljun_file(parent_repo)
  lines <- readLines(f, warn = FALSE)

  # Function definitions end before the PLOT divider; everything after it is
  # demo code. Fall back to the whole file only if the divider ever vanishes.
  cut_at <- grep("PLOT", lines)[1]
  if (is.na(cut_at)) cut_at <- length(lines) + 1L
  exprs <- tryCatch(parse(text = paste(lines[seq_len(cut_at - 1L)], collapse = "\n")),
                    error = function(e) stop("Could not parse Kljun FFP from ", f,
                                             ": ", conditionMessage(e), call. = FALSE))

  is_fun_def <- function(e) {
    is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
      is.name(e[[2]]) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))
  }
  for (e in exprs) if (is_fun_def(e)) eval(e, envir = envir)

  missing <- KLJUN_REQUIRED[!vapply(KLJUN_REQUIRED, exists, logical(1),
                                    envir = envir, inherits = FALSE)]
  if (length(missing)) {
    stop("Kljun FFP loaded from ", f, " but is missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  # Smoke test on a canned input: an empty or non-finite contour means the
  # loaded code is not the FFP we expect.
  ffp <- suppressWarnings(suppressMessages(get("calc_footprint_FFP", envir = envir)(
    zm = 20, z0 = 0.5, umean = 3, h = 1000, ol = -50, sigmav = 0.8,
    ustar = 0.4, r = 80, nx = 600, rslayer = 1, crop = 1)))
  xr <- tryCatch(ffp$xr[[1]], error = function(e) NULL)
  if (is.null(xr) || length(xr) < 3 || !all(is.finite(xr))) {
    stop("Kljun FFP loaded from ", f, " but its smoke test returned no ",
         "valid 80% contour.", call. = FALSE)
  }
  invisible(f)
}
