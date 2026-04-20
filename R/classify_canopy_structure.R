# classify_canopy_structure.R
# Classify sites by canopy height and sampling coverage, using tower attributes
# (DistZaxsCnpy, DistZaxsTow, DistZaxsLvlMeasTow, LvlMeasTow)
#
# This is a physics- and geometry-based classification that captures what
# actually matters for within-canopy profile analysis: how tall is the canopy
# and how well does the tower sample within vs above it?

#' Compute canopy metrics for one site from attr_df
#'
#' @param attr_df Data frame. Site attributes (one row per tower level)
#' @return Tibble with one row: site, canopy_height, tower_height, n_levels,
#'   n_below_canopy, n_above_canopy, lowest_height, highest_height,
#'   canopy_class, sampling_class
compute_canopy_metrics <- function(attr_df) {

  if (is.null(attr_df) || nrow(attr_df) == 0) return(tibble())

  site <- as.character(attr_df$Site[1])
  canopy_h <- as.numeric(attr_df$DistZaxsCnpy[1])
  tower_h  <- as.numeric(attr_df$DistZaxsTow[1])
  heights  <- as.numeric(attr_df$DistZaxsLvlMeasTow)
  n_lvl    <- length(heights)

  # Count levels below / above canopy
  n_below  <- sum(heights < canopy_h, na.rm = TRUE)
  n_above  <- sum(heights > canopy_h, na.rm = TRUE)

  # Canopy height bins (physics-based: matters for within-canopy structure)
  canopy_class <- case_when(
    is.na(canopy_h) | canopy_h < 1 ~ "No canopy (<1 m)",
    canopy_h <  5                   ~ "Short canopy (1-5 m)",
    canopy_h < 15                   ~ "Medium canopy (5-15 m)",
    canopy_h < 30                   ~ "Tall canopy (15-30 m)",
    TRUE                            ~ "Very tall canopy (>30 m)"
  )

  # Sampling coverage: can we resolve within-canopy structure?
  sampling_class <- case_when(
    n_below == 0                        ~ "Above-canopy only",
    n_below == 1                        ~ "Single below-canopy level",
    n_below >= 2 & n_above >= 1         ~ "Multi-level within-canopy",
    n_above == 0                        ~ "Below-canopy only"
  )

  tibble(
    site            = site,
    canopy_height   = canopy_h,
    tower_height    = tower_h,
    n_levels        = n_lvl,
    n_below_canopy  = n_below,
    n_above_canopy  = n_above,
    lowest_height   = min(heights, na.rm = TRUE),
    highest_height  = max(heights, na.rm = TRUE),
    canopy_class    = canopy_class,
    sampling_class  = sampling_class
  )
}


#' Compute canopy metrics for all sites
#'
#' @param all_data Named list (by site) of site data with $attr
#' @return Tibble with one row per site
compute_canopy_metrics_all <- function(all_data) {

  rows <- list()
  for (site in names(all_data)) {
    attr_df <- all_data[[site]]$attr
    if (is.null(attr_df)) next
    m <- compute_canopy_metrics(attr_df)
    if (nrow(m) > 0) rows[[site]] <- m
  }

  bind_rows(rows) %>%
    mutate(
      canopy_class = factor(canopy_class, levels = c(
        "No canopy (<1 m)",
        "Short canopy (1-5 m)",
        "Medium canopy (5-15 m)",
        "Tall canopy (15-30 m)",
        "Very tall canopy (>30 m)"
      )),
      sampling_class = factor(sampling_class, levels = c(
        "Above-canopy only",
        "Single below-canopy level",
        "Multi-level within-canopy",
        "Below-canopy only"
      ))
    )
}
