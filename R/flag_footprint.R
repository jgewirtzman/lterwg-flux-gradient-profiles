# flag_footprint.R
# Footprint-aware filtering and flagging
#
# Each tower level samples air from a different footprint. Lower levels have
# smaller, tower-proximal footprints; upper levels integrate over larger areas.
# Concentration differences between levels could reflect spatial heterogeneity
# rather than vertical processes. This module flags conditions where footprint
# effects are more or less problematic.

#' Flag footprint reliability conditions on aligned (paired) data
#'
#' Uses stability (z/L) and friction velocity (u*) already present in the
#' aligned concentration-flux data to assess when footprint overlap between
#' paired levels is likely good vs. poor.
#'
#' @param aligned_data Data frame. One gas from min9Diff.list (e.g., $CH4)
#' @param attr_df Data frame. Site attributes
#' @param gas Character. Gas name
#' @param zol_col Character. Column name for stability parameter z/L.
#'   Tries "zoL", then "MO.param" if not found.
#' @param ustar_col Character. Column name for friction velocity.
#'   Default: "ustar_interp"
#' @return Input dataframe with added columns:
#'   fp_overlap ("good"/"moderate"/"poor"), fp_height_ratio, fp_adjacent
flag_footprint_conditions <- function(aligned_data, attr_df, gas = "CH4",
                                       zol_col = NULL, ustar_col = "ustar_interp") {

  if (is.null(aligned_data) || nrow(aligned_data) == 0) return(aligned_data)

  # Find the stability column
  if (is.null(zol_col)) {
    if ("zoL" %in% names(aligned_data)) {
      zol_col <- "zoL"
    } else if ("MO.param" %in% names(aligned_data)) {
      zol_col <- "MO.param"
    } else {
      warning("No stability parameter column found; setting all flags to NA")
      aligned_data$fp_overlap <- NA_character_
      aligned_data$fp_height_ratio <- NA_real_
      aligned_data$fp_adjacent <- NA
      return(aligned_data)
    }
  }

  # Ensure tower heights are numeric
  aligned_data <- aligned_data %>%
    mutate(
      TowerHeight_A = as.numeric(TowerHeight_A),
      TowerHeight_B = as.numeric(TowerHeight_B)
    )

  aligned_data %>%
    mutate(
      # Height ratio: how different are the two measurement heights?
      # Ratio closer to 1 = more similar footprints
      fp_height_ratio = TowerHeight_A / TowerHeight_B,

      # Adjacent levels flag (from dLevelsAminusB)
      fp_adjacent = {
        parts <- strsplit(as.character(dLevelsAminusB), "_")
        sapply(parts, function(p) {
          if (length(p) == 2) abs(as.numeric(p[1]) - as.numeric(p[2])) == 1
          else NA
        })
      },

      # Stability-based overlap assessment
      zol_abs = abs(.data[[zol_col]]),
      ustar_val = .data[[ustar_col]],

      fp_overlap = case_when(
        is.na(zol_abs) | is.na(ustar_val) ~ NA_character_,
        # Well-mixed: near-neutral stability, decent turbulence
        zol_abs < 0.1 & ustar_val > 0.2   ~ "good",
        # Moderately mixed
        zol_abs < 0.5 & ustar_val > 0.1   ~ "moderate",
        # Poorly mixed: strongly stable or unstable, weak turbulence
        TRUE                                ~ "poor"
      ),
      fp_overlap = factor(fp_overlap, levels = c("good", "moderate", "poor"))
    ) %>%
    select(-zol_abs, -ustar_val)
}


#' Filter aligned data by footprint reliability
#'
#' @param flagged_data Data frame. Output of flag_footprint_conditions()
#' @param strictness Character. "strict", "moderate", or "relaxed"
#' @return Filtered data frame
filter_footprint_reliable <- function(flagged_data, strictness = "moderate") {

  switch(strictness,
    "strict" = flagged_data %>%
      filter(fp_overlap == "good", fp_adjacent == TRUE),
    "moderate" = flagged_data %>%
      filter(fp_overlap %in% c("good", "moderate")),
    "relaxed" = flagged_data,
    stop("strictness must be 'strict', 'moderate', or 'relaxed'")
  )
}
