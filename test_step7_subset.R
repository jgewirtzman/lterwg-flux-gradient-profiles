# test_step7_subset.R
# End-to-end sanity test of the pair-resolved footprint pipeline on a small
# multi-site subset (HARV/WREF/KONZ/TOOL) with mock-empty NWI. Validates:
#   - inventory + per-site checkpointing
#   - cross-site pair metrics
#   - regression at per-site, continental levels
#   - mechanism classification with NA Δf_wet
#   - plot generation for all figures that don't strictly need Δf_wet

source("config.R")

sites_to_process <- c("HARV", "WREF", "KONZ", "TOOL")
footprint_sample_n <- 500
do_footprint_mechanism <- TRUE
nwi_buffer_m <- 5000

# ---- Setup directories ----
nwi_dir <- file.path(proj_dir, "data", "nwi")
fp_dir  <- file.path(proj_dir, "output", "footprint")
fp_fig  <- file.path(fig_dir, "footprint")
dir.create(nwi_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fp_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(fp_fig,  showWarnings = FALSE, recursive = TRUE)

# ---- Load site data ----
message("=== Loading site data ===")
all_data <- load_all_sites(sites_to_process, data_dir)

# ---- Compute pair metrics for each site with empty NWI ----
message("=== Computing pair metrics (mock empty NWI) ===")
all_pair_metrics <- list()
for (site in sites_to_process) {
  fp_ckpt <- file.path(fp_dir, paste0(site, "_pair_metrics_CH4.rds"))
  # Use cached pair metrics if available (recompute is slow ~5 min/site)
  if (file.exists(fp_ckpt)) {
    pm <- readRDS(fp_ckpt)
    if (!is.null(pm) && nrow(pm) > 0) {
      all_pair_metrics[[site]] <- pm
      message(sprintf("  %s: %d cached pair rows", site, nrow(pm)))
      next
    }
  }

  if (is.null(all_data[[site]]$aligned) ||
      is.null(all_data[[site]]$aligned$CH4)) {
    message("  ", site, ": no aligned CH4 — skipping")
    next
  }
  lat <- as.numeric(all_data[[site]]$attr$LatTow[1])
  lon <- as.numeric(all_data[[site]]$attr$LonTow[1])
  empty_nwi <- sf::st_sf(geometry = sf::st_sfc(crs = pick_local_aea(lon, lat)))

  pm <- tryCatch(
    compute_footprint_pair_metrics(
      aligned_gas = all_data[[site]]$aligned$CH4,
      attr_df     = all_data[[site]]$attr,
      nwi_proj    = empty_nwi,
      sample_n    = footprint_sample_n,
      r_pct       = 80
    ),
    error = function(e) {
      message("  ", site, " failed: ", e$message); NULL
    }
  )
  if (!is.null(pm) && nrow(pm) > 0) {
    pm$site <- site
    saveRDS(pm, fp_ckpt)
    all_pair_metrics[[site]] <- pm
    message(sprintf("  %s: %d valid pair rows", site, nrow(pm)))
  }
}

pair_df <- bind_rows(all_pair_metrics)
message("\n=== Combined pair metrics ===")
message(sprintf("  total rows: %d   sites: %d",
                nrow(pair_df), length(unique(pair_df$site))))
message("  per-site row counts:")
print(pair_df %>% count(site))
message("  ΔArea summary (m^2):")
print(summary(pair_df$dArea))
message("  dConc summary (CH4 ppm equivalent units):")
print(summary(pair_df$dConc))

# ---- Filter to clean pairs ----
message("\n=== Filtering pairs ===")
# Don't drop NA dfwet since mock NWI gives Δf_wet = 0
pair_clean <- filter_pair_metrics(pair_df, adjacent_only = TRUE,
                                    drop_na_dfwet = FALSE)
message(sprintf("  rows after filtering: %d", nrow(pair_clean)))

# ---- Regression at multiple aggregations ----
message("\n=== Per-site regression ===")
per_site_fit <- fit_per_site_regression(pair_clean, gas = "CH4")
print(per_site_fit %>% select(site, n, slope_dArea, se_dArea, p_dArea, r2_dArea))

message("\n=== Continental pooled regression ===")
# Lmer needs variance in Δf_wet (which is 0 with mock NWI); fall back to dArea only
cont_fit_area <- tryCatch(
  lm(dConc ~ dArea + factor(site), data = pair_clean),
  error = function(e) { message("  failed: ", e$message); NULL }
)
if (!is.null(cont_fit_area)) {
  coefs <- summary(cont_fit_area)$coefficients
  message(sprintf("  ΔArea continental slope: %.3g  (SE %.2g, t=%.2f, p=%.3g)",
                  coefs["dArea", "Estimate"],
                  coefs["dArea", "Std. Error"],
                  coefs["dArea", "t value"],
                  coefs["dArea", "Pr(>|t|)"]))
  message(sprintf("  n=%d, sites=%d", nrow(pair_clean),
                  length(unique(pair_clean$site))))
}

# ---- Mechanism classification ----
message("\n=== Mechanism classification ===")
canopy_for_class <- compute_canopy_metrics_all(all_data)
classified <- classify_mechanism(per_site_fit, canopy_metrics = canopy_for_class)
print(classified %>% select(site, canopy_class, slope_dfwet, p_dfwet,
                             slope_dArea, p_dArea, mech_class))

# ---- Write outputs ----
write.csv(pair_clean %>% select(-any_of("geometry")),
          file.path(tab_dir, "test_pair_metrics_subset_CH4.csv"),
          row.names = FALSE)
write.csv(per_site_fit,
          file.path(tab_dir, "test_per_site_regression_subset_CH4.csv"),
          row.names = FALSE)
write.csv(classified,
          file.path(tab_dir, "test_mechanism_classification_subset_CH4.csv"),
          row.names = FALSE)

# ---- Plot generation tests ----
message("\n=== Plot generation tests ===")
mock_cont_fit <- tibble(
  gas = "CH4",
  slope_dfwet = NA, se_dfwet = NA, t_dfwet = NA,
  slope_dArea = coefs["dArea", "Estimate"],
  se_dArea = coefs["dArea", "Std. Error"],
  t_dArea = coefs["dArea", "t value"],
  n_total = nrow(pair_clean),
  n_sites = length(unique(pair_clean$site)),
  method = "lm(fixed effect site)"
)

# ΔArea scatter (this is the testable one without NWI)
tryCatch({
  p <- plot_dc_vs_dArea(pair_clean, canopy_for_class, gas = "CH4")
  ggsave(file.path(fp_fig, "TEST_dC_vs_dArea_CH4.png"),
         p, width = 10, height = 8, dpi = 130)
  message("  plot_dc_vs_dArea ✓")
}, error = function(e) message("  plot_dc_vs_dArea FAILED: ", e$message))

# Δf_wet scatter (should run but produce flat line since dfwet=0)
tryCatch({
  p <- plot_dc_vs_dfwet(pair_clean, canopy_for_class, mock_cont_fit, gas = "CH4")
  ggsave(file.path(fp_fig, "TEST_dC_vs_dfwet_CH4.png"),
         p, width = 10, height = 8, dpi = 130)
  message("  plot_dc_vs_dfwet ✓ (will be degenerate with mock NWI)")
}, error = function(e) message("  plot_dc_vs_dfwet FAILED: ", e$message))

# Per-site slope forest plot — works on dArea slopes, just labels say dfwet
# Need to manually plot dArea version
tryCatch({
  df <- classified %>%
    filter(!is.na(slope_dArea)) %>%
    arrange(slope_dArea) %>%
    mutate(site = factor(site, levels = site))
  p <- ggplot(df, aes(y = site, x = slope_dArea, color = canopy_class)) +
    geom_vline(xintercept = 0, color = "grey60", linetype = 2) +
    geom_errorbarh(aes(xmin = slope_dArea - se_dArea,
                        xmax = slope_dArea + se_dArea), height = 0) +
    geom_point(size = 3) +
    scale_color_viridis_d(option = "viridis", end = 0.9) +
    labs(x = "per-site slope of dC(CH4) ~ ΔA  (units/m²)",
         y = "site", color = "canopy class",
         title = "ΔArea-only mechanism-2 test, subset of 4 sites") +
    theme_minimal()
  ggsave(file.path(fp_fig, "TEST_per_site_dArea_slopes.png"),
         p, width = 8, height = 5, dpi = 130)
  message("  per-site ΔArea forest plot ✓")
}, error = function(e) message("  per-site ΔArea forest plot FAILED: ", e$message))

# Mechanism map (uses site_coords; classified has lat/lon? no — from nwi_inv)
tryCatch({
  site_coords <- pair_clean %>% select(site) %>% distinct() %>%
    rowwise() %>%
    mutate(lat = as.numeric(all_data[[site]]$attr$LatTow[1]),
           lon = as.numeric(all_data[[site]]$attr$LonTow[1])) %>%
    ungroup()
  p <- plot_mechanism_map(classified, site_coords, gas = "CH4")
  ggsave(file.path(fp_fig, "TEST_mechanism_map_CH4.png"),
         p, width = 12, height = 7, dpi = 130)
  message("  plot_mechanism_map ✓")
}, error = function(e) message("  plot_mechanism_map FAILED: ", e$message))

# f_wet vs height (will be flat at 0 with mock NWI but should still render)
tryCatch({
  p <- plot_fwet_vs_height(pair_clean, canopy_for_class, facet_by_class = FALSE)
  ggsave(file.path(fp_fig, "TEST_fwet_vs_height_CH4.png"),
         p, width = 10, height = 6, dpi = 130)
  message("  plot_fwet_vs_height ✓")
}, error = function(e) message("  plot_fwet_vs_height FAILED: ", e$message))

message("\n=== Test complete. Outputs ===")
message("  CSVs in: ", tab_dir)
message("  Plots in: ", fp_fig)
