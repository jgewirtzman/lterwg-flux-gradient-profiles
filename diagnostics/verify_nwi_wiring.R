# verify_nwi_wiring.R
# Narrow check that real NWI polygons now reach the footprint code and produce
# non-zero, height-varying f_wet. This verifies plumbing only — it is NOT a
# result. Footprints are still unrotated (no wind direction in the aligned
# data), so the Δf_wet magnitudes here are not physically meaningful.
#
# Usage:  Rscript verify_nwi_wiring.R [SITE]     (default HARV)

source("config.R")

site <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a[1] else "HARV" }
nwi_dir <- file.path(proj_dir, "data", "nwi")

message("=== Verifying NWI wiring for ", site, " ===")

nwi_proj <- load_nwi_for_site(site, nwi_dir)
if (is.null(nwi_proj)) stop("no NWI layer on disk for ", site)
message("  NWI polygons loaded: ", nrow(nwi_proj))

all_data <- load_all_sites(site, data_dir)
aligned  <- all_data[[site]]$aligned$CH4
attr_df  <- all_data[[site]]$attr
if (is.null(aligned)) stop("no aligned CH4 for ", site)

# Small sample — we only need to see f_wet come back non-zero and vary.
pm <- compute_footprint_pair_metrics(
  aligned_gas = aligned, attr_df = attr_df, nwi_proj = nwi_proj,
  sample_n = 150, r_pct = 80
)

message("\n=== Result ===")
message("  pair rows: ", nrow(pm))
if (nrow(pm) == 0) stop("no pair rows produced")

message("  f_wet_A summary:"); print(summary(pm$f_wet_A))
message("  f_wet_B summary:"); print(summary(pm$f_wet_B))
message("  dfwet summary:");   print(summary(pm$dfwet))

nz <- sum(pm$f_wet_A > 0 | pm$f_wet_B > 0, na.rm = TRUE)
message(sprintf("\n  rows with any non-zero f_wet: %d / %d", nz, nrow(pm)))
message(sprintf("  rows with dfwet != 0        : %d / %d",
                sum(abs(pm$dfwet) > 0, na.rm = TRUE), nrow(pm)))

# f_wet by measurement height — the quantity the mechanism-2 test leans on.
message("\n  mean f_wet by height (upper sensor of each pair):")
print(pm %>%
        group_by(z = round(z_A, 1)) %>%
        summarise(n = n(), mean_fwet = mean(f_wet_A, na.rm = TRUE),
                  mean_area_km2 = mean(area_A, na.rm = TRUE) / 1e6,
                  .groups = "drop") %>%
        arrange(z))

if (nz == 0) {
  message("\n  FAIL: f_wet still identically zero — polygons are not reaching the intersection.")
} else {
  message("\n  PASS: real wetland area is reaching the footprint intersection.")
}
