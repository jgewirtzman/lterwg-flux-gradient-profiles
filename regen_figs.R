suppressPackageStartupMessages({
  library(dplyr); library(lubridate); library(ggplot2); library(patchwork)
})

source("config.R")
for (f in list.files("R", pattern="[.]R$", full.names=TRUE)) source(f)
for (f in list.files("plots", pattern="[.]R$", full.names=TRUE)) source(f)

ckpt_dir <- "output/checkpoints"
fig_dir  <- "output/figures"
mw_dir   <- file.path(fig_dir, "method_workflow")

sites <- sub("[.]rds$", "", list.files(ckpt_dir, "[.]rds$"))
gases <- c("CH4", "CO2", "H2O")

all_fluxes      <- list()
all_profiles    <- list()
all_met         <- list()
all_winds       <- list()
all_source_sink <- list()

cat("Loading checkpoints...\n")
for (site in sites) {
  ckpt <- readRDS(file.path(ckpt_dir, paste0(site, ".rds")))
  all_fluxes[[site]]      <- ckpt$fluxes
  all_profiles[[site]]    <- ckpt$profiles
  all_met[[site]]         <- ckpt$met
  all_winds[[site]]       <- ckpt$winds
  all_source_sink[[site]] <- ckpt$source_sink
}
cat("Loaded", length(sites), "sites\n")

# Load canopy_metrics from the pre-saved CSV (written by full run_all)
canopy_metrics <- tryCatch(
  read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE),
  error = function(e) { cat("canopy_metrics CSV error:", e$message, "\n"); tibble() }
)
cat("canopy_metrics rows:", nrow(canopy_metrics), "\n")

# Load raw site data only for the sites needed by method_workflow (HARV)
mw_sites <- c("HARV")
cat("Loading raw data for method_workflow sites:", paste(mw_sites, collapse=", "), "\n")
all_data <- load_all_sites(mw_sites, data_dir)

cat("\nLayer source/sink figures...\n")
for (gas in gases) {
  gas_fluxes <- lapply(all_fluxes, function(f) f[[gas]])
  gas_fluxes <- gas_fluxes[!sapply(gas_fluxes, is.null)]
  gas_fluxes <- gas_fluxes[sapply(gas_fluxes, function(f) !is.null(f) && nrow(f) > 0)]
  gas_ss <- lapply(all_source_sink, function(x) x[[gas]])
  gas_ss <- gas_ss[!sapply(gas_ss, is.null)]
  gas_ss <- gas_ss[sapply(gas_ss, function(s) !is.null(s) && nrow(s) > 0)]
  if (length(gas_fluxes) == 0) { cat("  skip", gas, "\n"); next }
  tryCatch({
    p <- plot_layer_source_sink(gas_fluxes, canopy_metrics, gas,
                                source_sink_all = if (length(gas_ss) > 0) gas_ss else NULL)
    ggsave(file.path(fig_dir, paste0("multi_site_", gas, "_layer_source_sink.png")),
           p, width = 18, height = 16, dpi = 140, limitsize = FALSE)
    cat("  done:", gas, "\n")
  }, error = function(e) cat("  ERROR", gas, ":", e$message, "\n"))
}

cat("\nMethod workflow figures (HARV only)...\n")
for (gt in list(list("CH4", "Night"), list("CO2", "Day"))) {
  gas <- gt[[1]]; tt <- gt[[2]]
  tryCatch({
    p <- plot_method_workflow("HARV", gas,
                              all_profiles, all_met, all_winds,
                              all_fluxes, all_data,
                              filter_jja_night = TRUE, tod = tt,
                              source_sink_df = all_source_sink[["HARV"]][[gas]])
    ggsave(file.path(mw_dir, paste0("method_workflow_HARV_", gas, ".png")),
           p, width = 20, height = 11, dpi = 150)
    cat("  done: HARV", gas, "\n")
  }, error = function(e) cat("  ERROR HARV", gas, ":", e$message, "\n"))
}
cat("Done.\n")
