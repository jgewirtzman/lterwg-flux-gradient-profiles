# plot_footprint_maps.R
# Draw real footprint contours, rotated by their real wind direction, on top of
# the real NWI wetland map. Visual verification that the rotation places the
# footprint where the wind says it should, and a look at what each sensor height
# actually samples.
#
# Two figures per site:
#   1. many timestamps, contours coloured by wind direction (does the footprint
#      point INTO the wind, i.e. upwind, as it must?)
#   2. one timestamp, all sensor heights (how does the sampled patch grow?)
#
# Usage:  Rscript plot_footprint_maps.R [SITE ...]

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(lubridate) })

NWI_KM <- 10
nwi_dir <- file.path(proj_dir, "data", sprintf("nwi_%dkm", NWI_KM))
outdir  <- file.path(fig_dir, "footprint_maps")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

args  <- commandArgs(trailingOnly = TRUE)
sites <- if (length(args)) args else c("HARV", "LENO", "GRSM", "UNDE")
set.seed(3)

# Rotate a native Kljun contour onto the map. Verified to reproduce Kljun's own
# rotation exactly: native +x maps to compass bearing wind_dir (i.e. UPWIND,
# toward where the air came from), which as a math angle is 90 - wind_dir.
rotate_contour <- function(xr, yr, wind_dir) {
  phi <- (90 - wind_dir) * pi / 180
  list(x = xr * cos(phi) - yr * sin(phi),
       y = xr * sin(phi) + yr * cos(phi))
}

for (site in sites) {
  message("\n=== ", site)
  met <- load_fluxfoot_met(site); if (is.null(met)) { message("  no met"); next }
  nwi <- load_nwi_for_site(site, nwi_dir); if (is.null(nwi)) { message("  no NWI"); next }
  sd_ <- tryCatch(load_all_sites(site, data_dir), error = function(e) NULL)
  if (is.null(sd_) || is.null(sd_[[site]]$aligned$CH4)) { message("  no data"); next }

  attr_df <- sd_[[site]]$attr
  lat <- as.numeric(attr_df$LatTow[1]); lon <- as.numeric(attr_df$LonTow[1])
  crs <- pick_local_aea(lon, lat)
  nwi <- suppressWarnings(st_make_valid(st_transform(nwi, crs)))

  hl <- attr_df %>% select(TowerPosition, DistZaxsLvlMeasTow) %>% distinct() %>%
    mutate(TowerPosition = as.numeric(TowerPosition),
           z_m = as.numeric(DistZaxsLvlMeasTow)) %>%
    filter(!is.na(z_m), z_m > 0) %>% arrange(z_m)

  df <- sd_[[site]]$aligned$CH4 %>%
    mutate(time_round = round_date(timeMid, "30 minutes")) %>%
    inner_join(met %>% filter(!is.na(wind_dir)), by = "time_round")
  df <- df[which(is.finite(as.numeric(df$ustar_interp)) &
                 is.finite(as.numeric(df$L_obukhov)) &
                 is.finite(as.numeric(df$roughLength_calc))), , drop = FALSE]
  if (!nrow(df)) next

  # sample across wind-direction octants so every bearing is represented
  df$oct <- cut(df$wind_dir, breaks = seq(0, 360, 45), include.lowest = TRUE)
  samp <- df %>% group_by(oct) %>% slice_sample(n = 4) %>% ungroup()

  z_top <- max(hl$z_m); pos_top <- hl$TowerPosition[which.max(hl$z_m)]

  fp_top <- list()
  for (i in seq_len(nrow(samp))) {
    rec <- samp[i, , drop = FALSE]
    sv <- as.numeric(rec$sigma_v_meas[1])
    ust <- as.numeric(rec$ustar_interp[1]); L <- as.numeric(rec$L_obukhov[1])
    if (!is.finite(sv)) sv <- estimate_sigma_v(ust, as.numeric(rec$zoL[1]))
    d <- as.numeric(rec$z_displ_calc[1])
    ub <- paste0("ubar", pos_top)
    ffp <- compute_footprint_kljun(
      zm = pmax(z_top - ifelse(is.na(d), 0, d), 0.5),
      z0 = as.numeric(rec$roughLength_calc[1]),
      umean = if (ub %in% names(rec)) as.numeric(rec[[ub]][1]) else NA_real_,
      h = estimate_blh(L, ust, lat), ol = L, sigmav = sv, ustar = ust, r = 80)
    if (isTRUE(ffp$flag_err == 1) || length(ffp$xr) < 3) next
    rot <- rotate_contour(ffp$xr, ffp$yr, as.numeric(rec$wind_dir[1]))
    fp_top[[length(fp_top) + 1]] <- tibble(
      id = i, wd = as.numeric(rec$wind_dir[1]), x = rot$x, y = rot$y)
  }
  fp_top <- bind_rows(fp_top)
  if (!nrow(fp_top)) { message("  no valid footprints"); next }

  half <- min(NWI_KM * 1000, max(abs(c(fp_top$x, fp_top$y))) * 1.25)
  bb <- st_bbox(c(xmin = -half, ymin = -half, xmax = half, ymax = half), crs = crs)
  cr <- tryCatch(suppressWarnings(st_crop(nwi, bb)), error = function(e) NULL)
  wet_df <- NULL
  if (!is.null(cr) && nrow(cr) > 0) {
    cr <- tryCatch(suppressWarnings(st_collection_extract(cr, "POLYGON")),
                   error = function(e) cr)
    cr <- cr[!st_is_empty(cr), , drop = FALSE]
    if (nrow(cr) > 0) {
      co <- st_coordinates(st_cast(st_geometry(cr), "MULTIPOLYGON", warn = FALSE))
      lev <- intersect(c("L1", "L2", "L3"), colnames(co))
      wet_df <- tibble(x = co[, "X"], y = co[, "Y"],
                       grp = do.call(paste, lapply(lev, function(l) co[, l])))
    }
  }

  # ---------------- FIGURE 1: many wind directions ----------------
  p1 <- ggplot()
  if (!is.null(wet_df))
    p1 <- p1 + geom_polygon(data = wet_df, aes(x, y, group = grp),
                            fill = "#4a90c4", alpha = 0.75, color = NA)
  p1 <- p1 +
    geom_polygon(data = fp_top, aes(x, y, group = id, color = wd),
                 fill = NA, linewidth = 0.5) +
    geom_point(aes(x = 0, y = 0), shape = 4, size = 3, stroke = 1.2) +
    coord_equal(xlim = c(-half, half), ylim = c(-half, half)) +
    scale_color_gradientn(colours = c("#d73027", "#fdae61", "#4575b4",
                                      "#66bd63", "#d73027"),
                          limits = c(0, 360), breaks = seq(0, 360, 90),
                          name = "wind\nfrom (deg)") +
    labs(title = paste0(site, " — top sensor (", round(z_top, 1),
                        " m) footprints, real wind directions"),
         subtitle = "Blue = NWI wetland. x = tower. Each outline is the 80% footprint for one 30-min period.\nA footprint must point INTO the wind it came from: colour should match the side it lies on.",
         x = "east-west distance from tower (m)", y = "north-south (m)") +
    theme_minimal(base_size = 11)
  ggsave(file.path(outdir, paste0(site, "_footprints_by_winddir.png")), p1,
         width = 8.5, height = 7.5, dpi = 150)

  # ---------------- FIGURE 2: all heights, one timestamp ----------------
  rec <- samp[which.max(samp$ustar_interp), , drop = FALSE]
  sv <- as.numeric(rec$sigma_v_meas[1]); ust <- as.numeric(rec$ustar_interp[1])
  L <- as.numeric(rec$L_obukhov[1]); d <- as.numeric(rec$z_displ_calc[1])
  if (!is.finite(sv)) sv <- estimate_sigma_v(ust, as.numeric(rec$zoL[1]))
  hh <- list()
  for (k in seq_len(nrow(hl))) {
    ub <- paste0("ubar", hl$TowerPosition[k])
    ffp <- compute_footprint_kljun(
      zm = pmax(hl$z_m[k] - ifelse(is.na(d), 0, d), 0.5),
      z0 = as.numeric(rec$roughLength_calc[1]),
      umean = if (ub %in% names(rec)) as.numeric(rec[[ub]][1]) else NA_real_,
      h = estimate_blh(L, ust, lat), ol = L, sigmav = sv, ustar = ust, r = 80)
    if (isTRUE(ffp$flag_err == 1) || length(ffp$xr) < 3) next
    rot <- rotate_contour(ffp$xr, ffp$yr, as.numeric(rec$wind_dir[1]))
    hh[[length(hh) + 1]] <- tibble(z = hl$z_m[k], x = rot$x, y = rot$y)
  }
  hh <- bind_rows(hh)
  if (nrow(hh)) {
    h2 <- min(NWI_KM * 1000, max(abs(c(hh$x, hh$y))) * 1.3)
    p2 <- ggplot()
    if (!is.null(wet_df))
      p2 <- p2 + geom_polygon(data = wet_df, aes(x, y, group = grp),
                              fill = "#4a90c4", alpha = 0.75, color = NA)
    p2 <- p2 +
      geom_polygon(data = hh, aes(x, y, group = z, color = z), fill = NA,
                   linewidth = 0.7) +
      geom_point(aes(x = 0, y = 0), shape = 4, size = 3, stroke = 1.2) +
      coord_equal(xlim = c(-h2, h2), ylim = c(-h2, h2)) +
      scale_color_viridis_c(option = "inferno", end = 0.85, name = "sensor\nheight (m)") +
      labs(title = paste0(site, " — every sensor height, one 30-min period"),
           subtitle = sprintf("wind from %.0f deg, u* = %.2f m/s. Higher sensors sample a larger patch further upwind.",
                              as.numeric(rec$wind_dir[1]), ust),
           x = "east-west distance from tower (m)", y = "north-south (m)") +
      theme_minimal(base_size = 11)
    ggsave(file.path(outdir, paste0(site, "_footprints_by_height.png")), p2,
           width = 8, height = 7, dpi = 150)
  }
  message("  wrote 2 figures")
}
message("\nFigures: ", outdir)
