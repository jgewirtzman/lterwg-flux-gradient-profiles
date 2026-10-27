# make_footprint_summary_figs.R
# Footprint / wetland summary under the symmetric-wind (rotationally-averaged)
# assumption — see R/footprint_radial.R.
#
# Runs over every site that has an NWI layer. Per-site results are checkpointed
# to output/footprint/radial/, so re-runs are incremental and the figure code
# can be re-run without recomputing footprints.
#
# Outputs:
#   output/tables/footprint_fwet_by_height.csv   f_wet + reach per site x height
#   output/tables/footprint_fwet_summary.csv     one row per site (top-bottom)
#   SUMMARY_fwet_delta_by_site.png    cross-site: does the top sensor see more wetland?
#   SUMMARY_fwet_by_height_all.png    f_wet(z) for every site
#   SUMMARY_fwet_by_height.png        f_wet(z), exemplar sites, readable
#   SUMMARY_wetland_maps.png          wetland maps + footprint rings, exemplars
#   SUMMARY_radial_wetland_profile.png g(r), exemplars
#   SUMMARY_footprint_reach_by_height.png
#
# These are descriptive: they never touch dConc, so they are not the
# mechanism-2 test.
#
# Usage:  Rscript make_footprint_summary_figs.R [SITE ...]     (default: all)

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })

n_ts   <- 25     # timestamps sampled per site
n_r    <- 90     # radial bins
n_ang  <- 180    # angular samples per radius
SAT_HI <- 0.95   # f_wet above this at every height = saturated, cannot discriminate
SAT_LO <- 0.005  # f_wet below this at every height = effectively dry

# Buffer radius to use, in km. Footprints at tall towers reach past 5 km, and
# beyond the buffer wetland is silently counted as dry, so a wider buffer is
# the fix for the 14 sites flagged reach_beyond_buffer. Checkpoints and figures
# are kept separate per buffer so results are comparable rather than clobbered.
NWI_KM   <- as.numeric(Sys.getenv("NWI_KM", "5"))
R_CAP    <- NWI_KM * 1000
nwi_dir  <- if (NWI_KM == 5) file.path(proj_dir, "data", "nwi") else
            file.path(proj_dir, "data", sprintf("nwi_%dkm", round(NWI_KM)))
sfx      <- if (NWI_KM == 5) "" else sprintf("_%dkm", round(NWI_KM))
fp_fig   <- file.path(fig_dir, "footprint")
ckpt_dir <- file.path(proj_dir, "output", "footprint", paste0("radial", sfx))
dir.create(fp_fig,   showWarnings = FALSE, recursive = TRUE)
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)

# Exemplars get the detailed per-panel figures; 46 facets is unreadable.
exemplars <- c("BARR", "TOOL", "UNDE", "DELA", "HARV", "KONZ")

args <- commandArgs(trailingOnly = TRUE)
sites <- if (length(args) > 0) args else {
  have_nwi <- sub("_nwi\\.gpkg$", "",
                  list.files(nwi_dir, pattern = "_nwi\\.gpkg$"))
  sort(intersect(have_nwi,
                 list.dirs(data_dir, full.names = FALSE, recursive = FALSE)))
}
message("Sites to process: ", length(sites))

set.seed(42)

# ---------------------------------------------------------------- COMPUTE
for (site in sites) {
  ck <- file.path(ckpt_dir, paste0(site, ".rds"))
  if (file.exists(ck)) { message("[cached] ", site); next }
  message("\n=== ", site, " ===")

  nwi_proj <- load_nwi_for_site(site, nwi_dir)
  if (is.null(nwi_proj)) { message("  no NWI — skipping"); next }

  sd <- tryCatch(load_all_sites(site, data_dir), error = function(e) NULL)
  if (is.null(sd) || is.null(sd[[site]]$aligned$CH4)) {
    message("  no aligned CH4 — skipping"); next
  }
  attr_df <- sd[[site]]$attr
  lat <- as.numeric(attr_df$LatTow[1]); lon <- as.numeric(attr_df$LonTow[1])
  target_crs <- pick_local_aea(lon, lat)
  nwi_proj <- st_transform(nwi_proj, target_crs)

  hl <- attr_df %>%
    select(TowerPosition, DistZaxsLvlMeasTow) %>% distinct() %>%
    mutate(TowerPosition = as.numeric(TowerPosition),
           z_m = as.numeric(DistZaxsLvlMeasTow)) %>%
    filter(!is.na(z_m), z_m > 0) %>% arrange(z_m)
  if (nrow(hl) < 2) { message("  <2 heights — skipping"); next }

  df <- sd[[site]]$aligned$CH4
  df <- df[which(is.finite(as.numeric(df$ustar_interp)) &
                 is.finite(as.numeric(df$L_obukhov)) &
                 is.finite(as.numeric(df$roughLength_calc))), , drop = FALSE]
  if (nrow(df) == 0) { message("  no usable met rows — skipping"); next }
  df <- df[sample.int(nrow(df), min(n_ts, nrow(df))), , drop = FALSE]

  message("  footprints: ", nrow(df), " timestamps x ", nrow(hl), " heights")
  contours <- list()
  for (i in seq_len(nrow(df))) {
    fps <- tryCatch(
      compute_footprints_all_heights(met_row = df[i, , drop = FALSE],
                                     heights = hl$z_m,
                                     tower_positions = hl$TowerPosition,
                                     lat = lat, r = 80),
      error = function(e) NULL)
    if (is.null(fps)) next
    for (k in seq_len(nrow(fps))) {
      f <- fps$footprint[[k]]
      if (isTRUE(f$flag_err == 1) || length(f$xr) < 3) next
      contours[[length(contours) + 1]] <- list(z = fps$z_m[k],
                                               xr = f$xr, yr = f$yr)
    }
  }
  if (length(contours) == 0) { message("  no valid footprints"); next }

  reach <- max(vapply(contours, function(c) max(sqrt(c$xr^2 + c$yr^2)), 0))
  r_max <- min(R_CAP, max(reach * 1.05, 200))
  beyond <- reach > R_CAP
  if (beyond) message(sprintf("  NOTE: footprints reach %.0f m, beyond the %d m buffer",
                              reach, R_CAP))
  grid <- make_polar_grid(r_max, n_r = n_r, n_ang = n_ang)

  g <- wetland_radial_profile(nwi_proj, grid, target_crs)

  wl <- list()
  for (cc in contours) {
    w <- footprint_radial_weight(list(xr = cc$xr, yr = cc$yr, flag_err = 0L),
                                 grid, target_crs)
    if (is.null(w)) next
    key <- sprintf("%.2f", cc$z)
    wl[[key]] <- if (is.null(wl[[key]])) list(z = cc$z, sum = w, n = 1L) else
      list(z = cc$z, sum = wl[[key]]$sum + w, n = wl[[key]]$n + 1L)
  }
  if (length(wl) == 0) { message("  no radial weights"); next }

  ph <- bind_rows(lapply(wl, function(e) {
    w <- e$sum / e$n
    tibble(site = site, z = e$z, n_fp = e$n,
           f_wet = fwet_rotational(w, g, grid$r),
           r50   = footprint_radius_quantile(w, grid$r, 0.5),
           r90   = footprint_radius_quantile(w, grid$r, 0.9))
  })) %>% arrange(z)
  ph$f_wet_buffer  <- fwet_rotational(rep(1, length(grid$r)), g, grid$r)
  ph$reach_beyond_buffer <- beyond

  half <- min(R_CAP, max(ph$r90, na.rm = TRUE) * 1.35)
  cropped <- NULL
  if (site %in% exemplars) {
    bb <- st_bbox(c(xmin = -half, ymin = -half, xmax = half, ymax = half),
                  crs = target_crs)
    cr <- tryCatch(suppressWarnings(st_crop(st_make_valid(nwi_proj), bb)),
                   error = function(e) NULL)
    if (!is.null(cr) && nrow(cr) > 0) {
      # Cropping polygons can yield GEOMETRYCOLLECTIONs; keep polygons only so
      # the geometry has one type and can be flattened for plotting.
      cr <- tryCatch(suppressWarnings(st_collection_extract(cr, "POLYGON")),
                     error = function(e) cr)
      cr <- cr[!st_is_empty(cr), , drop = FALSE]
      if (nrow(cr) > 0) cropped <- st_cast(st_geometry(cr), "MULTIPOLYGON",
                                           warn = FALSE)
    }
  }

  saveRDS(list(ph = ph, rad = tibble(site = site, r = grid$r, g = g),
               half = half, nwi = cropped), ck)
  print(as.data.frame(ph %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
}

# ---------------------------------------------------------------- ASSEMBLE
cks <- list.files(ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(cks) == 0) stop("no per-site results")
res <- lapply(cks, readRDS)
names(res) <- sub("\\.rds$", "", basename(cks))

ph_all  <- bind_rows(lapply(res, `[[`, "ph"))
rad_all <- bind_rows(lapply(res, `[[`, "rad"))
half_of <- lapply(res, `[[`, "half")
nwi_of  <- lapply(res, `[[`, "nwi")

write.csv(ph_all, file.path(tab_dir, paste0("footprint_fwet_by_height", sfx, ".csv")),
          row.names = FALSE)

# One row per site: what the pair test would actually see, top vs bottom.
summ <- ph_all %>%
  filter(!is.na(f_wet)) %>%
  group_by(site) %>%
  arrange(z, .by_group = TRUE) %>%
  summarise(
    n_heights   = n(),
    z_bot       = first(z),      z_top      = last(z),
    fwet_bot    = first(f_wet),  fwet_top   = last(f_wet),
    n_fp_bot    = first(n_fp),   n_fp_top   = last(n_fp),
    r50_top     = last(r50),
    f_wet_buffer = first(f_wet_buffer),
    beyond_buffer = first(reach_beyond_buffer),
    .groups = "drop") %>%
  mutate(
    delta_fwet = fwet_top - fwet_bot,
    status = case_when(
      pmin(fwet_bot, fwet_top) > SAT_HI ~ "saturated wet",
      pmax(fwet_bot, fwet_top) < SAT_LO ~ "effectively dry",
      TRUE                              ~ "can discriminate"),
    direction = case_when(
      status != "can discriminate" ~ "n/a",
      delta_fwet > 0               ~ "top sees more wetland",
      TRUE                         ~ "top sees less wetland")) %>%
  arrange(desc(delta_fwet))
write.csv(summ, file.path(tab_dir, paste0("footprint_fwet_summary", sfx, ".csv")),
          row.names = FALSE)

message("\n=== Cross-site summary ===")
print(as.data.frame(summ %>% select(site, z_bot, z_top, fwet_bot, fwet_top,
                                    delta_fwet, status, direction) %>%
                    mutate(across(where(is.numeric), ~round(.x, 4)))))
message("\nstatus counts:"); print(table(summ$status))
message("direction among discriminating sites:")
print(table(summ$direction[summ$status == "can discriminate"]))

lab_of <- setNames(sprintf("%s  (buffer f_wet %.2f)", summ$site, summ$f_wet_buffer),
                   summ$site)

# ------------------------------------------------------- FIG: cross-site delta
sd_df <- summ %>% mutate(site = factor(site, levels = rev(site)))
p0 <- ggplot(sd_df, aes(x = delta_fwet, y = site, color = status)) +
  geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  geom_segment(aes(x = 0, xend = delta_fwet, yend = site), linewidth = 0.5) +
  geom_point(aes(size = f_wet_buffer)) +
  scale_color_manual(values = c("can discriminate" = "#1b6ca8",
                                "saturated wet"    = "#c0392b",
                                "effectively dry"  = "#999999")) +
  scale_size_continuous(range = c(1.2, 4.5), name = "buffer\nf_wet") +
  labs(title = "Does the top sensor see more wetland than the bottom one?",
       subtitle = paste("f_wet(top) - f_wet(bottom), rotationally averaged over wind direction.",
                        "\nRed = wetland everywhere at every height; grey = essentially no wetland. Neither can test the mechanism."),
       x = "delta f_wet  (top - bottom sensor)", y = NULL, color = NULL) +
  theme_minimal(base_size = 10)
ggsave(file.path(fp_fig, paste0("SUMMARY_fwet_delta_by_site", sfx, ".png")), p0,
       width = 10, height = max(6, 0.22 * nrow(summ) + 2.5), dpi = 150)
message("  wrote SUMMARY_fwet_delta_by_site.png")

# ------------------------------------------------------- FIG: all sites f_wet(z)
ph_all$site_lab <- lab_of[ph_all$site]
p_all <- ggplot(ph_all, aes(x = f_wet, y = z)) +
  geom_path(color = "grey55") +
  geom_point(aes(fill = f_wet), shape = 21, size = 2, color = "grey25") +
  geom_vline(aes(xintercept = f_wet_buffer), linetype = 2,
             color = "firebrick", linewidth = 0.4) +
  facet_wrap(~ site_lab, scales = "free", ncol = 6) +
  scale_fill_viridis_c(option = "mako", direction = -1, guide = "none") +
  labs(title = "Wetland fraction seen by the footprint, by measurement height — all sites",
       subtitle = "Rotationally averaged over wind direction. Dashed line = wetland fraction of the whole 5 km buffer.",
       x = "f_wet", y = "measurement height z (m)") +
  theme_minimal(base_size = 8) +
  theme(strip.text = element_text(face = "bold", size = 7))
ggsave(file.path(fp_fig, paste0("SUMMARY_fwet_by_height_all", sfx, ".png")), p_all,
       width = 17, height = max(8, 2.1 * ceiling(nrow(summ) / 6)), dpi = 140)
message("  wrote SUMMARY_fwet_by_height_all.png")

# ------------------------------------------------------- FIGS: exemplar detail
ex <- intersect(exemplars, unique(ph_all$site))
ph_ex  <- ph_all  %>% filter(site %in% ex)
rad_ex <- rad_all %>% filter(site %in% ex) %>% mutate(site_lab = lab_of[site])

if (nrow(ph_ex) > 0) {
  p1 <- ggplot(ph_ex, aes(x = f_wet, y = z)) +
    geom_path(color = "grey55") +
    geom_point(aes(fill = f_wet), shape = 21, size = 3.4, color = "grey20") +
    geom_vline(aes(xintercept = f_wet_buffer), linetype = 2, color = "firebrick") +
    facet_wrap(~ site_lab, scales = "free") +
    scale_fill_viridis_c(option = "mako", direction = -1, guide = "none") +
    labs(title = "Wetland fraction seen by the footprint, by measurement height",
         subtitle = "Rotationally averaged over wind direction. Dashed line = wetland fraction of the whole 5 km buffer.",
         x = "f_wet  (footprint-weighted wetland fraction)",
         y = "measurement height z (m)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fp_fig, paste0("SUMMARY_fwet_by_height", sfx, ".png")), p1,
         width = 11, height = 7, dpi = 150)

  p2 <- ggplot(ph_ex, aes(x = r50, y = z, color = site)) +
    geom_path() + geom_point(size = 2.6) +
    geom_errorbarh(aes(xmin = r50, xmax = r90), height = 0, alpha = 0.45) +
    scale_color_viridis_d(option = "turbo", end = 0.9) +
    labs(title = "Footprint reach grows with measurement height",
         subtitle = "Point = median radius of the 80% footprint; bar extends to the 90th percentile radius",
         x = "distance from tower (m)", y = "measurement height z (m)", color = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(fp_fig, paste0("SUMMARY_footprint_reach_by_height", sfx, ".png")), p2,
         width = 9, height = 6, dpi = 150)

  p3 <- ggplot(rad_ex, aes(x = r, y = g)) +
    geom_area(fill = "#3b6ea5", alpha = 0.32) +
    geom_line(color = "#25496b") +
    geom_vline(data = ph_ex %>% mutate(site_lab = lab_of[site]),
               aes(xintercept = r50, color = z), linewidth = 0.6) +
    facet_wrap(~ site_lab, scales = "free") +
    scale_color_viridis_c(option = "inferno", end = 0.85, name = "sensor\nheight (m)") +
    labs(title = "Where the wetland is, relative to the tower",
         subtitle = "Shaded: fraction of the circle at radius r that is wetland. Vertical lines: median footprint radius per sensor height.",
         x = "distance from tower (m)", y = "wetland fraction of circle  g(r)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(fp_fig, paste0("SUMMARY_radial_wetland_profile", sfx, ".png")), p3,
         width = 11, height = 7, dpi = 150)
  message("  wrote exemplar detail figures")
}

# ------------------------------------------------------- FIG: maps (exemplars)
have_maps <- names(nwi_of)[!vapply(nwi_of, is.null, logical(1))]
if (length(have_maps) > 0) {
  th <- seq(0, 2 * pi, length.out = 181)
  ring_df <- list()
  for (site in have_maps) {
    rs <- ph_all %>% filter(site == !!site)
    for (k in seq_len(nrow(rs))) {
      for (kind in c("median", "90th pct")) {
        rr <- if (kind == "median") rs$r50[k] else rs$r90[k]
        if (!is.finite(rr)) next
        ring_df[[length(ring_df) + 1]] <- tibble(
          site = site, z = rs$z[k], kind = kind,
          grp = paste(site, rs$z[k], kind),
          x = rr * cos(th), y = rr * sin(th))
      }
    }
  }
  ring_df <- bind_rows(ring_df)

  # Site windows differ by an order of magnitude, so each panel is drawn in
  # units of its own half-width. That keeps circles circular under coord_equal,
  # which free facet scales would forbid. True half-width goes in the label.
  map_lab <- vapply(have_maps, function(s)
    sprintf("%s - buffer f_wet %.2f - window +/-%.0f m",
            s, summ$f_wet_buffer[summ$site == s][1], half_of[[s]]), character(1))
  scale_xy <- function(d) {
    h <- unlist(half_of[d$site]); d$x <- d$x / h; d$y <- d$y / h; d
  }
  ring_df <- scale_xy(ring_df); ring_df$site_lab <- map_lab[ring_df$site]

  nwi_df <- list()
  for (site in have_maps) {
    co <- st_coordinates(nwi_of[[site]])
    lev <- intersect(c("L1", "L2", "L3"), colnames(co))
    nwi_df[[site]] <- tibble(
      site = site, x = co[, "X"], y = co[, "Y"],
      grp = do.call(paste, c(list(site), lapply(lev, function(l) co[, l]))))
  }
  nwi_df <- scale_xy(bind_rows(nwi_df)); nwi_df$site_lab <- map_lab[nwi_df$site]

  p4 <- ggplot() +
    geom_polygon(data = nwi_df, aes(x = x, y = y, group = grp),
                 fill = "#4a90c4", color = NA, alpha = 0.75) +
    geom_path(data = ring_df,
              aes(x = x, y = y, group = grp, color = z, linetype = kind),
              linewidth = 0.7) +
    geom_point(data = distinct(ring_df, site_lab), aes(x = 0, y = 0),
               shape = 4, size = 2.4, stroke = 1.1, color = "black") +
    facet_wrap(~ site_lab) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1)) +
    scale_color_viridis_c(option = "inferno", end = 0.85, name = "sensor\nheight (m)") +
    scale_linetype_manual(values = c("median" = "solid", "90th pct" = "22"),
                          name = "footprint\nradius") +
    labs(title = "NWI wetlands around each tower, with per-height footprint rings",
         subtitle = paste("Radius of the 80% footprint at each measurement height, averaged over wind direction. x = tower.",
                          "\nSolid = median radius, dashed = 90th percentile. Each panel is scaled to its own window (see panel title)."),
         x = "distance from tower (fraction of panel window)",
         y = "distance from tower (fraction of panel window)") +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
  ggsave(file.path(fp_fig, paste0("SUMMARY_wetland_maps", sfx, ".png")), p4,
         width = 12, height = 8, dpi = 150)
  message("  wrote SUMMARY_wetland_maps.png")
}

message("\nFigures in: ", fp_fig)
