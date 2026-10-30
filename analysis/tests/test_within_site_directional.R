# test_within_site_directional.R
#
# The within-site mechanism test, which only becomes possible with directional
# footprints.
#
# Design: for each adjacent sensor pair at each 30-min timestamp, does the
# concentration difference between the two heights respond to the difference in
# wetland fraction between their two footprints?
#
#   dConc ~ delta_fwet + delta_Area   (+ met controls), site fixed effects
#
# Prediction if footprint composition drives the apparent height structure: the
# upper sensor seeing more wetland raises its concentration relative to the
# lower one, so the delta_fwet slope should be POSITIVE.
#
# Why this works now and did not before:
#   * dConc is measured per pair per timestamp, so it has none of the QC-overlap
#     problem that makes top-vs-bottom FLUX pairing impossible (median 19
#     matched timestamps per site).
#   * Under symmetric-wind averaging delta_fwet was ~80% footprint size, i.e.
#     nearly collinear with delta_Area, so the two could not be separated.
#     Directional f_wet cuts that to ~30%.
#
# The symmetric version of the same regression is run alongside as a control:
# if directional finds signal where symmetric cannot, that is the identification
# working rather than a new artefact.
#
# Usage:  Rscript test_within_site_directional.R

source("config.R")
suppressPackageStartupMessages({ library(lubridate); library(ggplot2) })

dir_ck <- file.path(proj_dir, "output", "footprint", "directional")
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE)

files <- list.files(dir_ck, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no directional checkpoints")

pairs <- list()
for (fp in files) {
  d <- readRDS(fp)
  site <- d$site[1]
  if (!"time_round" %in% names(d)) {
    message(site, ": checkpoint predates time_round — re-run analysis/footprint/compute_directional_fwet.R")
    next
  }
  sd_ <- tryCatch(load_all_sites(site, data_dir), error = function(e) NULL)
  if (is.null(sd_) || is.null(sd_[[site]]$aligned$CH4)) next
  al <- sd_[[site]]$aligned$CH4 %>%
    mutate(time_round = round_date(timeMid, "30 minutes"),
           pos_A = as.numeric(TowerPosition_A),
           pos_B = as.numeric(TowerPosition_B)) %>%
    filter(is.finite(as.numeric(dConc))) %>%
    select(time_round, pos_A, pos_B, dConc, ustar = ustar_interp, zoL,
           dLevelsAminusB) %>%
    distinct(time_round, pos_A, pos_B, .keep_all = TRUE)

  # adjacent pairs from the per-height footprint table
  w <- d %>% select(time_round, TowerPosition, f_wet_dir, f_wet_sym, area)
  jp <- w %>% rename(pos_A = TowerPosition, fwA_dir = f_wet_dir,
                     fwA_sym = f_wet_sym, areaA = area) %>%
    inner_join(w %>% rename(pos_B = TowerPosition, fwB_dir = f_wet_dir,
                            fwB_sym = f_wet_sym, areaB = area),
               by = "time_round", relationship = "many-to-many") %>%
    filter(pos_A == pos_B + 1) %>%          # adjacent, A above B
    mutate(dfwet_dir = fwA_dir - fwB_dir,
           dfwet_sym = fwA_sym - fwB_sym,
           dArea     = areaA - areaB)
  if (!nrow(jp)) next
  jj <- jp %>% inner_join(al, by = c("time_round", "pos_A", "pos_B"))
  if (nrow(jj) < 20) next
  jj$site <- site
  pairs[[site]] <- jj
}
p <- bind_rows(pairs)
if (!nrow(p)) stop("no joined pairs")
p <- p %>% left_join(cm %>% select(site, canopy_height), by = "site") %>%
  mutate(dConc = as.numeric(dConc), ustar = as.numeric(ustar), zoL = as.numeric(zoL))
write.csv(p, file.path(tab_dir, "within_site_directional_pairs.csv"), row.names = FALSE)

sink(file.path(tab_dir, "within_site_directional_report.txt"), split = TRUE)
cat("==================================================================\n")
cat(" WITHIN-SITE test: does dConc respond to the wetland difference\n")
cat(" between the two sensors' footprints?\n")
cat("==================================================================\n")
cat(sprintf("\n  pairs: %d | sites: %d\n", nrow(p), length(unique(p$site))))

cat("\n--- how collinear is delta_fwet with delta_Area? ---\n")
cat(sprintf("  symmetric  : r = %.3f\n", cor(p$dfwet_sym, p$dArea, use = "complete.obs")))
cat(sprintf("  directional: r = %.3f\n", cor(p$dfwet_dir, p$dArea, use = "complete.obs")))
cat("  (the whole point: directional should be far less collinear)\n")

shw <- function(lbl, fm, dat = p) {
  fit <- tryCatch(lm(fm, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(invisible())
  s <- summary(fit)
  keep <- setdiff(rownames(s$coefficients), grep("^factor\\(site\\)", rownames(s$coefficients), value = TRUE))
  cat("\n  ", lbl, "  (n = ", nobs(fit), ")\n", sep = "")
  print(round(s$coefficients[keep, , drop = FALSE], 6))
  cat(sprintf("   adj R^2 = %.4f\n", s$adj.r.squared))
}

cat("\n--- pooled, site fixed effects ---\n")
shw("dConc ~ dfwet_dir + factor(site)",                  dConc ~ dfwet_dir + factor(site))
shw("dConc ~ dfwet_dir + dArea + factor(site)",          dConc ~ dfwet_dir + dArea + factor(site))
shw("dConc ~ dfwet_dir + dArea + ustar + zoL + factor(site)",
    dConc ~ dfwet_dir + dArea + ustar + zoL + factor(site))
cat("\n  CONTROL — the same model with the SYMMETRIC (non-identifying) f_wet:\n")
shw("dConc ~ dfwet_sym + dArea + factor(site)",          dConc ~ dfwet_sym + dArea + factor(site))

cat("\n--- per-site slopes on dfwet_dir (with dArea controlled) ---\n")
ps <- list()
for (s in unique(p$site)) {
  d <- p %>% filter(site == s)
  if (nrow(d) < 30) next
  f <- tryCatch(lm(dConc ~ dfwet_dir + dArea, data = d), error = function(e) NULL)
  if (is.null(f) || !"dfwet_dir" %in% rownames(summary(f)$coefficients)) next
  cf <- summary(f)$coefficients["dfwet_dir", ]
  ps[[length(ps) + 1]] <- tibble(site = s, n = nrow(d), slope = cf[1],
                                 se = cf[2], p = cf[4],
                                 canopy_height = d$canopy_height[1])
}
ps <- bind_rows(ps)
if (nrow(ps)) {
  print(as.data.frame(ps %>% arrange(p) %>%
    mutate(slope = signif(slope, 3), se = signif(se, 3), p = signif(p, 3))),
    row.names = FALSE)
  cat(sprintf("\n  sites with positive slope: %d/%d | significant at p<0.05: %d\n",
              sum(ps$slope > 0), nrow(ps), sum(ps$p < 0.05)))
  tt <- t.test(ps$slope)
  cat(sprintf("  mean per-site slope %+.4g, 95%% CI [%+.4g, %+.4g], p = %.4f\n",
              tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
}
cat("\n--- caveats ---\n")
cat("  Wetland is static in NWI, so delta_fwet varies only through the footprint;\n")
cat("  a directional signal is compositional but still footprint-model dependent.\n")
cat("  Rotation convention assumes angZaxsErth is meteorological (direction FROM).\n")
sink()

if (nrow(ps)) {
  ps <- ps %>% mutate(site = reorder(site, slope))
  g <- ggplot(ps, aes(x = slope, y = site)) +
    geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
    geom_errorbarh(aes(xmin = slope - se, xmax = slope + se), height = 0, color = "grey60") +
    geom_point(aes(color = p < 0.05), size = 2.6) +
    scale_color_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#8a8a8a"),
                       labels = c("n.s.", "p < 0.05"), name = NULL) +
    labs(title = "Within-site response of dConc to the footprint wetland difference",
         subtitle = "Directional f_wet, delta_Area controlled. Positive = upper sensor seeing more wetland raises its concentration.",
         x = "slope of dConc on delta f_wet", y = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(fig_dir, "footprint", "WITHIN_site_dfwet_slopes.png"), g,
         width = 9, height = 6.5, dpi = 150)
}
message("\nReport: ", file.path(tab_dir, "within_site_directional_report.txt"))
