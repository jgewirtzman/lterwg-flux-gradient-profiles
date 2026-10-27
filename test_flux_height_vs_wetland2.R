# test_flux_height_vs_wetland2.R
# Corrected sign-contrast test.
#
# Two problems with the first pass:
#   1. dF_topbot (top level minus bottom level) is not normalised by the height
#      range spanned. Tall towers span more height and so mechanically produce a
#      bigger difference. Use the per-metre slope instead.
#   2. The wetland groups are confounded with canopy: the "top sees MORE
#      wetland" sites average 16.7 m canopy vs 8.7 m for "no wetland gradient",
#      and tower top height correlates with canopy height at r = 0.98. So the
#      group contrast partly compares tall forests against short grassland.
#      Control for canopy height and height span explicitly.
#
# Also splits by stability, because the first pass hinted the trend reverses in
# convective conditions — which would point at the flux-gradient method rather
# than at any surface source.

source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })

GAS   <- "CH4"
MIN_N <- 30
MIN_LEV <- 3
NULL_BAND <- 0.01

ck_dir <- file.path(proj_dir, "output", "checkpoints")
fw <- read.csv(file.path(tab_dir, "footprint_fwet_summary.csv"), stringsAsFactors = FALSE)
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE) %>%
  select(site, canopy_height, canopy_class)

sites <- intersect(sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$"))),
                   fw$site)

trend_of <- function(lev) {
  # per-metre slope of median flux on height, weighted by sample size
  if (nrow(lev) < MIN_LEV) return(c(NA_real_, NA_real_, NA_real_))
  fit <- tryCatch(lm(F ~ z, data = lev, weights = n), error = function(e) NULL)
  if (is.null(fit)) return(c(NA_real_, NA_real_, NA_real_))
  s <- summary(fit)
  c(unname(coef(fit)[2]), unname(coef(s)[2, 2]),
    lev$F[which.max(lev$z)] - lev$F[which.min(lev$z)])
}

rows <- list(); stab_rows <- list()
for (site in sites) {
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(site, ".rds"))), error = function(e) NULL)
  if (is.null(x) || is.null(x$fluxes[[GAS]])) next
  f <- x$fluxes[[GAS]] %>%
    filter(is.finite(F_FG), F_outlier %in% c(FALSE, 0), !is.na(z_mid))

  lev <- f %>% group_by(z = z_mid) %>%
    summarise(n = n(), F = median(F_FG), .groups = "drop") %>%
    filter(n >= MIN_N) %>% arrange(z)
  tr <- trend_of(lev)
  if (is.na(tr[1])) next

  rows[[site]] <- tibble(site = site, n_levels = nrow(lev),
    z_bot = min(lev$z), z_top = max(lev$z), z_span = max(lev$z) - min(lev$z),
    F_bot = lev$F[which.min(lev$z)], F_top = lev$F[which.max(lev$z)],
    slope = tr[1], slope_se = tr[2], dF_topbot = tr[3], n_obs = sum(lev$n))

  for (st in c("unstable", "neutral", "stable")) {
    lv <- f %>% filter(Stability_500 == st) %>%
      group_by(z = z_mid) %>%
      summarise(n = n(), F = median(F_FG), .groups = "drop") %>%
      filter(n >= MIN_N) %>% arrange(z)
    t2 <- trend_of(lv)
    if (is.na(t2[1])) next
    stab_rows[[length(stab_rows) + 1]] <- tibble(
      site = site, stab = st, slope = t2[1], dF_topbot = t2[3],
      F_bot = lv$F[which.min(lv$z)], F_top = lv$F[which.max(lv$z)],
      n_obs = sum(lv$n))
  }
}

d <- bind_rows(rows) %>%
  left_join(fw %>% select(site, delta_fwet, f_wet_buffer, status), by = "site") %>%
  left_join(cm, by = "site") %>%
  mutate(wet_group = factor(case_when(
      status %in% c("effectively dry", "saturated wet") ~ "no wetland gradient",
      delta_fwet >  NULL_BAND ~ "top sees MORE wetland",
      delta_fwet < -NULL_BAND ~ "top sees LESS wetland",
      TRUE ~ "no wetland gradient"),
    levels = c("top sees MORE wetland", "no wetland gradient", "top sees LESS wetland")))
sd_ <- bind_rows(stab_rows) %>% left_join(d %>% select(site, wet_group), by = "site")

write.csv(d,   file.path(tab_dir, "flux_height_vs_wetland_v2.csv"), row.names = FALSE)
write.csv(sd_, file.path(tab_dir, "flux_height_by_stability.csv"),  row.names = FALSE)

sink(file.path(tab_dir, "flux_height_vs_wetland_report_v2.txt"), split = TRUE)
cat("==============================================================\n")
cat(" CH4 flux height trend vs footprint wetland - corrected test\n")
cat("==============================================================\n")
cat("\nPrimary metric: slope of median F_FG on height (umol m-2 s-1 per m),\n")
cat("which unlike the raw top-minus-bottom difference does not scale with\n")
cat("how much height the tower spans.\n\n")

cat("--- 0. The phenomenon ---\n")
cat(sprintf("  sites: %d | uptake at bottom (F_bot<0): %d | slope>0: %d\n",
            nrow(d), sum(d$F_bot < 0), sum(d$slope > 0)))
tt <- t.test(d$slope)
cat(sprintf("  mean slope %+.3g  95%% CI [%+.3g, %+.3g]  p=%.3g\n",
            tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))

cat("\n--- 1. Slope by wetland group ---\n")
for (g in levels(d$wet_group)) {
  v <- d$slope[d$wet_group == g]; v <- v[is.finite(v)]
  if (length(v) < 2) next
  t1 <- t.test(v)
  cat(sprintf("  %-24s n=%2d  mean %+.3g  CI [%+.3g, %+.3g]  p=%.3g  %s\n",
              g, length(v), t1$estimate, t1$conf.int[1], t1$conf.int[2], t1$p.value,
              ifelse(t1$p.value < 0.05, "SIGNIFICANT", "n.s.")))
}
cat("\n  Wetland-footprint predicts: first > 0, second = 0, third < 0.\n")

cat("\n--- 2. Confounding: what else differs between these groups? ---\n")
print(as.data.frame(d %>% group_by(wet_group) %>%
  summarise(n = n(), canopy_h = round(mean(canopy_height, na.rm = TRUE), 1),
            z_top = round(mean(z_top), 1), z_span = round(mean(z_span), 1),
            .groups = "drop")), row.names = FALSE)

cat("\n--- 3. Regressions on slope ---\n")
show <- function(lbl, fm) {
  fit <- tryCatch(lm(fm, data = d), error = function(e) NULL)
  if (is.null(fit)) return(invisible())
  cat("\n ", lbl, "\n"); print(round(summary(fit)$coefficients, 6))
  cat(sprintf("  adj R^2 = %.3f, n = %d\n", summary(fit)$adj.r.squared, nobs(fit)))
}
show("slope ~ delta_fwet",                    slope ~ delta_fwet)
show("slope ~ canopy_height",                 slope ~ canopy_height)
show("slope ~ delta_fwet + canopy_height",    slope ~ delta_fwet + canopy_height)
show("slope ~ delta_fwet + canopy_height + z_span", slope ~ delta_fwet + canopy_height + z_span)

cat("\n--- 4. Stability split (the method check) ---\n")
print(as.data.frame(sd_ %>% group_by(stab) %>%
  summarise(n_sites = n(), mean_slope = signif(mean(slope), 3),
            median_slope = signif(median(slope), 3),
            frac_pos = round(mean(slope > 0), 2),
            mean_Fbot = signif(mean(F_bot), 3),
            mean_Ftop = signif(mean(F_top), 3), .groups = "drop")), row.names = FALSE)
for (st in c("unstable", "neutral", "stable")) {
  v <- sd_$slope[sd_$stab == st]; v <- v[is.finite(v)]
  if (length(v) < 3) next
  t2 <- t.test(v)
  cat(sprintf("  %-9s n=%2d mean %+.3g CI [%+.3g,%+.3g] p=%.3g\n",
              st, length(v), t2$estimate, t2$conf.int[1], t2$conf.int[2], t2$p.value))
}
cat("\n  If the height trend only appears in stable conditions, that points at\n")
cat("  the flux-gradient method (roughness sublayer, decoupling, storage)\n")
cat("  rather than at any surface or canopy source.\n")

cat("\n--- 5. Stability split within wetland groups ---\n")
print(as.data.frame(sd_ %>% group_by(stab, wet_group) %>%
  summarise(n = n(), mean_slope = signif(mean(slope), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = wet_group, values_from = mean_slope)), row.names = FALSE)
sink()

# ---- figures ---------------------------------------------------------------
fp_fig <- file.path(fig_dir, "footprint")
pal <- c("top sees MORE wetland" = "#1b6ca8", "no wetland gradient" = "#8a8a8a",
         "top sees LESS wetland" = "#c0392b")

p1 <- ggplot(d, aes(x = delta_fwet, y = slope, color = wet_group)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey85",
              linewidth = 0.5, inherit.aes = FALSE, aes(x = delta_fwet, y = slope)) +
  geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = site), size = 2.4, max.overlaps = 20,
                           show.legend = FALSE) +
  scale_color_manual(values = pal, name = NULL) +
  labs(title = "CH4 flux height trend vs footprint wetland gradient",
       subtitle = "Slope of median F_FG on height, per metre. Wetland-footprint predicts a rising line through the origin.",
       x = "delta f_wet (top - bottom sensor)",
       y = "dF_FG/dz  (umol m-2 s-1 per m)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp_fig, "TEST2_slope_vs_dfwet.png"), p1, width = 10, height = 7, dpi = 150)

p2 <- ggplot(sd_, aes(x = stab, y = slope, color = wet_group)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_boxplot(aes(group = stab), outlier.shape = NA, width = 0.55,
               fill = NA, color = "grey40") +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.85) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_discrete(limits = c("unstable", "neutral", "stable")) +
  labs(title = "The height trend depends on stability",
       subtitle = "Slope of median F_FG on height, computed separately within each stability class",
       x = NULL, y = "dF_FG/dz  (umol m-2 s-1 per m)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp_fig, "TEST2_slope_by_stability.png"), p2, width = 10, height = 6, dpi = 150)

message("\nReport: ", file.path(tab_dir, "flux_height_vs_wetland_report_v2.txt"))
