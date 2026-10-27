# test_flux_height_vs_wetland.R
# Sign-contrast test: is the height trend in CH4 flux explained by footprint
# wetland composition, or by something else (canopy/stem source, or method)?
#
# The logic. We know the sign of df_wet/dz at every site. The two hypotheses
# predict different things, site by site:
#
#   wetland-footprint    dF/dz should track delta_fwet: positive where the top
#                        sensor sees more wetland, NEGATIVE where it sees less,
#                        and ZERO where there is no wetland gradient at all.
#   canopy/stem source   dF/dz > 0 regardless of delta_fwet, scaling instead
#                        with canopy structure.
#
# So the decisive groups are the sites where delta_fwet is <= 0. If the flux
# still rises with height there, footprint wetness cannot be the cause.
#
# This uses only the SIGN of a large, stable, site-level property, so it is
# robust to the things we cannot currently pin down: no wind direction, the
# 5 km NWI buffer truncation, and absolute f_wet being biased.
#
# Usage:  Rscript test_flux_height_vs_wetland.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })

GAS       <- "CH4"
MIN_N     <- 30      # min observations for a site x level median to count
MIN_LEV   <- 3       # min levels for a site to give a height trend
NULL_BAND <- 0.01    # |delta_fwet| below this = no usable wetland gradient

ck_dir <- file.path(proj_dir, "output", "checkpoints")

# ---- wetland gradient per site -------------------------------------------
fw <- read.csv(file.path(tab_dir, "footprint_fwet_summary.csv"),
               stringsAsFactors = FALSE)

# ---- canopy structure (mechanism-1 predictor) ----------------------------
cm_path <- file.path(tab_dir, "canopy_metrics.csv")
cm <- if (file.exists(cm_path)) read.csv(cm_path, stringsAsFactors = FALSE) else NULL
if (!is.null(cm)) {
  hcol <- intersect(c("canopy_height", "canopy_height_m", "height_canopy"), names(cm))
  ccol <- intersect(c("canopy_class", "structure_class"), names(cm))
  cm <- cm %>% select(site = 1,
                      canopy_height = all_of(if (length(hcol)) hcol[1] else names(cm)[2]),
                      any_of(setNames(ccol[1], "canopy_class")))
}

# ---- per-site flux height trend ------------------------------------------
sites <- sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$")))
sites <- intersect(sites, fw$site)
message("Sites with both flux checkpoints and wetland gradient: ", length(sites))

lev_rows <- list(); site_rows <- list(); layer_rows <- list()

for (site in sites) {
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(site, ".rds"))),
                error = function(e) NULL)
  if (is.null(x) || is.null(x$fluxes[[GAS]])) next
  f <- x$fluxes[[GAS]]
  if (!all(c("F_FG", "z_mid", "F_outlier") %in% names(f))) next

  f <- f %>%
    filter(!is.na(F_FG), is.finite(F_FG), !isTRUE(F_outlier),
           F_outlier %in% c(FALSE, 0), !is.na(z_mid))

  # Per level: median flux (CH4 fluxes are small and heavy-tailed, so median
  # rather than mean), plus a stable/unstable split as a sensitivity check.
  lev <- f %>%
    group_by(z = z_mid) %>%
    summarise(n = n(), F_med = median(F_FG), F_mean = mean(F_FG),
              F_se = sd(F_FG) / sqrt(n()),
              F_med_unstable = median(F_FG[which(Stability_500 == "unstable")]),
              n_unstable = sum(Stability_500 == "unstable", na.rm = TRUE),
              .groups = "drop") %>%
    filter(n >= MIN_N) %>% arrange(z)
  if (nrow(lev) < MIN_LEV) next
  lev$site <- site
  lev_rows[[site]] <- lev

  # Height trend: bottom-to-top difference, and a slope across levels.
  slope <- tryCatch(unname(coef(lm(F_med ~ z, data = lev, weights = n))[2]),
                    error = function(e) NA_real_)
  slope_uns <- tryCatch(unname(coef(lm(F_med_unstable ~ z,
                                       data = lev %>% filter(is.finite(F_med_unstable)),
                                       weights = n_unstable))[2]),
                        error = function(e) NA_real_)
  site_rows[[site]] <- tibble(
    site = site, n_levels = nrow(lev),
    z_bot_f = first(lev$z), z_top_f = last(lev$z),
    F_bot = first(lev$F_med), F_top = last(lev$F_med),
    dF_topbot = last(lev$F_med) - first(lev$F_med),
    dF_slope = slope, dF_slope_unstable = slope_uns,
    n_obs = sum(lev$n))

  # Layer-resolved S_net (dF + storage) as the physically cleaner quantity.
  ss <- x$source_sink[[GAS]]
  if (!is.null(ss) && nrow(ss) > 0 && "S_net" %in% names(ss)) {
    layer_rows[[site]] <- ss %>%
      filter(is.finite(S_net)) %>%
      group_by(z_lower, z_upper) %>%
      summarise(n = n(), S_net_med = median(S_net),
                dF_med = median(dF, na.rm = TRUE), .groups = "drop") %>%
      filter(n >= MIN_N) %>% mutate(site = site)
  }
}

lev_all  <- bind_rows(lev_rows)
site_all <- bind_rows(site_rows)
if (nrow(site_all) == 0) stop("no sites produced a flux height trend")

# ---- join wetland gradient + canopy ---------------------------------------
d <- site_all %>%
  left_join(fw %>% select(site, delta_fwet, fwet_bot, fwet_top,
                          f_wet_buffer, status, beyond_buffer), by = "site")
if (!is.null(cm)) d <- d %>% left_join(cm, by = "site")

d <- d %>%
  mutate(wet_group = case_when(
    status %in% c("effectively dry", "saturated wet") ~ "no wetland gradient",
    delta_fwet >  NULL_BAND                           ~ "top sees MORE wetland",
    delta_fwet < -NULL_BAND                           ~ "top sees LESS wetland",
    TRUE                                              ~ "no wetland gradient"),
    wet_group = factor(wet_group,
      levels = c("top sees MORE wetland", "no wetland gradient",
                 "top sees LESS wetland")))

write.csv(d, file.path(tab_dir, "flux_height_vs_wetland.csv"), row.names = FALSE)
write.csv(lev_all, file.path(tab_dir, "flux_by_level_by_site.csv"), row.names = FALSE)

# ---- the test --------------------------------------------------------------
sink(file.path(tab_dir, "flux_height_vs_wetland_report.txt"), split = TRUE)

cat("=====================================================================\n")
cat(" Does the CH4 flux height trend follow footprint wetland composition?\n")
cat("=====================================================================\n\n")
cat("Sites:", nrow(d), " | flux units: umol m-2 s-1 (F_FG)\n")
cat("dF_topbot = median F_FG at highest level - median F_FG at lowest level\n\n")

cat("--- 0. Is the phenomenon there at all? ---\n")
cat(sprintf("  sites with F_bot < 0 (uptake at the bottom): %d / %d\n",
            sum(d$F_bot < 0, na.rm = TRUE), nrow(d)))
cat(sprintf("  sites with dF_topbot > 0 (flux rises with height): %d / %d\n",
            sum(d$dF_topbot > 0, na.rm = TRUE), nrow(d)))
bt <- tryCatch(t.test(d$dF_topbot), error = function(e) NULL)
if (!is.null(bt)) cat(sprintf("  mean dF_topbot = %.4g  (95%% CI %.4g to %.4g, p = %.3g)\n",
                              bt$estimate, bt$conf.int[1], bt$conf.int[2], bt$p.value))

cat("\n--- 1. Height trend by wetland group ---\n")
grp <- d %>% group_by(wet_group) %>%
  summarise(n_sites = n(),
            mean_dF = mean(dF_topbot, na.rm = TRUE),
            median_dF = median(dF_topbot, na.rm = TRUE),
            frac_positive = mean(dF_topbot > 0, na.rm = TRUE),
            .groups = "drop")
print(as.data.frame(grp), row.names = FALSE)

cat("\n  The decisive rows are the last two. If the wetland-footprint\n")
cat("  explanation were driving this, 'no wetland gradient' should sit at\n")
cat("  zero and 'top sees LESS wetland' should be NEGATIVE.\n\n")

for (g in levels(d$wet_group)) {
  v <- d$dF_topbot[d$wet_group == g]
  v <- v[is.finite(v)]
  if (length(v) < 2) { cat(sprintf("  %-24s n=%d (too few to test)\n", g, length(v))); next }
  tt <- t.test(v)
  cat(sprintf("  %-24s n=%2d  mean %+.4g  95%% CI [%+.4g, %+.4g]  p=%.3g  %s\n",
              g, length(v), tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value,
              ifelse(tt$p.value < 0.05,
                     ifelse(tt$estimate > 0, "significantly POSITIVE", "significantly NEGATIVE"),
                     "not distinguishable from 0")))
}

cat("\n--- 2. Does dF track delta_fwet across sites? ---\n")
fit1 <- tryCatch(lm(dF_topbot ~ delta_fwet, data = d), error = function(e) NULL)
if (!is.null(fit1)) {
  s1 <- summary(fit1)
  cat(sprintf("  dF_topbot ~ delta_fwet : slope %+.4g (SE %.3g), p = %.3g, R^2 = %.3f, n = %d\n",
              coef(s1)[2, 1], coef(s1)[2, 2], coef(s1)[2, 4], s1$r.squared,
              nobs(fit1)))
  cat("  (wetland-footprint predicts a clearly POSITIVE slope here)\n")
}
dd <- d %>% filter(status == "can discriminate")
fit1b <- tryCatch(lm(dF_topbot ~ delta_fwet, data = dd), error = function(e) NULL)
if (!is.null(fit1b)) {
  s1b <- summary(fit1b)
  cat(sprintf("  same, discriminating sites only: slope %+.4g, p = %.3g, R^2 = %.3f, n = %d\n",
              coef(s1b)[2, 1], coef(s1b)[2, 4], s1b$r.squared, nobs(fit1b)))
}

cat("\n--- 3. Canopy structure (mechanism 1) ---\n")
if ("canopy_height" %in% names(d) && sum(is.finite(d$canopy_height)) > 5) {
  fit2 <- tryCatch(lm(dF_topbot ~ canopy_height, data = d), error = function(e) NULL)
  if (!is.null(fit2)) {
    s2 <- summary(fit2)
    cat(sprintf("  dF_topbot ~ canopy_height : slope %+.4g, p = %.3g, R^2 = %.3f\n",
                coef(s2)[2, 1], coef(s2)[2, 4], s2$r.squared))
  }
  fit3 <- tryCatch(lm(dF_topbot ~ delta_fwet + canopy_height, data = d),
                   error = function(e) NULL)
  if (!is.null(fit3)) {
    s3 <- summary(fit3)
    cat("  both together:\n"); print(round(coef(s3), 5))
  }
} else cat("  canopy_height unavailable\n")

cat("\n--- 4. Sensitivity: unstable (convective) periods only ---\n")
du <- d %>% filter(is.finite(dF_slope_unstable))
if (nrow(du) > 3) {
  cat(sprintf("  sites: %d, mean slope %+.4g\n", nrow(du), mean(du$dF_slope_unstable)))
  g2 <- du %>% group_by(wet_group) %>%
    summarise(n = n(), mean_slope = mean(dF_slope_unstable), .groups = "drop")
  print(as.data.frame(g2), row.names = FALSE)
}

cat("\n--- 5. Caveats carried into this table ---\n")
cat(sprintf("  sites whose footprints exceed the 5 km NWI buffer: %d\n",
            sum(d$beyond_buffer, na.rm = TRUE)))
cat("  delta_fwet is a symmetric-wind estimate (no wind direction).\n")
cat("  F_FG at the lowest level may sit in the roughness sublayer where\n")
cat("  K-theory is unreliable — a method artifact would also produce a\n")
cat("  height trend, and this test does not rule that out.\n")
sink()

# ---- figures ---------------------------------------------------------------
fp_fig <- file.path(fig_dir, "footprint")
pal <- c("top sees MORE wetland" = "#1b6ca8",
         "no wetland gradient"   = "#8a8a8a",
         "top sees LESS wetland" = "#c0392b")

p1 <- ggplot(d, aes(x = delta_fwet, y = dF_topbot, color = wet_group)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey50", linetype = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "grey30",
              fill = "grey80", linewidth = 0.5, inherit.aes = FALSE,
              aes(x = delta_fwet, y = dF_topbot)) +
  geom_point(size = 2.6) +
  ggrepel::geom_text_repel(aes(label = site), size = 2.5, max.overlaps = 12,
                           show.legend = FALSE) +
  scale_color_manual(values = pal, name = NULL) +
  labs(title = "Does the CH4 flux height trend follow footprint wetland composition?",
       subtitle = paste("If wetland composition drove it, points would lie on a rising line through the origin.",
                        "\nRed points (top sensor sees LESS wetland) should sit BELOW zero."),
       x = "delta f_wet   (top - bottom sensor)",
       y = "delta F_FG   (top - bottom level, umol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp_fig, "TEST_dF_vs_dfwet.png"), p1, width = 10, height = 7, dpi = 150)

p2 <- ggplot(d, aes(x = wet_group, y = dF_topbot, color = wet_group)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = 2) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  geom_jitter(width = 0.13, height = 0, size = 2.2, alpha = 0.85) +
  scale_color_manual(values = pal, guide = "none") +
  labs(title = "Flux height trend, grouped by what the footprint wetness does with height",
       subtitle = "Wetland-footprint explanation predicts: left > 0, middle = 0, right < 0",
       x = NULL, y = "delta F_FG  (top - bottom level, umol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp_fig, "TEST_dF_by_wetland_group.png"), p2, width = 9, height = 6, dpi = 150)

lev_all2 <- lev_all %>% left_join(d %>% select(site, wet_group), by = "site")
p3 <- ggplot(lev_all2, aes(x = F_med, y = z, group = site, color = wet_group)) +
  geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  geom_path(alpha = 0.75) + geom_point(size = 1.3) +
  facet_wrap(~ wet_group) +
  scale_color_manual(values = pal, guide = "none") +
  labs(title = "CH4 flux profiles by site, grouped by the wetland gradient",
       subtitle = "Median F_FG at each level. Negative = uptake.",
       x = "median F_FG (umol m-2 s-1)", y = "height z (m)") +
  theme_minimal(base_size = 11) +
  coord_cartesian(xlim = quantile(lev_all2$F_med, c(0.02, 0.98), na.rm = TRUE))
ggsave(file.path(fp_fig, "TEST_flux_profiles_by_group.png"), p3, width = 12, height = 6, dpi = 150)

message("\nWrote figures + tables. Report: ",
        file.path(tab_dir, "flux_height_vs_wetland_report.txt"))
