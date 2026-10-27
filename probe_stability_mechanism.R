# probe_stability_mechanism.R
# Is stability actually driving the flux height trend, or is it a proxy for
# something else (tall towers, site type)?
#
# Four checks:
#
#  1. WITHIN-SITE. The cross-site comparison could be confounded by which sites
#     contribute to which class. Compute the trend for the SAME site in stable
#     and unstable conditions and pair them. If a site flips sign against
#     itself, no between-site confound can explain it.
#
#  2. HEIGHT DEPENDENCE. Is the flip just "tall towers"? Regress the size of the
#     flip on tower height and canopy height.
#
#  3. MULTI-GAS. CO2 and H2O have completely different source distributions and
#     diel cycles from CH4. If all three flip the same way with stability, the
#     cause is shared transport/method, not CH4 biology.
#
#  4. K vs GRADIENT. F_FG = -K(z) * dC/dz. If the true flux were constant with
#     height, F_FG would STILL vary with height whenever K(z) is misspecified,
#     and K errors are stability-dependent. Decompose the height dependence of
#     F into the height dependence of K and of dC/dz.
#
# Usage:  Rscript probe_stability_mechanism.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })

MIN_N <- 30; MIN_LEV <- 3
ck_dir <- file.path(proj_dir, "output", "checkpoints")
sites  <- sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$")))
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE) %>%
  select(site, canopy_height, canopy_class)

wslope <- function(df, ycol) {
  if (nrow(df) < MIN_LEV) return(NA_real_)
  f <- tryCatch(lm(reformulate("z", ycol), data = df, weights = df$n),
                error = function(e) NULL)
  if (is.null(f)) NA_real_ else unname(coef(f)[2])
}

trend <- list(); prof <- list()
for (site in sites) {
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(site, ".rds"))), error = function(e) NULL)
  if (is.null(x)) next
  for (gas in c("CH4", "CO2", "H2O")) {
    f <- x$fluxes[[gas]]
    if (is.null(f) || !all(c("F_FG","z_mid","Stability_500","K","dCdz") %in% names(f))) next
    f <- f %>% filter(is.finite(F_FG), F_outlier %in% c(FALSE, 0), !is.na(z_mid),
                      !is.na(Stability_500))
    if (nrow(f) < 100) next
    scale_F <- median(abs(f$F_FG), na.rm = TRUE)   # for cross-gas comparability
    if (!is.finite(scale_F) || scale_F <= 0) next

    for (st in c("unstable", "neutral", "stable")) {
      lv <- f %>% filter(Stability_500 == st) %>%
        group_by(z = z_mid) %>%
        summarise(n = n(), F = median(F_FG), K = median(K, na.rm = TRUE),
                  dCdz = median(dCdz, na.rm = TRUE), .groups = "drop") %>%
        filter(n >= MIN_N) %>% arrange(z)
      if (nrow(lv) < MIN_LEV) next
      lv$site <- site; lv$gas <- gas; lv$stab <- st; lv$scale_F <- scale_F
      prof[[length(prof) + 1]] <- lv

      trend[[length(trend) + 1]] <- tibble(
        site = site, gas = gas, stab = st,
        slope     = wslope(lv, "F"),
        slope_rel = wslope(lv, "F") / scale_F,   # per metre, per unit typical flux
        slope_K   = wslope(lv, "K") / max(median(lv$K, na.rm = TRUE), 1e-12),
        slope_g   = wslope(lv %>% mutate(absg = abs(dCdz)), "absg") /
                    max(median(abs(lv$dCdz), na.rm = TRUE), 1e-12),
        n_lev = nrow(lv))
    }
  }
}
tr <- bind_rows(trend); pf <- bind_rows(prof)
write.csv(tr, file.path(tab_dir, "stability_mechanism_trends.csv"), row.names = FALSE)

sink(file.path(tab_dir, "stability_mechanism_report.txt"), split = TRUE)
cat("=================================================================\n")
cat(" Is stability really driving the flux height trend?\n")
cat("=================================================================\n")

# ---- 1. within-site pairing ----
cat("\n--- 1. WITHIN-SITE: does the same site flip against itself? ---\n")
for (gas in c("CH4", "CO2", "H2O")) {
  w <- tr %>% filter(gas == !!gas, stab %in% c("stable", "unstable")) %>%
    select(site, stab, slope_rel) %>%
    pivot_wider(names_from = stab, values_from = slope_rel) %>%
    filter(is.finite(stable), is.finite(unstable))
  if (nrow(w) < 5) next
  w$diff <- w$stable - w$unstable
  tt <- t.test(w$diff)
  cat(sprintf("\n  %s  (n = %d sites with both classes)\n", gas, nrow(w)))
  cat(sprintf("    unstable slope < 0 : %2d/%d      stable slope > 0 : %2d/%d\n",
              sum(w$unstable < 0), nrow(w), sum(w$stable > 0), nrow(w)))
  cat(sprintf("    sites flipping sign (unstable<0 AND stable>0): %d/%d (%.0f%%)\n",
              sum(w$unstable < 0 & w$stable > 0), nrow(w),
              100 * mean(w$unstable < 0 & w$stable > 0)))
  cat(sprintf("    paired stable-minus-unstable: mean %+.4g, CI [%+.4g, %+.4g], p = %.3g\n",
              tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
  assign(paste0("w_", gas), w)
}

# ---- 2. is it just tall towers? ----
cat("\n--- 2. Is the flip explained by tower/canopy height? ---\n")
w <- get("w_CH4") %>%
  left_join(pf %>% filter(gas == "CH4") %>% group_by(site) %>%
              summarise(z_top = max(z), .groups = "drop"), by = "site") %>%
  left_join(cm, by = "site")
for (v in c("z_top", "canopy_height")) {
  fit <- tryCatch(lm(reformulate(v, "diff"), data = w), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)
  cat(sprintf("  diff ~ %-14s slope %+.4g (p = %.3g), R^2 = %.3f\n",
              v, coef(s)[2,1], coef(s)[2,4], s$r.squared))
}
cat("  Also: does the flip appear at SHORT towers too?\n")
shortw <- w %>% filter(z_top < median(z_top, na.rm = TRUE))
tallw  <- w %>% filter(z_top >= median(z_top, na.rm = TRUE))
for (nm in c("short", "tall")) {
  v <- if (nm == "short") shortw$diff else tallw$diff
  v <- v[is.finite(v)]
  if (length(v) < 3) next
  tt <- t.test(v)
  cat(sprintf("    %-5s towers (n=%2d): mean flip %+.4g, p = %.3g\n",
              nm, length(v), tt$estimate, tt$p.value))
}

# ---- 3. multi-gas ----
cat("\n--- 3. MULTI-GAS: do CO2 and H2O flip the same way? ---\n")
cat("   (different sources and diel cycles from CH4 — a shared flip implicates\n")
cat("    transport/method rather than CH4 biology)\n\n")
print(as.data.frame(tr %>% group_by(gas, stab) %>%
  summarise(n = n(), mean_rel_slope = signif(mean(slope_rel, na.rm = TRUE), 3),
            frac_negative = round(mean(slope_rel < 0, na.rm = TRUE), 2),
            .groups = "drop")), row.names = FALSE)

# ---- 4. K vs gradient ----
cat("\n--- 4. Does the height dependence come from K or from dC/dz? ---\n")
cat("   Relative per-metre slope of each term, by stability (CH4):\n\n")
print(as.data.frame(tr %>% filter(gas == "CH4") %>% group_by(stab) %>%
  summarise(n = n(),
            F_slope_rel = signif(mean(slope_rel, na.rm = TRUE), 3),
            K_slope_rel = signif(mean(slope_K, na.rm = TRUE), 3),
            absgrad_slope_rel = signif(mean(slope_g, na.rm = TRUE), 3),
            .groups = "drop")), row.names = FALSE)
cat("\n   K rises with height in every class (it must — K grows with z).\n")
cat("   The question is whether F tracks K or tracks the gradient.\n")
sink()

# ---- figures ----
fp <- file.path(fig_dir, "footprint")
wl <- get("w_CH4") %>% mutate(site = reorder(site, diff))
p1 <- ggplot(wl) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_segment(aes(x = site, xend = site, y = unstable, yend = stable),
               color = "grey70") +
  geom_point(aes(x = site, y = unstable, color = "unstable"), size = 2) +
  geom_point(aes(x = site, y = stable,   color = "stable"),   size = 2) +
  coord_flip() +
  scale_color_manual(values = c(unstable = "#e07b39", stable = "#2c5f9e"), name = NULL) +
  labs(title = "Each site against itself: convective vs stable",
       subtitle = "Relative height trend in CH4 flux (per metre, per unit typical flux). A line crossing zero is a within-site sign flip.",
       x = NULL, y = "relative dF/dz") +
  theme_minimal(base_size = 9)
ggsave(file.path(fp, "TEST4_within_site_stability_flip.png"), p1,
       width = 8, height = 9, dpi = 150)

p2 <- ggplot(tr %>% filter(is.finite(slope_rel)),
             aes(x = factor(stab, levels = c("unstable","neutral","stable")),
                 y = slope_rel)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "grey40", width = 0.55) +
  geom_jitter(width = 0.15, size = 1.4, alpha = 0.5) +
  facet_wrap(~ gas, scales = "free_y") +
  labs(title = "The same stability reversal appears in CO2 and H2O",
       subtitle = "Relative height trend in gradient flux, per metre. CO2 and H2O have different sources and diel cycles from CH4.",
       x = NULL, y = "relative dF/dz") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "TEST4_multigas_stability.png"), p2, width = 11, height = 5, dpi = 150)

message("\nReport: ", file.path(tab_dir, "stability_mechanism_report.txt"))
