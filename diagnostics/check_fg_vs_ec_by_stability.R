# check_fg_vs_ec_by_stability.R
# The flux-gradient height trend reverses sign between convective and stable
# conditions. Two readings of that: (a) the method degrades in stable
# conditions, (b) the sources really do differ day vs night. Stability and
# day/night are confounded, so the reversal alone cannot separate them.
#
# The eddy-covariance flux at tower top is an INDEPENDENT measurement of the
# same quantity. If F_FG tracks F_EC in convective conditions but diverges in
# stable ones, that is the method failing, not the surface changing.
#
# Usage:  Rscript check_fg_vs_ec_by_stability.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2) })

GAS <- "CO2"   # NEON provides no EC methane flux; CO2 tests the METHOD, which is gas-independent
ck_dir <- file.path(proj_dir, "output", "checkpoints")
sites <- sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$")))

rows <- list(); pts <- list()
for (site in sites) {
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(site, ".rds"))), error = function(e) NULL)
  if (is.null(x) || is.null(x$fluxes[[GAS]])) next
  f <- x$fluxes[[GAS]]
  if (!all(c("F_FG", "F_EC", "z_mid", "Stability_500") %in% names(f))) next

  # top level only — that is where EC is measured
  ztop <- max(f$z_mid, na.rm = TRUE)
  ff <- f %>% filter(is.finite(F_FG), is.finite(F_EC),
                     F_outlier %in% c(FALSE, 0),
                     abs(z_mid - ztop) < 1e-6, !is.na(Stability_500))
  if (nrow(ff) < 50) next
  ff$site <- site
  pts[[site]] <- ff %>% select(site, Stability_500, F_FG, F_EC)

  for (st in c("unstable", "neutral", "stable")) {
    g <- ff %>% filter(Stability_500 == st)
    if (nrow(g) < 30) next
    fit <- tryCatch(lm(F_FG ~ F_EC, data = g), error = function(e) NULL)
    rows[[length(rows) + 1]] <- tibble(
      site = site, stab = st, n = nrow(g),
      r = suppressWarnings(cor(g$F_FG, g$F_EC)),
      slope = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]),
      med_FG = median(g$F_FG), med_EC = median(g$F_EC),
      bias = median(g$F_FG - g$F_EC))
  }
}
d <- bind_rows(rows)
write.csv(d, file.path(tab_dir, "fg_vs_ec_by_stability.csv"), row.names = FALSE)

sink(file.path(tab_dir, "fg_vs_ec_by_stability.txt"), split = TRUE)
cat("=========================================================\n")
cat(" Does the gradient flux agree with EC, and does that\n")
cat(" agreement depend on stability?  (CO2, top level)\n")
cat("=========================================================\n\n")
print(as.data.frame(d %>% group_by(stab) %>%
  summarise(n_sites = n(), median_r = round(median(r, na.rm = TRUE), 3),
            frac_r_gt_0.3 = round(mean(r > 0.3, na.rm = TRUE), 2),
            median_slope = signif(median(slope, na.rm = TRUE), 3),
            median_bias = signif(median(bias, na.rm = TRUE), 3),
            med_FG = signif(median(med_FG), 3),
            med_EC = signif(median(med_EC), 3), .groups = "drop")),
  row.names = FALSE)

cat("\n  r        = per-site correlation between F_FG and F_EC\n")
cat("  slope    = regression of F_FG on F_EC; 1.0 would be perfect agreement\n")
cat("  bias     = median(F_FG - F_EC)\n")
cat("\nIf r and slope collapse in stable conditions while EC itself stays\n")
cat("well behaved, the stable-condition height trend is a method artifact.\n")

for (st in c("unstable", "neutral", "stable")) {
  v <- d$r[d$stab == st]; v <- v[is.finite(v)]
  if (length(v) < 3) next
  tt <- t.test(v)
  cat(sprintf("\n  %-9s median r = %+.3f, mean %+.3f, CI [%+.3f, %+.3f]",
              st, median(v), tt$estimate, tt$conf.int[1], tt$conf.int[2]))
}
cat("\n")
sink()

pt <- bind_rows(pts)
q <- quantile(c(pt$F_FG, pt$F_EC), c(0.01, 0.99), na.rm = TRUE)
p <- ggplot(pt, aes(x = F_EC, y = F_FG)) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.06, size = 0.5) +
  geom_smooth(method = "lm", color = "#1b6ca8", fill = "#9ecae1", linewidth = 0.7) +
  facet_wrap(~ factor(Stability_500, levels = c("unstable", "neutral", "stable"))) +
  coord_cartesian(xlim = q, ylim = q) +
  labs(title = "Gradient flux vs eddy covariance, by stability (CO2, top level)",
       subtitle = "Dashed red = 1:1. Agreement in convective conditions and collapse in stable ones would mark a method failure.",
       x = "F_EC (umol m-2 s-1)", y = "F_FG (umol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "footprint", "TEST3_FG_vs_EC_by_stability.png"), p,
       width = 11, height = 4.6, dpi = 150)

message("\nReport: ", file.path(tab_dir, "fg_vs_ec_by_stability.txt"))
