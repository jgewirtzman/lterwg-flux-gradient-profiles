# test_vertical_change_vs_wetland.R
#
# THE question: CH4 flux is negative at the lowest level and becomes less
# negative or positive higher up. Is that vertical CHANGE caused by the upper
# footprint containing more wetland, or by a source at height (canopy/stems)?
#
# So the dependent variable is the top-minus-bottom flux difference, NOT the
# absolute flux at the top. Regressing absolute top-level flux on wetland only
# shows that wetter landscapes emit more, which does not address the vertical
# structure at all.
#
# dF_column = sum of per-layer dF across adjacent layers, within a timestamp
#           = F(top level) - F(bottom level)
# Using adjacent layers rather than pairing the two extreme levels directly
# keeps thousands of timestamps per site instead of a handful, because adjacent
# levels co-occur far more often than the extremes do.
#
# Predictors:
#   delta_fwet     - wetland fraction seen by top sensor minus bottom sensor.
#                    Mechanism "footprint wetland" predicts a POSITIVE slope.
#   canopy_height  - mechanism "source at height" predicts a POSITIVE slope.
#
# Basis: PAR day/night, 10 km NWI buffer, existing flux QC.
#
# Usage:  Rscript test_vertical_change_vs_wetland.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(lubridate) })

NWI_KM <- 10
ck_dir <- file.path(proj_dir, "output", "checkpoints")
fw <- read.csv(file.path(tab_dir, sprintf("footprint_fwet_summary_%dkm.csv", NWI_KM)),
               stringsAsFactors = FALSE)
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE)
sites <- intersect(sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$"))),
                   fw$site)

rows <- list()
for (s in sites) {
  chm <- cm$canopy_height[cm$site == s]
  if (!length(chm) || is.na(chm)) next
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(s, ".rds"))), error = function(e) NULL)
  if (is.null(x) || is.null(x$source_sink$CH4) || is.null(x$met)) next
  if (!"PAR" %in% names(x$met)) next
  par <- x$met %>% select(t = timeMid, PAR) %>% filter(!is.na(PAR))

  ss <- x$source_sink$CH4 %>%
    filter(is.finite(dF)) %>%
    mutate(Canopy_L1 = paste0(ifelse(z_upper > chm, "A", "W"),
                              ifelse(z_lower > chm, "A", "W"))) %>%
    inner_join(par, by = c("time_round" = "t")) %>%
    mutate(tod = ifelse(PAR > 10, "Day", "Night"))
  if (nrow(ss) < 200) next

  for (excl in c(FALSE, TRUE)) {
    a <- if (excl) ss %>% filter(Canopy_L1 != "WW") else ss
    if (nrow(a) < 100) next
    for (td in c("Day", "Night")) {
      b <- a %>% filter(tod == td)
      if (nrow(b) < 100) next
      col <- b %>% group_by(time_round) %>%
        summarise(dF_col = sum(dF), n_lay = n(), .groups = "drop")
      if (nrow(col) < 50) next
      rows[[length(rows) + 1]] <- tibble(
        site = s, tod = td, excl_WW = excl,
        dF_col = median(col$dF_col), n_ts = nrow(col),
        z_lo = min(b$z_lower), z_hi = max(b$z_upper),
        canopy_height = chm)
    }
  }
}
d <- bind_rows(rows) %>%
  left_join(fw %>% select(site, delta_fwet, fwet_top, fwet_bot, f_wet_buffer,
                          status, beyond_buffer), by = "site") %>%
  mutate(z_span = z_hi - z_lo,
         dF_per_m = dF_col / z_span)
write.csv(d, file.path(tab_dir, "vertical_change_vs_wetland.csv"), row.names = FALSE)

sink(file.path(tab_dir, "vertical_change_report.txt"), split = TRUE)
cat("==================================================================\n")
cat(" Is the BOTTOM-TO-TOP change in CH4 flux explained by wetland,\n")
cat(" or by a source at height?\n")
cat("==================================================================\n")
cat("\n dF_col = F(top level) - F(bottom level), nmol m-2 s-1\n")
cat(" Mechanism predictions: wetland -> positive slope on delta_fwet\n")
cat("                        source at height -> positive slope on canopy_height\n")

shw <- function(lbl, fm, dat) {
  dat <- dat %>% filter(is.finite(dF_col), is.finite(delta_fwet))
  if (nrow(dat) < 10) { cat("\n  ", lbl, ": too few sites\n", sep = ""); return(invisible()) }
  fit <- tryCatch(lm(fm, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(invisible())
  s <- summary(fit)
  cat("\n  ", lbl, "  (n = ", nobs(fit), ")\n", sep = "")
  print(round(s$coefficients, 6))
  cat(sprintf("   adj R^2 = %.3f\n", s$adj.r.squared))
}

for (ex in c(FALSE, TRUE)) {
  cat("\n\n######## ", if (ex) "EXCLUDING within-canopy (WW) layers"
                       else "ALL layers", " ########\n", sep = "")
  for (td in c("Night", "Day")) {
    sub <- d %>% filter(excl_WW == ex, tod == td)
    cat("\n--- ", td, " (", nrow(sub), " sites) ---\n", sep = "")
    cat(sprintf("  median dF_col = %+.4g   | positive at %d/%d sites\n",
                median(sub$dF_col, na.rm = TRUE),
                sum(sub$dF_col > 0, na.rm = TRUE), nrow(sub)))
    tt <- tryCatch(t.test(sub$dF_col), error = function(e) NULL)
    if (!is.null(tt)) cat(sprintf("  mean %+.4g, 95%% CI [%+.4g, %+.4g], p = %.3g\n",
                                  tt$estimate, tt$conf.int[1], tt$conf.int[2], tt$p.value))
    shw("dF_col ~ delta_fwet",                  dF_col ~ delta_fwet, sub)
    shw("dF_col ~ canopy_height",               dF_col ~ canopy_height, sub)
    shw("dF_col ~ delta_fwet + canopy_height",  dF_col ~ delta_fwet + canopy_height, sub)
    shw("dF_per_m ~ delta_fwet + canopy_height", dF_per_m ~ delta_fwet + canopy_height, sub)
  }
}
cat("\n--- caveats ---\n")
cat("  delta_fwet is a symmetric-wind estimate. Unlike the absolute-flux\n")
cat("  regression, this one CANNOT fall back on buffer-scale wetness as a\n")
cat("  footprint-free control: delta_fwet is intrinsically a footprint\n")
cat("  quantity, and under rotational averaging it is largely driven by the\n")
cat("  footprint SIZE difference between the two heights.\n")
sink()

pd <- d %>% filter(tod == "Night", is.finite(dF_col), is.finite(delta_fwet)) %>%
  mutate(basis = ifelse(excl_WW, "excluding WW layers", "all layers"))
p <- ggplot(pd, aes(x = delta_fwet, y = dF_col)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  geom_smooth(method = "lm", color = "#1b6ca8", fill = "#cfe3f3", linewidth = 0.6) +
  geom_point(aes(color = canopy_height, size = canopy_height)) +
  ggrepel::geom_text_repel(aes(label = site), size = 2.2, max.overlaps = 14) +
  facet_wrap(~ basis, scales = "free_y") +
  scale_color_viridis_c(option = "viridis", end = 0.9, name = "canopy (m)") +
  scale_size_continuous(range = c(1.4, 4.5), guide = "none") +
  labs(title = "Bottom-to-top change in CH4 flux vs change in footprint wetland",
       subtitle = "Night, PAR-based, 10 km buffer. Wetland mechanism predicts a rising line through the origin.",
       x = "delta f_wet  (top sensor - bottom sensor)",
       y = "F(top) - F(bottom)   (nmol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "footprint", "FINAL_dFcolumn_vs_dfwet.png"), p,
       width = 12, height = 6, dpi = 150)

message("\nReport: ", file.path(tab_dir, "vertical_change_report.txt"))
