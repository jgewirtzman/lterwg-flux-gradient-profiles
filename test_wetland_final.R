# test_wetland_final.R
# Does the ecosystem-level CH4 source scale with footprint wetland?
#
# Run on the corrected basis established over this analysis:
#   - PAR-based day/night          (time_round is UTC; hour-based splits are wrong)
#   - 10 km NWI buffer             (at 5 km, 14 sites had footprints running past it)
#   - paired within-timestamp      (never difference two independently-computed medians)
#   - Canopy_L1 == "WW" pairs both included and excluded, reported side by side,
#     because excluding them removes the below-canopy profile the ground-level
#     uptake signal depends on, while improving the CO2 control.
#
# Flux QC (ustar >= 0.1, dConcSNR >= 3) is already applied upstream in
# compute_fluxes.R and matches the parent repos' filter_fluxes().
#
# Usage:  Rscript test_wetland_final.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(lubridate) })

NWI_KM <- 10
ck_dir <- file.path(proj_dir, "output", "checkpoints")
fw <- read.csv(file.path(tab_dir, sprintf("footprint_fwet_summary_%dkm.csv", NWI_KM)),
               stringsAsFactors = FALSE)
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE)

sites <- intersect(sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$"))),
                   fw$site)
message("sites: ", length(sites))

rows <- list()
for (s in sites) {
  chm <- cm$canopy_height[cm$site == s]
  if (!length(chm) || is.na(chm)) next
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(s, ".rds"))), error = function(e) NULL)
  if (is.null(x) || is.null(x$fluxes$CH4) || is.null(x$met)) next
  if (!"PAR" %in% names(x$met)) next
  par <- x$met %>% select(t = timeMid, PAR) %>% filter(!is.na(PAR))

  f <- x$fluxes$CH4 %>%
    filter(F_outlier %in% c(FALSE, 0), is.finite(F_FG), !is.na(z_mid)) %>%
    mutate(Canopy_L1 = paste0(ifelse(z_above > chm, "A", "W"),
                              ifelse(z      > chm, "A", "W"))) %>%
    inner_join(par, by = c("time_round" = "t")) %>%
    mutate(tod = ifelse(PAR > 10, "Day", "Night"))
  if (nrow(f) < 200) next

  for (excl in c(FALSE, TRUE)) {
    a <- if (excl) f %>% filter(Canopy_L1 != "WW") else f
    if (nrow(a) < 100) next
    for (td in c("Day", "Night")) {
      b <- a %>% filter(tod == td)
      if (nrow(b) < 100) next
      zt <- max(b$z_mid); zb <- min(b$z_mid)
      if (zt <= zb) next
      # ecosystem-level flux: median at the top level
      Feco <- median(b$F_FG[b$z_mid == zt], na.rm = TRUE)
      # canopy layer: paired within timestamp, then median
      pr <- b %>% filter(z_mid %in% c(zb, zt)) %>%
        select(time_round, z_mid, F_FG) %>%
        tidyr::pivot_wider(names_from = z_mid, values_from = F_FG,
                           values_fn = median) %>%
        rename(lo = 2, hi = 3) %>% filter(is.finite(lo), is.finite(hi))
      lay <- if (nrow(pr) >= 30) median(pr$hi - pr$lo) else NA_real_
      rows[[length(rows) + 1]] <- tibble(
        site = s, tod = td, excl_WW = excl, z_top = zt, z_bot = zb,
        F_eco = Feco, canopy_layer = lay, n_pair = nrow(pr),
        canopy_height = chm)
    }
  }
}
d <- bind_rows(rows) %>%
  left_join(fw %>% select(site, fwet_top, delta_fwet, f_wet_buffer,
                          status, beyond_buffer), by = "site")
write.csv(d, file.path(tab_dir, "wetland_final_results.csv"), row.names = FALSE)

sink(file.path(tab_dir, "wetland_final_report.txt"), split = TRUE)
cat("=================================================================\n")
cat(" Ecosystem-level CH4 vs footprint wetland — corrected basis\n")
cat(" PAR day/night | 10 km NWI buffer | paired differences\n")
cat("=================================================================\n")

shw <- function(lbl, fm, dat) {
  dat <- dat %>% filter(is.finite(F_eco), is.finite(fwet_top))
  if (nrow(dat) < 8) { cat("\n  ", lbl, ": too few sites\n", sep = ""); return(invisible()) }
  fit <- tryCatch(lm(fm, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(invisible())
  s <- summary(fit)
  cat("\n  ", lbl, "  (n = ", nobs(fit), ")\n", sep = "")
  print(round(s$coefficients, 5))
  cat(sprintf("   adj R^2 = %.3f\n", s$adj.r.squared))
}

for (ex in c(FALSE, TRUE)) {
  cat("\n\n############ ", if (ex) "EXCLUDING within-canopy (WW) pairs"
                          else "ALL pairs (WW included)", " ############\n", sep = "")
  for (td in c("Night", "Day")) {
    sub <- d %>% filter(excl_WW == ex, tod == td)
    cat("\n--- ", td, " ---\n", sep = "")
    cat(sprintf("  sites: %d | median F_eco %+.4g | median canopy layer %+.4g\n",
                nrow(sub), median(sub$F_eco, na.rm = TRUE),
                median(sub$canopy_layer, na.rm = TRUE)))
    shw("F_eco ~ fwet_top", F_eco ~ fwet_top, sub)
    shw("F_eco ~ fwet_top + canopy_height", F_eco ~ fwet_top + canopy_height, sub)
    shw("canopy_layer ~ delta_fwet + canopy_height",
        canopy_layer ~ delta_fwet + canopy_height, sub)
  }
}

cat("\n\n--- caveats ---\n")
cat(sprintf("  sites still exceeding the 10 km buffer: %d\n",
            sum(d$beyond_buffer[!duplicated(d$site)], na.rm = TRUE)))
cat("  f_wet is a symmetric-wind estimate — see the wind-direction diagnostic;\n")
cat("  under rotational averaging f_wet is 83-96%% footprint SIZE, so a positive\n")
cat("  fwet_top coefficient cannot yet be cleanly separated from footprint size.\n")
sink()

# figure: the headline relationship, both bases, night
pd <- d %>% filter(tod == "Night", is.finite(F_eco), is.finite(fwet_top)) %>%
  mutate(basis = ifelse(excl_WW, "excluding WW pairs", "all pairs"))
p <- ggplot(pd, aes(x = fwet_top, y = F_eco)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_smooth(method = "lm", color = "#1b6ca8", fill = "#cfe3f3", linewidth = 0.6) +
  geom_point(aes(color = canopy_height, size = canopy_height)) +
  ggrepel::geom_text_repel(aes(label = site), size = 2.2, max.overlaps = 14) +
  facet_wrap(~ basis, scales = "free_y") +
  scale_color_viridis_c(option = "viridis", end = 0.9, name = "canopy (m)") +
  scale_size_continuous(range = c(1.4, 4.5), guide = "none") +
  labs(title = "Ecosystem-level CH4 vs wetland fraction in the upper footprint",
       subtitle = "Night, PAR-based. 10 km NWI buffer. Left: all pairs. Right: within-canopy pairs excluded.",
       x = "f_wet seen by the top sensor",
       y = "F_FG at ecosystem level (nmol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "footprint", "FINAL_Feco_vs_fwet.png"), p,
       width = 12, height = 6, dpi = 150)

message("\nReport: ", file.path(tab_dir, "wetland_final_report.txt"))
