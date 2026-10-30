# test_ecosystem_flux_vs_wetland.R
# Does the ecosystem-level CH4 source scale with footprint wetland?
#
# Framed the way the project already frames it, and on the project's own
# filtered basis (JJA Night, via aggregate_flux_profiles), rather than the
# all-conditions stability stratification used in the earlier probe.
#
# Established pattern at tall sites:
#   CO2 : positive at ground level (Rs) -> negative at ecosystem level (GPP-Rs)
#   CH4 : negative at ground level (soil uptake) -> POSITIVE at ecosystem level
#         (source unknown — this is the question)
#
# So the target quantity is the above-canopy CH4 flux, and the direct test is
# whether it scales with the wetland fraction the upper sensor actually sees.
# Note this asks about f_wet at the top sensor, not delta f_wet: if the
# ecosystem-level source is wetland in the footprint, the flux should track how
# much wetland the upper footprint contains.
#
# CO2 is carried through as a STRUCTURAL reference, not a bias reference: the
# gases genuinely have different vertical source distributions, so subtracting
# one from the other would be wrong. CO2's job here is to confirm the layer
# decomposition reproduces the known Rs / NEE structure.
#
# Usage:  Rscript test_ecosystem_flux_vs_wetland.R

source("config.R")
source(file.path(proj_dir, "plots", "plot_cross_site_synthesis.R"))
suppressPackageStartupMessages({ library(ggplot2); library(lubridate) })

ck_dir <- file.path(proj_dir, "output", "checkpoints")
sites  <- sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$")))

fw <- read.csv(file.path(tab_dir, "footprint_fwet_summary.csv"), stringsAsFactors = FALSE)
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE)

# ---- per-site JJA-Night flux profiles, using the project's own aggregator ----
prof <- list()
for (gas in c("CH4", "CO2")) {
  lf_all <- list()
  for (site in sites) {
    x <- tryCatch(readRDS(file.path(ck_dir, paste0(site, ".rds"))), error = function(e) NULL)
    if (is.null(x) || is.null(x$fluxes[[gas]])) next
    f <- x$fluxes[[gas]]
    if (!all(c("time_round", "z", "F_FG") %in% names(f))) next
    lf_all[[site]] <- f %>% filter(F_outlier %in% c(FALSE, 0))
  }
  p <- aggregate_flux_profiles(lf_all, filter_jja_night = TRUE)
  if (nrow(p) == 0) next
  p$gas <- gas
  prof[[gas]] <- p
}
prof <- bind_rows(prof)
write.csv(prof, file.path(tab_dir, "jja_night_flux_profiles.csv"), row.names = FALSE)

# ---- ground level vs ecosystem level -------------------------------------
# "ecosystem level" = highest sensor (above canopy where one exists);
# "ground level"    = lowest sensor.
ge <- prof %>%
  filter(n >= 20) %>%
  group_by(gas, site) %>%
  arrange(z, .by_group = TRUE) %>%
  summarise(n_lev = n(),
            z_ground = first(z), z_eco = last(z),
            F_ground = first(F_med), F_eco = last(F_med),
            dF = last(F_med) - first(F_med),
            n_ground = first(n), n_eco = last(n),
            .groups = "drop") %>%
  left_join(cm %>% select(site, canopy_height, canopy_class), by = "site") %>%
  mutate(tall = canopy_height >= 5)

sink(file.path(tab_dir, "ecosystem_flux_vs_wetland.txt"), split = TRUE)
cat("======================================================================\n")
cat(" Ecosystem-level CH4 source vs footprint wetland   (JJA Night)\n")
cat("======================================================================\n")

cat("\n--- 0. Does the established pattern reproduce? (tall sites, canopy >= 5 m) ---\n")
for (g in c("CO2", "CH4")) {
  t <- ge %>% filter(gas == g, tall)
  if (nrow(t) == 0) next
  cat(sprintf("\n  %s  (n = %d tall sites)\n", g, nrow(t)))
  cat(sprintf("    ground level : median %+.4g   | positive at %d/%d sites\n",
              median(t$F_ground), sum(t$F_ground > 0), nrow(t)))
  cat(sprintf("    ecosystem lvl: median %+.4g   | positive at %d/%d sites\n",
              median(t$F_eco), sum(t$F_eco > 0), nrow(t)))
}
cat("\n  Expected: CO2 positive at ground (Rs), negative aloft (GPP-Rs).\n")
cat("            CH4 negative at ground (soil uptake), positive aloft (?).\n")

# ---- the test -------------------------------------------------------------
d <- ge %>% filter(gas == "CH4") %>%
  left_join(fw %>% select(site, fwet_top, fwet_bot, delta_fwet,
                          f_wet_buffer, status, beyond_buffer), by = "site")

cat("\n--- 1. Does the ECOSYSTEM-level CH4 flux scale with what the upper\n")
cat("       footprint actually sees? ---\n")
shw <- function(lbl, fm, dat = d) {
  fit <- tryCatch(lm(fm, data = dat), error = function(e) NULL)
  if (is.null(fit)) return(invisible())
  s <- summary(fit)
  cat("\n  ", lbl, "  (n = ", nobs(fit), ")\n", sep = "")
  print(round(s$coefficients, 5))
  cat(sprintf("   adj R^2 = %.3f\n", s$adj.r.squared))
}
shw("F_eco ~ fwet_top",                       F_eco ~ fwet_top)
shw("F_eco ~ fwet_top + canopy_height",       F_eco ~ fwet_top + canopy_height)
shw("F_eco ~ f_wet_buffer",                   F_eco ~ f_wet_buffer)
shw("dF (eco - ground) ~ delta_fwet",         dF ~ delta_fwet)
shw("dF ~ delta_fwet + canopy_height",        dF ~ delta_fwet + canopy_height)

cat("\n--- 2. Restricted to tall sites (canopy >= 5 m), where the\n")
cat("       ground/ecosystem contrast is meaningful ---\n")
dt <- d %>% filter(tall)
shw("F_eco ~ fwet_top  [tall only]",          F_eco ~ fwet_top,  dt)
shw("F_eco ~ fwet_top + canopy_height [tall]", F_eco ~ fwet_top + canopy_height, dt)

cat("\n--- 3. Excluding sites whose footprints outrun the 5 km NWI buffer ---\n")
dc <- d %>% filter(!beyond_buffer)
shw("F_eco ~ fwet_top  [buffer-clean]",       F_eco ~ fwet_top,  dc)

cat("\n--- 4. Sites with positive ecosystem-level CH4: are they the wet ones? ---\n")
pos <- d %>% filter(is.finite(F_eco))
if (nrow(pos) > 4) {
  pos <- pos %>% mutate(eco_source = F_eco > 0)
  print(as.data.frame(pos %>% group_by(eco_source) %>%
    summarise(n = n(), median_fwet_top = round(median(fwet_top, na.rm = TRUE), 4),
              median_canopy = round(median(canopy_height, na.rm = TRUE), 1),
              .groups = "drop")), row.names = FALSE)
  tt <- tryCatch(t.test(fwet_top ~ eco_source, data = pos), error = function(e) NULL)
  if (!is.null(tt)) cat(sprintf("\n  f_wet_top, source vs sink sites: p = %.3g\n", tt$p.value))
}

cat("\n--- 5. Caveats ---\n")
cat("  JJA Night only, matching the project's other flux figures.\n")
cat("  f_wet is a symmetric-wind estimate; see the wind-direction diagnostic.\n")
cat(sprintf("  %d of %d sites have footprints exceeding the 5 km NWI buffer.\n",
            sum(d$beyond_buffer, na.rm = TRUE), nrow(d)))
sink()

# ---- figures ---------------------------------------------------------------
fp <- file.path(fig_dir, "footprint")

p1 <- ggplot(d %>% filter(is.finite(F_eco), is.finite(fwet_top)),
             aes(x = fwet_top, y = F_eco)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_smooth(method = "lm", color = "#1b6ca8", fill = "#cfe3f3", linewidth = 0.6) +
  geom_point(aes(size = canopy_height, color = canopy_height)) +
  ggrepel::geom_text_repel(aes(label = site), size = 2.4, max.overlaps = 18) +
  scale_color_viridis_c(option = "viridis", end = 0.9, name = "canopy (m)") +
  scale_size_continuous(range = c(1.5, 5), guide = "none") +
  labs(title = "Ecosystem-level CH4 flux vs wetland fraction in the upper footprint",
       subtitle = "JJA Night. Positive values are the unexplained above-canopy source.",
       x = "f_wet seen by the top sensor",
       y = "F_FG at ecosystem level (nmol m-2 s-1)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "TEST5_Feco_vs_fwet_top.png"), p1, width = 9.5, height = 6.5, dpi = 150)

gg <- ge %>% filter(tall, is.finite(F_ground), is.finite(F_eco)) %>%
  select(gas, site, F_ground, F_eco) %>%
  tidyr::pivot_longer(c(F_ground, F_eco), names_to = "level", values_to = "F") %>%
  mutate(level = factor(ifelse(level == "F_ground", "ground", "ecosystem"),
                        levels = c("ground", "ecosystem")))
p2 <- ggplot(gg, aes(x = level, y = F, group = site)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_line(alpha = 0.4) + geom_point(size = 1.6, alpha = 0.8) +
  facet_wrap(~ gas, scales = "free_y") +
  labs(title = "Ground vs ecosystem level, tall sites (JJA Night)",
       subtitle = "CO2 goes positive->negative (Rs then GPP-Rs). CH4 goes negative->positive: the open question.",
       x = NULL, y = "median F_FG") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "TEST5_ground_vs_ecosystem.png"), p2, width = 9, height = 5.5, dpi = 150)

message("\nReport: ", file.path(tab_dir, "ecosystem_flux_vs_wetland.txt"))
