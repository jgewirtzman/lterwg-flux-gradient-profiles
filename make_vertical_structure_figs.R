# make_vertical_structure_figs.R
# Graphs for the core question: CH4 flux is more positive at the top than the
# bottom — is that mostly tall sites, and does the day/night split matter?
#
# Every panel shows the same quantity computed three ways (all data, day, night)
# and on two bases (all layers vs within-canopy layers excluded), so the effect
# of each analytical choice is visible rather than asserted.
#
# Usage:  Rscript make_vertical_structure_figs.R

source("config.R")
suppressPackageStartupMessages({ library(ggplot2); library(lubridate); library(tidyr) })

ck_dir <- file.path(proj_dir, "output", "checkpoints")
cm <- read.csv(file.path(tab_dir, "canopy_metrics.csv"), stringsAsFactors = FALSE)
sites <- sub("\\.rds$", "", basename(list.files(ck_dir, pattern = "\\.rds$")))
fp <- file.path(fig_dir, "footprint")

col_rows <- list(); prof_rows <- list(); ww_rows <- list()

for (s in sites) {
  chm <- cm$canopy_height[cm$site == s]
  if (!length(chm) || is.na(chm)) next
  x <- tryCatch(readRDS(file.path(ck_dir, paste0(s, ".rds"))), error = function(e) NULL)
  if (is.null(x) || is.null(x$met) || !"PAR" %in% names(x$met)) next
  par <- x$met %>% select(t = timeMid, PAR) %>% filter(!is.na(PAR))

  # ---- column change: F(top) - F(bottom), summed over adjacent layers ----
  ss <- x$source_sink$CH4
  if (!is.null(ss) && nrow(ss) > 100) {
    ss <- ss %>% filter(is.finite(dF)) %>%
      mutate(WW = (z_upper <= chm) & (z_lower <= chm)) %>%
      left_join(par, by = c("time_round" = "t")) %>%
      mutate(tod = ifelse(is.na(PAR), NA, ifelse(PAR > 10, "Day", "Night")))
    for (ex in c(FALSE, TRUE)) {
      a <- if (ex) ss %>% filter(!WW) else ss
      for (td in c("All", "Day", "Night")) {
        b <- if (td == "All") a else a %>% filter(tod == td)
        if (nrow(b) < 100) next
        col <- b %>% group_by(time_round) %>%
          summarise(dF_col = sum(dF), .groups = "drop")
        if (nrow(col) < 50) next
        col_rows[[length(col_rows) + 1]] <- tibble(
          site = s, canopy_height = chm, excl_WW = ex, tod = td,
          dF_col = median(col$dF_col), n_ts = nrow(col),
          z_span = max(b$z_upper) - min(b$z_lower))
      }
    }
  }

  # ---- flux profile by relative height ----
  f <- x$fluxes$CH4
  if (!is.null(f) && nrow(f) > 100) {
    f <- f %>% filter(F_outlier %in% c(FALSE, 0), is.finite(F_FG), !is.na(z_mid)) %>%
      left_join(par, by = c("time_round" = "t")) %>%
      mutate(tod = ifelse(is.na(PAR), NA, ifelse(PAR > 10, "Day", "Night")),
             WW  = (z_above <= chm) & (z <= chm))
    for (td in c("All", "Day", "Night")) {
      b <- if (td == "All") f else f %>% filter(tod == td)
      if (nrow(b) < 100) next
      prof_rows[[length(prof_rows) + 1]] <- b %>%
        group_by(z_mid) %>%
        summarise(n = n(), F = median(F_FG), frac_WW = mean(WW), .groups = "drop") %>%
        filter(n >= 20) %>%
        mutate(site = s, tod = td, canopy_height = chm,
               z_rel = z_mid / chm)
    }
    ww_rows[[length(ww_rows) + 1]] <- f %>%
      group_by(z_mid) %>% summarise(n = n(), frac_WW = mean(WW), .groups = "drop") %>%
      mutate(site = s, canopy_height = chm, z_rel = z_mid / chm)
  }
}

d  <- bind_rows(col_rows)
pr <- bind_rows(prof_rows)
write.csv(d, file.path(tab_dir, "vertical_structure_columns.csv"), row.names = FALSE)

lvl <- function(v) factor(v, levels = c("All", "Day", "Night"))
d$tod  <- lvl(d$tod);  pr$tod <- lvl(pr$tod)
d$basis  <- factor(ifelse(d$excl_WW, "within-canopy layers excluded", "all layers"),
                   levels = c("all layers", "within-canopy layers excluded"))

# ============================ FIG 1 ============================
# Is it mostly tall sites? And does day/night change the answer?
p1 <- ggplot(d, aes(x = canopy_height, y = dF_col)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_smooth(method = "lm", color = "#1b6ca8", fill = "#cfe3f3", linewidth = 0.6) +
  geom_point(aes(color = dF_col > 0), size = 1.9, alpha = 0.9) +
  facet_grid(basis ~ tod, scales = "free_y") +
  scale_color_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2c6fad"),
                     labels = c("top more negative", "top more positive"), name = NULL) +
  labs(title = "Is the top-minus-bottom CH4 flux difference a tall-site phenomenon?",
       subtitle = "Each point is a site. Positive = flux more positive at the top level than the bottom.",
       x = "canopy height (m)", y = "F(top) - F(bottom)   (nmol m-2 s-1)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(fp, "VS1_dFcolumn_vs_canopy.png"), p1, width = 12, height = 7, dpi = 150)

# ============================ FIG 2 ============================
# The vertical profile itself, in units of canopy height.
pr2 <- pr %>% filter(canopy_height >= 5, is.finite(z_rel), z_rel <= 3)
bins <- pr2 %>%
  mutate(zb = cut(z_rel, breaks = c(0, .25, .5, .75, 1, 1.5, 2, 3),
                  labels = c("0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1",
                             "1-1.5", "1.5-2", "2-3"))) %>%
  filter(!is.na(zb)) %>%
  group_by(tod, zb) %>%
  summarise(n_site = n_distinct(site), med = median(F),
            q25 = quantile(F, .25), q75 = quantile(F, .75), .groups = "drop")
p2 <- ggplot(bins, aes(x = med, y = zb, group = tod)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = 2) +
  geom_errorbarh(aes(xmin = q25, xmax = q75), height = 0, color = "grey65") +
  geom_path(color = "grey40") +
  geom_point(aes(fill = med > 0), shape = 21, size = 3, color = "grey20") +
  facet_wrap(~ tod) +
  scale_fill_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2c6fad"),
                    labels = c("uptake", "emission"), name = NULL) +
  labs(title = "CH4 flux profile in units of canopy height (tall sites, canopy >= 5 m)",
       subtitle = "Median across sites, bars = inter-site IQR. Bins at or below 1 are within the canopy.",
       x = "median F_FG (nmol m-2 s-1)", y = "height / canopy height") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "VS2_profile_relative_height.png"), p2, width = 11, height = 5.5, dpi = 150)

# ============================ FIG 3 ============================
# Does the day/night split matter? Paired, same site.
w <- d %>% filter(tod != "All") %>%
  select(site, canopy_height, basis, tod, dF_col) %>%
  pivot_wider(names_from = tod, values_from = dF_col) %>%
  filter(is.finite(Day), is.finite(Night))
p3 <- ggplot(w, aes(x = Day, y = Night)) +
  geom_hline(yintercept = 0, color = "grey55", linetype = 2) +
  geom_vline(xintercept = 0, color = "grey55", linetype = 2) +
  geom_abline(slope = 1, intercept = 0, color = "firebrick", linetype = 3) +
  geom_point(aes(color = canopy_height), size = 2.3) +
  facet_wrap(~ basis, scales = "free") +
  scale_color_viridis_c(option = "viridis", end = .9, name = "canopy (m)") +
  labs(title = "Does the day/night split matter?",
       subtitle = paste("Same site, same quantity, computed on daytime vs nighttime data.",
                        "\nPoints off the dotted 1:1 line mean the split changes the answer."),
       x = "F(top) - F(bottom), Day", y = "F(top) - F(bottom), Night") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "VS3_day_vs_night.png"), p3, width = 11, height = 5.5, dpi = 150)

# ============================ FIG 4 ============================
# How much of the profile is within-canopy?
ww <- bind_rows(ww_rows) %>% filter(canopy_height >= 5, is.finite(z_rel), z_rel <= 3)
p4 <- ggplot(ww, aes(x = z_rel, y = frac_WW)) +
  geom_vline(xintercept = 1, color = "firebrick", linetype = 2) +
  geom_point(aes(size = n), alpha = .35, color = "#2c6fad") +
  scale_size_continuous(range = c(.5, 4), guide = "none") +
  labs(title = "Which measurement levels sit inside the canopy?",
       subtitle = "Red line = canopy top. Points left of it are levels whose pairs are largely within-canopy, where gradient theory does not hold.",
       x = "measurement height / canopy height",
       y = "fraction of that level's pairs classed within-canopy") +
  theme_minimal(base_size = 11)
ggsave(file.path(fp, "VS4_within_canopy_fraction.png"), p4, width = 9, height = 5, dpi = 150)

# ---- numbers to accompany the figures ----
sink(file.path(tab_dir, "vertical_structure_summary.txt"), split = TRUE)
cat("=== F(top) - F(bottom): does the analytical choice change the answer? ===\n\n")
print(as.data.frame(d %>% group_by(basis, tod) %>%
  summarise(n_sites = n(), median = signif(median(dF_col), 3),
            frac_positive = round(mean(dF_col > 0), 2),
            p_vs_zero = signif(t.test(dF_col)$p.value, 3), .groups = "drop")),
  row.names = FALSE)
cat("\n=== correlation with canopy height ===\n\n")
print(as.data.frame(d %>% group_by(basis, tod) %>%
  summarise(n = n(), r = round(cor(dF_col, canopy_height, use = "complete.obs"), 3),
            p = signif(cor.test(dF_col, canopy_height)$p.value, 3), .groups = "drop")),
  row.names = FALSE)
sink()

message("\nFigures: ", fp)
