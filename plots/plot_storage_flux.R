# plot_storage_flux.R
# Diel patterns of storage flux (dC/dt) by tower level

#' Plot diel pattern of storage flux at each height
#'
#' @param storage_summary Tibble. Output of summarize_storage_flux_diel()
#' @param gas Character. Gas name
#' @param site_name Character. Site code for title
#' @return ggplot object
plot_storage_flux_diel <- function(storage_summary, gas = "CH4",
                                    site_name = NULL) {

  if (is.null(storage_summary) || nrow(storage_summary) == 0) {
    return(ggplot() + theme_void() + labs(title = "No data"))
  }

  df <- storage_summary %>%
    mutate(height_label = paste0(round(height_m, 1), " m"))

  title_text <- paste0(gas, " Storage Flux (dC/dt) by Height")
  if (!is.null(site_name)) title_text <- paste0(site_name, ": ", title_text)

  ggplot(df, aes(x = hour, y = dCdt_mean,
                  color = reorder(height_label, height_m),
                  group = height_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = dCdt_mean - dCdt_se, ymax = dCdt_mean + dCdt_se,
                    fill = reorder(height_label, height_m)),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.2) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    scale_color_viridis_d(option = "D", end = 0.9, name = "Height") +
    scale_fill_viridis_d(option = "D", end = 0.9, guide = "none") +
    facet_wrap(~ season, scales = "free_y") +
    labs(
      x = "Hour of day",
      y = paste0("dC/dt (", gas_units(gas), "/hr)"),
      title = title_text,
      subtitle = "Positive = accumulation, Negative = depletion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}
