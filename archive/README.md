# Archived scripts

Superseded scripts, kept for the record. Do not cite their outputs. Each
was replaced because of a specific flaw, listed below, so do not revive one
without fixing it first.

| script | flaw | superseded by |
|---|---|---|
| `test_flux_height_vs_wetland.R` | The top-minus-bottom flux difference was not normalised by the height range the tower spans. Taller towers span more height, so they mechanically gave bigger differences, and the "more wetland" sites were also the taller ones. | `analysis/tests/test_flux_height_vs_wetland2.R` (per-metre slope) |
| `test_vertical_change_vs_wetland.R` | Estimated F(top) − F(bottom) as the median over timestamps of the summed layer dF. That sum covers only the layers that passed QC at each timestamp, and the full column is present at only about 1% of timestamps (13 of ~5,700 at HARV). The estimate therefore came from a short, shifting part of the column and was about a third of the true magnitude. Its conclusion, that canopy height predicts the vertical change and wetland does not, was reversed under a correct estimator. | Stratified difference of per-level medians (within season × day/night cells); see the project memory notes |
| `test_ecosystem_flux_vs_wetland.R` | Built on `aggregate_flux_profiles()` while it still classified day and night by UTC hour. About half of "night" records were daylight, up to 96% at PUUM. | `analysis/tests/test_wetland_final.R` (PAR-based day/night) |
| `test_step7_subset.R` | A plumbing test with an empty mock wetland layer, written before NWI could be fetched. Δf_wet is identically 0 by design, so its regression outputs mean nothing. | The real NWI pipeline in `analysis/footprint/` |
