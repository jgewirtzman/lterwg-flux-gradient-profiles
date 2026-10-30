# lterwg-flux-gradient-profiles

Height-resolved analysis of NEON flux-gradient data: how CH4 (and CO2, H2O)
concentrations and fluxes vary with measurement height, and whether that
vertical structure reflects real sources and sinks — soil uptake below,
emissions from vegetation above — or simply the different footprints that
different heights see.

The project builds on the LTER flux-gradient working group code
([lterwg-flux-gradient](https://github.com/lter/lterwg-flux-gradient),
[-eval](https://github.com/lter/lterwg-flux-gradient-eval),
[-methane](https://github.com/lter/lterwg-flux-gradient-methane)).
Those repos produce one validated ecosystem-scale flux per site. This one
computes fluxes at **every** height.

## Status

Work in progress. All results so far use the original **Aug 2021 – Jun 2024**
record, and most were computed on subsamples. They should be re-run on the
version 3, 2021–2026 record (see [Data](#data)) before they are cited.
Conclusions about the cause of the vertical structure have changed more than
once as estimators and conventions were corrected. Report them with that
sensitivity shown.

Two intended papers:

- **Ecosystem fluxes and their footprints.** How the validated ensemble flux
  at each site reflects the vegetation and wetland in its footprint. This
  uses the LTER v3 ensemble, so the footprint should be computed at the
  heights of each site's selected sensor pair, not at the top of the tower.
- **Vertical structure.** Flux between each pair of heights compared with the
  composition of each height's footprint. This uses this project's own
  per-height fluxes.

## Layout

```
config.R          paths, packages, site metadata; every script sources it
run_all.R         main pipeline: load, fit alpha, per-height fluxes, figures
regen_figs.R      regenerate figures from checkpoints without recomputing
R/                functions
plots/            plotting functions
analysis/
  footprint/      NWI wetland fetch, Kljun footprints, footprint wetland fraction
  vertical/       flux profiles in units of canopy height
  tests/          hypothesis tests (current versions)
diagnostics/      method checks: gradient vs EC, stability, wind direction
archive/          superseded scripts; see archive/README.md for why
```

Run every script **from the project root**, for example
`Rscript analysis/footprint/fetch_all_nwi.R`. Scripts find `config.R` and
`R/` relative to the working directory.

## Setup

R ≥ 4.3 with: tidyverse, lubridate, sf, httr, jsonlite, googledrive, scales,
patchwork, ggrepel, broom, and optionally lme4.

Two things are read from a checkout of
[lterwg-flux-gradient](https://github.com/lter/lterwg-flux-gradient), expected
as a sibling folder:

- `metadata/Ameriflux_NEON field-sites.csv`: IGBP class, climate and
  coordinates for the 47 sites. This is the same table the methane workflow
  uses, so site classes match the methane paper.
- `exploratory/fluxfootprint_kljun.R`: the Kljun et al. (2015) footprint
  model. It isn't copied into this repo because the file carries a copyright
  notice but no license allowing redistribution. Please cite
  Kljun, N., Calanca, P., Rotach, M. W. and Schmid, H. P. (2015), *Geosci.
  Model Dev.* 8, 3695–3713.

Nothing else is loaded from the parent repo.

### Paths

Each path can be overridden with an environment variable:

| variable | meaning | default |
|---|---|---|
| `FLUXPROFILES_DIR` | project root | working directory |
| `FLUXGRADIENT_REPO` | lterwg-flux-gradient checkout | sibling folder |
| `DATA_DIRS` | data directories, separated by `;` | `./data` |
| `FLUXFOOT_DIR` | NEON 30-min products (wind direction) | `~/NEON/fluxgrad_min30` |
| `SITE_METADATA_FILE` | site metadata table | in `FLUXGRADIENT_REPO` |
| `KLJUN_FFP_FILE` | Kljun FFP source | in `FLUXGRADIENT_REPO` |

Keep data outside a Google Drive–synced folder if you can. Downloads written
into the mounted Drive folder have come out as 0-byte files with no error.

## Data

Each data directory holds per-site folders: `<dir>/<SITE>/<SITE>_<suffix>.RData`.

| file | object | contents |
|---|---|---|
| `_9min` | `min9.list` | per-level concentrations |
| `_aligned_conc_flux_9min` | `min9Diff.list` | adjacent-pair concentration differences with met and EC flux |
| `_attr` | `attr.df` | tower and inlet heights, canopy height, coordinates |
| `_30min` | `min30.list` | includes `FluxFoot$angZaxsErth` (wind direction) and `veloYaxsHorSd` (σv) |

**Version 3 (October 2026).** The 2021–2026 record comes in two download
periods, `NEON_Aligned_Concentrations` and
`NEON_Aligned_Concentrations_202407_202608`. List both in `DATA_DIRS` and
the loader stacks them, dropping duplicate windows at the boundary with the
same rule as `stack.periods()` in the upstream flux workflow. Version 3 lives
on the Malone Lab server and isn't on the shared Google Drive, which was last
modified on 2026-06-30. `R/download_data.R` can therefore fetch only the
original 2021–2024 data.

Version 3 changed several things upstream that affect this project:
`L_obukhov` now uses k = 0.4 and is reset to the stability-correction L
inside the flux step, so the aligned files still carry the old value; zeta
is now computed per pair; and the ensemble flux is assigned to the half-hour
of the sampling midpoint. See the upstream README for the full list.

Wetland: [NWI](https://www.fws.gov/program/national-wetlands-inventory)
polygons, fetched by `analysis/footprint/fetch_all_nwi.R` into
`data/nwi/` (5 km) or `data/nwi_10km/`. The methane workflow uses WAD2M
(0.5°) and MODIS MCD12C1 (0.05°). Both are too coarse for tower footprints,
which cover 0.1 to a few km², so NWI is the only footprint-scale product
covering every NEON site. NWI is static.

## Conventions and pitfalls

- **Timestamps are UTC.** Never classify day and night or bin diel curves
  with `hour()`. Use `add_tod_par()` where PAR is available. Otherwise use
  `add_tod(df, lon = …)`, `local_hour()` or `is_core_night()`, which work in
  local solar time. Before this was fixed, about half of "night" records
  were daylight.
- **Three definitions of night are in use**: 18:00–06:00 in `add_tod()`,
  22:00–06:00 in `plot_cross_site_synthesis.R`, and 22:00–05:00 in
  `filter_pair_metrics()`. They haven't been reconciled.
- **IGBP class and canopy height disagree at some sites**: NIWO (ENF,
  0.6 m canopy), CLBJ (GRA, 16 m) and BLAN (DBF, 4 m). IGBP describes the
  landscape and `DistZaxsCnpy` describes the canopy at the tower. Choose the
  one that matches the question.
- **Units.** CH4 fluxes are computed in µmol m⁻² s⁻¹ and plotted ×1000 as
  nmol m⁻² s⁻¹.
- **Pairing levels.** The top and bottom inlets almost never pass QC in the
  same half-hour (a median of 19 shared timestamps per site), so
  top-minus-bottom estimates can't be paired. Stratify by season × day/night
  instead. Never difference two medians computed over different samples, and
  never sum per-timestamp layer differences without checking how often the
  full column is present.
- **Footprint rotation.** The Kljun footprint is the upwind source area.
  Native +x points toward the bearing the wind comes from, at a math angle of
  90 − wind_dir. This has been checked against Kljun's own rotation to
  0.00 m. If you change it, test against an *asymmetric* synthetic target: a
  symmetric test passes even when the rotation is reflected.
