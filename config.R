# config.R
# Central configuration for height-resolved concentration profile analysis

# -- Paths --
proj_dir    <- "/Users/jongewirtzman/My Drive/Research/Flux Gradient/height-conc-analysis"
parent_repo <- "/Users/jongewirtzman/My Drive/Research/Flux Gradient/lterwg-flux-gradient"
data_dir    <- file.path(proj_dir, "data")
fig_dir     <- file.path(proj_dir, "output", "figures")
tab_dir     <- file.path(proj_dir, "output", "tables")

# -- Google Drive --
gdrive_folder_id <- "1Q99CT77DnqMl2mrUtuikcY47BFpckKw3"

# -- Packages --
library(tidyverse)
library(lubridate)
library(googledrive)
library(scales)
library(patchwork)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)

# -- Source this project's R/ functions --
proj_functions <- list.files(file.path(proj_dir, "R"),
                              pattern = "\\.R$", full.names = TRUE)
invisible(lapply(proj_functions, source))

# -- Source this project's plot functions --
plot_functions <- list.files(file.path(proj_dir, "plots"),
                              pattern = "\\.R$", full.names = TRUE)
invisible(lapply(plot_functions, source))

# -- Load the Kljun (2015) footprint model --
# Nothing else is taken from the parent repo. This project used to source every
# file in its functions/ directory -- 54 functions, none of which are called
# here -- which only added package loads and name-collision risk.
load_kljun(parent_repo)

# -- Site metadata and site lists --
# Derived from the AmeriFlux/NEON table the methane workflow uses, so site
# classes match the methane paper. This replaces hand-coded all_sites,
# forest_sites and ecosystem_types lists; the latter cited NLCD but matched
# neither NLCD nor the measured canopy (see R/site_metadata.R).
SITE_METADATA <- load_site_metadata(parent_repo)
all_sites <- SITE_METADATA$site

# Grouping used for the by-ecotype figures: IGBP class, as in the methane paper.
ecosystem_types <- setNames(paste0(SITE_METADATA$igbp, " - ", SITE_METADATA$igbp_label),
                            SITE_METADATA$site)
