# config.R
# Central configuration for height-resolved concentration profile analysis

# -- Paths --
# Every path can be overridden with an environment variable, following the
# convention in the LTER flux-gradient repos (Sys.getenv("LOCALDIR_...")).
# Scripts source this file from the project root, so that is the default.
#
#   FLUXPROFILES_DIR   project root (default: the working directory)
#   FLUXGRADIENT_REPO  lterwg-flux-gradient checkout (default: sibling folder)
#   DATA_DIRS          data directories, separated by ";". Give several to stack
#                      consecutive NEON download periods into one record, e.g.
#                      .../NEON_Aligned_Concentrations;.../NEON_Aligned_Concentrations_202407_202608
#
# Keep data outside a Google Drive-synced folder if you can: downloads written
# into the mounted Drive folder have silently come out as 0-byte files.
proj_dir <- normalizePath(Sys.getenv("FLUXPROFILES_DIR", unset = getwd()),
                          mustWork = FALSE)
if (!file.exists(file.path(proj_dir, "config.R"))) {
  stop("config.R not found in ", proj_dir, ". Run from the project root or set ",
       "FLUXPROFILES_DIR.", call. = FALSE)
}
parent_repo <- Sys.getenv("FLUXGRADIENT_REPO",
                          unset = file.path(dirname(proj_dir), "lterwg-flux-gradient"))
data_dir <- {
  dd <- Sys.getenv("DATA_DIRS", unset = "")
  if (nzchar(dd)) trimws(strsplit(dd, ";", fixed = TRUE)[[1]]) else file.path(proj_dir, "data")
}
missing_dirs <- data_dir[!dir.exists(data_dir)]
if (length(missing_dirs)) {
  warning("Data director", if (length(missing_dirs) > 1) "ies" else "y",
          " not found: ", paste(missing_dirs, collapse = ", "), call. = FALSE)
}
fig_dir <- file.path(proj_dir, "output", "figures")
tab_dir <- file.path(proj_dir, "output", "tables")

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
