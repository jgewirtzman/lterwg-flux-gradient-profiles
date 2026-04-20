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

# -- Site list (all NEON sites with flux-gradient data) --
all_sites <- c(
  "ABBY", "BARR", "BART", "BLAN", "BONA", "CLBJ", "CPER", "DCFS",
  "DEJU", "DELA", "DSNY", "GRSM", "GUAN", "HARV", "HEAL", "JERC",
  "JORN", "KONA", "KONZ", "LAJA", "LENO", "MLBS", "MOAB", "NIWO",
  "NOGP", "OAES", "ONAQ", "ORNL", "OSBS", "PUUM", "RMNP", "SCBI",
  "SERC", "SJER", "SOAP", "SRER", "STEI", "STER", "TALL", "TEAK",
  "TOOL", "TREE", "UKFS", "UNDE", "WOOD", "WREF", "YELL"
)

# -- Forested NEON sites (tall canopy ecosystems where within-canopy
#    profile structure is scientifically interesting) --
forest_sites <- c(
  # Eastern deciduous / mixed
  "HARV", "BART", "SCBI", "SERC", "GRSM", "MLBS", "ORNL",
  "BLAN", "STEI", "TREE", "UKFS", "UNDE",
  # Southeast pine / bottomland
  "DELA", "LENO", "TALL", "OSBS", "JERC",
  # Western evergreen
  "WREF", "ABBY", "SOAP", "TEAK", "SJER",
  # Subalpine / montane
  "NIWO", "RMNP", "YELL",
  # Boreal
  "BONA", "DEJU", "HEAL",
  # Tropical
  "GUAN", "LAJA", "PUUM"
)

# -- Ecosystem type assignments for all NEON terrestrial sites --
# Source: NEON site characteristics, NLCD dominant land cover
ecosystem_types <- c(
  # FOREST - Eastern deciduous / mixed
  HARV = "Forest - Eastern Deciduous",
  BART = "Forest - Eastern Deciduous",
  SCBI = "Forest - Eastern Deciduous",
  SERC = "Forest - Eastern Deciduous",
  GRSM = "Forest - Eastern Deciduous",
  MLBS = "Forest - Eastern Deciduous",
  ORNL = "Forest - Eastern Deciduous",
  BLAN = "Forest - Eastern Deciduous",
  STEI = "Forest - Eastern Deciduous",
  TREE = "Forest - Eastern Deciduous",
  UKFS = "Forest - Eastern Deciduous",
  UNDE = "Forest - Eastern Deciduous",
  # FOREST - Southeast pine / bottomland
  DELA = "Forest - SE Pine/Bottomland",
  LENO = "Forest - SE Pine/Bottomland",
  TALL = "Forest - SE Pine/Bottomland",
  OSBS = "Forest - SE Pine/Bottomland",
  JERC = "Forest - SE Pine/Bottomland",
  # FOREST - Western evergreen
  WREF = "Forest - Western Evergreen",
  ABBY = "Forest - Western Evergreen",
  SOAP = "Forest - Western Evergreen",
  TEAK = "Forest - Western Evergreen",
  SJER = "Forest - Western Evergreen",
  # FOREST - Subalpine / montane
  NIWO = "Forest - Subalpine/Montane",
  RMNP = "Forest - Subalpine/Montane",
  YELL = "Forest - Subalpine/Montane",
  # FOREST - Boreal
  BONA = "Forest - Boreal",
  DEJU = "Forest - Boreal",
  HEAL = "Forest - Boreal",
  # FOREST - Tropical
  GUAN = "Forest - Tropical",
  LAJA = "Forest - Tropical",
  PUUM = "Forest - Tropical",
  # GRASSLAND - Prairie
  CPER = "Grassland",
  KONZ = "Grassland",
  KONA = "Grassland",
  NOGP = "Grassland",
  DCFS = "Grassland",
  WOOD = "Grassland",
  STER = "Grassland",
  DSNY = "Grassland",
  OAES = "Grassland",
  # SHRUBLAND / DESERT
  SRER = "Shrubland/Desert",
  JORN = "Shrubland/Desert",
  MOAB = "Shrubland/Desert",
  ONAQ = "Shrubland/Desert",
  # TUNDRA
  BARR = "Tundra",
  TOOL = "Tundra",
  # WOODLAND / SAVANNA
  CLBJ = "Woodland/Savanna"
)

# -- Packages --
library(tidyverse)
library(lubridate)
library(googledrive)
library(scales)
library(patchwork)
if (requireNamespace("ggrepel", quietly = TRUE)) library(ggrepel)

# -- Source parent repo functions (safely -- some may have missing dependencies) --
parent_functions <- list.files(file.path(parent_repo, "functions"),
                               pattern = "\\.R$", full.names = TRUE)
for (f in parent_functions) {
  tryCatch(source(f), error = function(e) {
    message("  Skipping parent function (missing dependency): ", basename(f))
  })
}

# -- Source this project's R/ functions --
proj_functions <- list.files(file.path(proj_dir, "R"),
                              pattern = "\\.R$", full.names = TRUE)
invisible(lapply(proj_functions, source))

# -- Source this project's plot functions --
plot_functions <- list.files(file.path(proj_dir, "plots"),
                              pattern = "\\.R$", full.names = TRUE)
invisible(lapply(plot_functions, source))
