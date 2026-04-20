# download_data.R
# Download NEON flux-gradient data from Google Drive

#' Download data files for a single site from Google Drive
#'
#' @param site Character. NEON site code (e.g., "HARV")
#' @param gdrive_folder_id Character. Google Drive folder ID containing site subfolders
#' @param local_dir Character. Local directory to save data
#' @param file_types Character vector. Which file types to download.
#'   Options: "9min" (raw per-level), "aligned_conc_flux_9min" (paired dConc),
#'   "attr" (site attributes), "aligned_conc_flux_30min"
#' @param overwrite Logical. Re-download if file already exists locally?
#' @return Named list of paths to extracted RData files
download_site_data <- function(site,
                                gdrive_folder_id,
                                local_dir,
                                file_types = c("9min", "aligned_conc_flux_9min", "attr"),
                                overwrite = FALSE) {

  # Create site directory
  site_dir <- file.path(local_dir, site)
  dir.create(site_dir, showWarnings = FALSE, recursive = TRUE)

  # Get Google Drive folder contents
  drive_url <- googledrive::as_id(
    paste0("https://drive.google.com/drive/folders/", gdrive_folder_id)
  )
  data_folder <- googledrive::drive_ls(path = drive_url)

  # Navigate into site subfolder
  site_folder_id <- data_folder$id[data_folder$name == site]
  if (length(site_folder_id) == 0) {
    warning("Site folder not found on Google Drive: ", site)
    return(NULL)
  }
  site_folder <- googledrive::drive_ls(path = site_folder_id)

  # Download each file type
  paths <- list()
  for (ftype in file_types) {
    zip_name <- paste0(site, "_", ftype, ".zip")

    # Check if already downloaded
    rdata_candidates <- list.files(site_dir, pattern = ftype, full.names = TRUE)
    rdata_candidates <- rdata_candidates[grepl("\\.RData$|\\.Rdata$|\\.rdata$", rdata_candidates)]
    if (length(rdata_candidates) > 0 && !overwrite) {
      message("  Already exists: ", zip_name, " -> skipping")
      paths[[ftype]] <- rdata_candidates[1]
      next
    }

    # Find file on Drive
    file_row <- site_folder[site_folder$name == zip_name, ]
    if (nrow(file_row) == 0) {
      warning("  File not found on Drive: ", zip_name)
      next
    }

    # Download
    zip_path <- file.path(site_dir, zip_name)
    message("  Downloading: ", zip_name)
    googledrive::drive_download(file = file_row$id, path = zip_path, overwrite = TRUE)

    # Unzip
    if (grepl("\\.zip$", zip_name)) {
      utils::unzip(zip_path, exdir = site_dir)
      file.remove(zip_path)  # clean up zip after extracting

      # Some zips contain nested paths (e.g., data/SITE/file.Rdata)
      # Move any deeply nested RData files up to site_dir
      nested <- list.files(site_dir, pattern = "\\.RData$|\\.Rdata$|\\.rdata$",
                           recursive = TRUE, full.names = TRUE)
      for (nf in nested) {
        if (dirname(nf) != site_dir) {
          file.rename(nf, file.path(site_dir, basename(nf)))
        }
      }
      # Clean up empty nested dirs
      nested_dirs <- list.dirs(site_dir, recursive = TRUE, full.names = TRUE)
      nested_dirs <- setdiff(nested_dirs, site_dir)
      for (nd in rev(nested_dirs)) {  # reverse to remove deepest first
        if (length(list.files(nd, all.files = TRUE, no.. = TRUE)) == 0) {
          unlink(nd, recursive = TRUE)
        }
      }
    }

    # Record path to extracted file
    rdata_candidates <- list.files(site_dir, pattern = ftype, full.names = TRUE)
    rdata_candidates <- rdata_candidates[grepl("\\.RData$|\\.Rdata$|\\.rdata$", rdata_candidates)]
    if (length(rdata_candidates) > 0) {
      paths[[ftype]] <- rdata_candidates[1]
    }
  }

  paths
}


#' Download data for multiple sites
#'
#' @param sites Character vector of NEON site codes
#' @param gdrive_folder_id Character. Google Drive folder ID
#' @param local_dir Character. Local directory to save data
#' @param file_types Character vector. File types to download
#' @param overwrite Logical. Re-download existing files?
#' @return Named list (by site) of file path lists
download_all_sites <- function(sites,
                                gdrive_folder_id,
                                local_dir,
                                file_types = c("9min", "aligned_conc_flux_9min", "attr"),
                                overwrite = FALSE) {

  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  results <- list()

  for (site in sites) {
    message("Downloading data for: ", site)
    results[[site]] <- tryCatch(
      download_site_data(site, gdrive_folder_id, local_dir, file_types, overwrite),
      error = function(e) {
        warning("Failed to download ", site, ": ", e$message)
        NULL
      }
    )
  }

  results
}
