# HEADER --------------------------------------------
#
# Author: Laura Joan Feyrer
# Email:  ljfeyrer@dal.ca
# Date Updated: 2026-01-15
#
# Script Name: 09_create_kde_contours.R
#
# Description:
## Standalone script to generate multi-quantile contour shapefiles
## from existing KDE rasters. Run AFTER main KDE analysis.
## Keeps main script fast and clean while providing multi-contour option.
#
# Usage:
## Rscript 09_create_kde_contours.R
## Or source("09_create_kde_contours.R")
#

cat("\n=== BATCH MULTI-CONTOUR GENERATION ===\n\n")

# Load libraries
library(sf)
library(terra)
library(dplyr)
library(parallel)

# ========================
# CONFIGURATION
# ========================

# Input/output directories (adjust as needed)
RASTER_DIR <- "output/tif/"
OUTPUT_DIR <- "output/shapes/multi/"

# Contour quantiles to create
QUANTILES <- c(0.85, 0.90, 0.95)

# Parallel processing
USE_PARALLEL <- TRUE
DETECTED_CORES <- detectCores()
N_CORES <- if (is.na(DETECTED_CORES)) 1 else max(1, DETECTED_CORES - 1)

# Overwrite existing files?
OVERWRITE <- F

# ========================
# HELPER FUNCTION
# ========================

create_multi_contour_shapefile <- function(kde_raster, output_path, 
                                           quantiles = c(0.85, 0.90, 0.95),
                                           verbose = TRUE) {
  
  if(verbose) cat("  Processing:", basename(output_path), "\n")
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Calculate all quantile thresholds at once (efficient!)
  kde_values <- values(kde_raster, mat = FALSE, na.rm = TRUE)
  thresholds <- quantile(kde_values, probs = quantiles, type = 6, na.rm = TRUE)
  
  if(verbose) {
    cat("    Thresholds:", paste(round(thresholds, 6), collapse = ", "), "\n")
  }
  
  # Create contour for each quantile
  contours_list <- list()
  
  for(i in seq_along(quantiles)) {
    q <- quantiles[i]
    threshold <- thresholds[i]
    
    # Binary raster for this quantile
    kde_binary <- kde_raster
    kde_binary[kde_binary < threshold] <- NA
    kde_binary[!is.na(kde_binary)] <- 1
    
    # Convert to polygon
    kde_poly <- as.polygons(kde_binary, dissolve = TRUE)
    
    if(!is.null(kde_poly)) {
      kde_sf <- st_as_sf(kde_poly)
      kde_sf <- st_make_valid(kde_sf)
      
      # Add quantile attributes
      kde_sf <- kde_sf %>%
        mutate(
          quantile = q,
          quantile_pct = q * 100,
          quantile_label = paste0("q", q * 100),
          threshold = threshold
        ) %>%
        select(quantile, quantile_pct, quantile_label, threshold, geometry)
      
      contours_list[[i]] <- kde_sf
      
      if(verbose) cat("      q", q*100, ": ", nrow(kde_sf), " polygon(s)\n", sep="")
    }
  }
  
  if(length(contours_list) == 0) {
    if(verbose) cat("    WARNING: No contours created\n")
    return(invisible(NULL))
  }
  
  # Combine all quantiles into single shapefile
  multi_contour_sf <- bind_rows(contours_list)
  
  # Save
  write_sf(multi_contour_sf, output_path, delete_dsn = TRUE)
  
  if(verbose) cat("    Saved:", nrow(multi_contour_sf), "total features\n")
  
  return(invisible(multi_contour_sf))
}

# ========================
# BATCH PROCESSING
# ========================

batch_create_multi_contours <- function(raster_dir, 
                                        output_dir, 
                                        quantiles = c(0.85, 0.90, 0.95),
                                        use_parallel = TRUE,
                                        n_cores = NULL,
                                        overwrite = FALSE,
                                        pattern = "\\.tif$") {
  
  # Create output directory
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Find all KDE rasters
  tif_files <- list.files(raster_dir, pattern = pattern, 
                          full.names = TRUE, recursive = TRUE)
  
  if(length(tif_files) == 0) {
    stop("No raster files found in ", raster_dir)
  }
  
  cat("Found", length(tif_files), "raster files\n")
  cat("Quantiles:", paste(quantiles, collapse = ", "), "\n")
  
  if(use_parallel) {
    if(is.null(n_cores)) n_cores <- max(1, detectCores() - 1)
    cat("Using parallel processing with", n_cores, "cores\n")
  } else {
    cat("Using serial processing\n")
  }
  
  cat("\n")
  
  # Function to process one raster
  process_one <- function(tif_file) {
    # Generate output filename
    base_name <- tools::file_path_sans_ext(basename(tif_file))
    output_file <- file.path(output_dir, paste0(base_name, "_multi.shp"))
    
    # Skip if already exists (unless overwrite = TRUE)
    if(!overwrite && file.exists(output_file)) {
      cat("  Skipping (exists):", basename(output_file), "\n")
      return(list(file = basename(tif_file), status = "skipped"))
    }
    
    tryCatch({
      # Load raster
      kde_raster <- rast(tif_file)
      
      # Create multi-contour shapefile
      result <- create_multi_contour_shapefile(
        kde_raster, 
        output_file, 
        quantiles,
        verbose = !use_parallel  # Suppress verbose in parallel mode
      )
      
      return(list(file = basename(tif_file), status = "success"))
      
    }, error = function(e) {
      cat("  ERROR processing", basename(tif_file), ":", e$message, "\n")
      return(list(file = basename(tif_file), status = "error", message = e$message))
    })
  }
  
  # Process files
  start_time <- Sys.time()
  
  if(use_parallel) {
    results <- mclapply(tif_files, process_one, mc.cores = n_cores)
  } else {
    results <- lapply(tif_files, process_one)
  }
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Summarize results
  statuses <- sapply(results, function(x) x$status)
  
  cat("\n=== SUMMARY ===\n")
  cat("Total files:", length(tif_files), "\n")
  cat("  Success:", sum(statuses == "success"), "\n")
  cat("  Skipped:", sum(statuses == "skipped"), "\n")
  cat("  Errors:", sum(statuses == "error"), "\n")
  cat("Time elapsed:", round(elapsed, 1), "seconds\n")
  
  if(sum(statuses == "error") > 0) {
    cat("\nFiles with errors:\n")
    error_files <- sapply(results[statuses == "error"], function(x) x$file)
    cat(paste0("  - ", error_files, collapse = "\n"), "\n")
  }
  
  return(invisible(results))
}

# ========================
# RUN BATCH PROCESSING
# ========================

cat("Starting batch multi-contour generation...\n\n")
cat("Configuration:\n")
cat("  Input:", RASTER_DIR, "\n")
cat("  Output:", OUTPUT_DIR, "\n")
cat("  Quantiles:", paste(QUANTILES, collapse = ", "), "\n")
cat("  Parallel:", USE_PARALLEL, "\n")
if(USE_PARALLEL) cat("  Cores:", N_CORES, "\n")
cat("  Overwrite:", OVERWRITE, "\n")
cat("\n")

results <- batch_create_multi_contours(
  raster_dir = RASTER_DIR,
  output_dir = OUTPUT_DIR,
  quantiles = QUANTILES,
  use_parallel = USE_PARALLEL,
  n_cores = N_CORES,
  overwrite = OVERWRITE
)

cat("\n=== COMPLETE ===\n")
cat("Multi-contour shapefiles saved to:", OUTPUT_DIR, "\n\n")

# ========================
# USAGE EXAMPLES
# ========================

# Example 1: Load and filter a multi-contour shapefile
# contours <- st_read("output/shapes/baseline_match_v3/multi/pam_baleen_Sei_Whale_fixed_bw10000_multi.shp")
# q90_only <- contours %>% filter(quantile == 0.90)
# q95_only <- contours %>% filter(quantile == 0.95)

# Example 2: Plot multiple contours
# library(ggplot2)
# ggplot() +
#   geom_sf(data = contours, aes(color = quantile_label), fill = NA) +
#   scale_color_manual(values = c("q85" = "blue", "q90" = "red", "q95" = "green"))

# Example 3: Run for specific parameter folder
# batch_create_multi_contours(
#   raster_dir = "output/tif/my_parameters/",
#   output_dir = "output/shapes/my_parameters/multi/",
#   quantiles = c(0.50, 0.75, 0.90, 0.95),
#   use_parallel = TRUE,
#   overwrite = TRUE
# )
