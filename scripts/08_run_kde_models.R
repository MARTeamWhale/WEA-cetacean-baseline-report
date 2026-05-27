# HEADER --------------------------------------------
#
# Author: Laura Joan Feyrer
# Email:  ljfeyrer@dal.ca
# Most Recent Date Updated: 2026-03-15
#
# Script Name: run_all_kdes_optimized_SIMPLIFIED_v6.R
#
# Description:
## SIMPLIFIED version matching KDES_CE_Collab.R exactly.
## Fixed coastline extent to match CE_Collab plotting.
## Grouped sightings KDEs: beaked whales only, and Deep Divers
## (beaked + Physeter sperm whales).
#
# Changes from previous version (v5):
# - Refactored to remove repeated boilerplate; no functional changes.
# - Added run_kde_group() helper: wraps proximity weight cache logic +
#   performKDE() call into one function, eliminating ~8 copy-paste blocks.
# - Added clean_output_files() helper: replaces 3 identical vapply/grepl
#   cleanup loops (tif, shp, png) with one call each.
# - performKDE() calls now use a base_kde_args list for the 7 shared
#   global arguments (window, sigma, resolution, etc.) so each call
#   only specifies what differs.
# - Result: ~700 lines with identical outputs.
#

# SET OPTIONS ---------------------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 999)
options(encoding = "UTF-8")

# RESOLUTION - MUST MATCH CE_COLLAB EXACTLY
RESOLUTION <- 500  # 500m cells

# ...............
# PERFORMANCE SETTINGS------
# ...............

ENABLE_CACHING  <- FALSE
FAST_MODE       <- FALSE
USE_PARALLEL    <- TRUE
CLEAR_CACHE     <- FALSE
CACHE_DIR       <- "cache/"

# OUTPUT ORGANIZATION
USE_PARAMETER_FOLDERS <- FALSE
PARAMETER_LABEL       <- ""

# ...............
# DATASET AND PARAMETER SELECTION------
# ...............

RUN_SIGHTINGS <- TRUE
RUN_PAM       <- TRUE

species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")

# ANALYSIS TYPE OPTIONS
RUN_INDIVIDUAL_SPECIES          <- TRUE
RUN_FAMILY_GROUPS               <- TRUE
RUN_INDIVIDUAL_PAM              <- TRUE
RUN_GROUPED_PAM                 <- TRUE
RUN_BEAKED_PAM                  <- TRUE
RUN_GROUPED_BEAKED_PAM          <- TRUE
RUN_GROUPED_BEAKED_SIGHTINGS    <- TRUE   # grouped beaked whale sightings KDE
RUN_GROUPED_DEEP_DIVERS_SIGHTINGS <- TRUE # beaked + Physeter sperm whales sightings KDE
RUN_ANY_BEAKED_PAM              <- RUN_BEAKED_PAM || RUN_GROUPED_BEAKED_PAM

# BANDWIDTH OPTIONS
USE_BW_DIGGLE_INDIVIDUAL <- FALSE
USE_BW_DIGGLE_FAMILY     <- FALSE

# PROXIMITY WEIGHTING OPTIONS
USE_PROXIMITY_WEIGHTING_INDIVIDUAL <- TRUE
USE_PROXIMITY_WEIGHTING_FAMILY     <- FALSE

# BANDWIDTH CONFIGURATION
DEFAULT_SIGHTINGS_SIGMA <- 10000
DEFAULT_PAM_SIGMA       <- 10000
SPECIES_BANDWIDTH       <- list()

# PROXIMITY WEIGHTING CONFIGURATION
PROXIMITY_K_NEIGHBORS  <- 8
PROXIMITY_MAX_DISTANCE <- 5000
PROXIMITY_MIN_WEIGHT   <- 0.05

# VISUALIZATION CONFIGURATION
CONTOUR_QUANTILE <- 0.90
buffer_percent   <- 0.05

# ...............
# SPECIES GROUPING DEFINITIONS------
# ...............

# Common names in cetacean_sf treated as beaked whales.
# Adjust to match your data's common_name values.
beaked_species_sightings <- c(
  "Northern Bottlenose Whale",
  "Sowerby's Beaked Whale",
  "True's Beaked Whale",
  "Gervais' Beaked Whale",
  "Blainville's Beaked Whale",
  "Cuvier's Beaked Whale"
)

# Regex applied to common_name (case-insensitive) to identify Physeter sperm whales.
# Catches "Sperm Whale", "sperm whale", "SPERM WHALE", "Physeter macrocephalus", etc.
# Kogia (dwarf/pygmy sperm whales) excluded deliberately.
# Consistent with CE_Collab Deep_divers definition:
#   "Physeter|SPERM\\s+WHALE|Mesoplodon|Hyperoodon|Ziphius|..."
PHYSETER_REGEX <- "(?i)(physeter|sperm\\s+whale)"

# PAM buffer: stations within this distance of study area boundary are included
PAM_BUFFER_DIST_M <- 10000  # 10 km; matches PAM summaries and comparison maps

# Sightings date filter
DATE_FILTER_SIGHTINGS <- as.Date("2010-01-01")

# Minimum sightings records per species to include in individual KDE
MIN_RECORDS_SPECIES <- 25

# ................................
# LIBRARIES------
# ................................

pacman::p_load(sf, tidyverse, readxl, here, leaflet,
               scales, terra, ggrepel, viridis, ggspatial, spatstat,
               dplyr, patchwork, digest, lubridate)
suppressWarnings(source(here::here("scripts/00_load_helpers.R")))

if (USE_PARALLEL) {
  pacman::p_load(parallel)
  detected_cores <- parallel::detectCores()
  n_cores <- if (is.na(detected_cores)) 1 else max(1, detected_cores - 1)
  cat("Parallel processing enabled with", n_cores, "cores\n")
}

# Projections
UTM20 <- SPATIAL_CRS_UTM20

# ................................
# OUTPUT DIRECTORIES------
# ................................

if (USE_PARAMETER_FOLDERS) {
  param_suffix <- if (PARAMETER_LABEL != "") PARAMETER_LABEL else
    paste0("bw", DEFAULT_PAM_SIGMA[1], "_q",
           gsub("\\.", "", as.character(CONTOUR_QUANTILE[1])))
  raster_dir    <- file.path("output/tif",    param_suffix, "")
  plot_dir      <- file.path("output/figs/KDE_maps", param_suffix, "")
  shapefile_dir <- file.path("output/shapes", param_suffix, "")
  cat("\nUsing parameter folder:", param_suffix, "\n")
} else {
  raster_dir    <- "output/tif/"
  plot_dir      <- "output/figs/KDE_maps/"
  shapefile_dir <- "output/shapes/"
}

if (ENABLE_CACHING && !dir.exists(CACHE_DIR))
  dir.create(CACHE_DIR, recursive = TRUE)

# ................................
# CACHING FUNCTIONS------
# ................................

save_to_cache <- function(obj, name, subdir = NULL) {
  if (!ENABLE_CACHING) return(invisible(NULL))
  cache_path <- if (!is.null(subdir)) file.path(CACHE_DIR, subdir) else CACHE_DIR
  if (!dir.exists(cache_path)) dir.create(cache_path, recursive = TRUE)
  saveRDS(obj, file.path(cache_path, paste0(name, ".rds")))
}

load_from_cache <- function(name, subdir = NULL) {
  if (!ENABLE_CACHING) return(NULL)
  cache_path <- if (!is.null(subdir)) file.path(CACHE_DIR, subdir) else CACHE_DIR
  f <- file.path(cache_path, paste0(name, ".rds"))
  if (file.exists(f)) { cat("  Loading from cache:", f, "\n"); return(readRDS(f)) }
  return(NULL)
}

get_cache_hash <- function(...) digest::digest(list(...), algo = "md5")

# ................................
# CLEANUP HELPER------
# ................................

# Removes files in dir matching any of the given prefixes against basename(file).
# pattern: glob passed to list.files (e.g. "*.tif")
clean_output_files <- function(dir, pattern, prefixes) {
  if (!dir.exists(dir)) { dir.create(dir, recursive = TRUE); return(invisible(NULL)) }
  files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0 || length(prefixes) == 0) return(invisible(NULL))
  to_rm <- files[vapply(files, function(f)
    any(vapply(prefixes, function(p) grepl(p, basename(f)), logical(1))),
    logical(1))]
  if (length(to_rm) > 0) {
    file.remove(to_rm)
    cat("Removed", length(to_rm), "files matching enabled analyses from", dir, "\n")
  }
}

# ................................
# CLEANUP------
# ................................

if (!FAST_MODE) {
  cat("\n=== CLEANING OUTPUT FOLDERS ===\n\n")
  
  pfx <- c()
  if (RUN_SIGHTINGS) {
    if (RUN_INDIVIDUAL_SPECIES)           pfx <- c(pfx, "sightings_species")
    if (RUN_FAMILY_GROUPS)                pfx <- c(pfx, "sightings_Odontocetes", "sightings_Mysticetes")
    if (RUN_GROUPED_BEAKED_SIGHTINGS)     pfx <- c(pfx, "sightings_beaked_grouped")
    if (RUN_GROUPED_DEEP_DIVERS_SIGHTINGS) pfx <- c(pfx, "sightings_deep_divers_grouped")
  }
  if (RUN_PAM) {
    if (RUN_INDIVIDUAL_PAM) pfx <- c(pfx, "pam_baleen_")
    if (RUN_GROUPED_PAM)    pfx <- c(pfx, "pam_grouped_baleen", "grouped_baleen_pam")
  }
  if (RUN_BEAKED_PAM)         pfx <- c(pfx, "pam_beaked_")
  if (RUN_GROUPED_BEAKED_PAM) pfx <- c(pfx, "pam_grouped_beaked", "grouped_beaked_pam")
  
  if (!ENABLE_CACHING || CLEAR_CACHE) {
    clean_output_files(raster_dir[1],    "\\.tif$",                       pfx)
    clean_output_files(shapefile_dir[1], "\\.(shp|shx|dbf|prj|cpg)$",    pfx)
  } else {
    cat("Keeping cached raster files\n")
    if (!dir.exists(raster_dir[1])) dir.create(raster_dir[1], recursive = TRUE)
  }
  clean_output_files(plot_dir[1], "\\.png$", pfx)
  
  if (CLEAR_CACHE && dir.exists(CACHE_DIR)) {
    cache_files <- list.files(CACHE_DIR, full.names = TRUE, recursive = TRUE)
    if (length(cache_files) > 0) { file.remove(cache_files); cat("Cleared cache\n") }
  }
}

# ................................
# CALCULATE PROXIMITY WEIGHTS------
# ................................

calculate_proximity_weights <- function(data_sf,
                                        k_neighbors  = PROXIMITY_K_NEIGHBORS,
                                        max_distance = PROXIMITY_MAX_DISTANCE,
                                        min_weight   = PROXIMITY_MIN_WEIGHT) {
  cat("  Calculating proximity weights: k =", k_neighbors,
      "| max_dist =", max_distance, "m | min_weight =", min_weight, "\n")
  coords   <- st_coordinates(data_sf)
  dist_mat <- as.matrix(dist(coords))
  prox     <- apply(dist_mat, 1, function(row) {
    mean(sort(row)[2:min(k_neighbors + 1, length(row))])
  })
  w <- 1 / (1 + (prox / max_distance))
  w <- (w - min(w)) / (max(w) - min(w))
  w <- pmax(w, min_weight)
  cat("  Weights - min:", round(min(w), 3), "| mean:", round(mean(w), 3),
      "| max:", round(max(w), 3), "\n")
  return(w)
}

# ................................
# RUN KDE GROUP HELPER------
# ................................

# Handles proximity weight caching and calls performKDE().
# use_prox:    logical - whether to apply proximity weighting
# cache_label: short string for the weights cache key
# species_col / species_list: passed through to performKDE()
# group_label: data_source_label for plot titles
# prefix:      output_prefix for performKDE(); "_proximity" appended automatically if use_prox=TRUE
# sigma / bw_diggle: bandwidth args
# base_args:   list of shared performKDE() args (window, resolution, coastline, etc.)

run_kde_group <- function(data_sf, species_col, species_list,
                          use_prox, cache_label,
                          prefix, group_label,
                          sigma, bw_diggle = FALSE,
                          base_args) {
  
  weight_col   <- NULL
  output_prefix <- prefix
  
  if (use_prox) {
    wt_key <- paste0(cache_label, "_weights_",
                     get_cache_hash(nrow(data_sf), PROXIMITY_K_NEIGHBORS,
                                    PROXIMITY_MAX_DISTANCE, PROXIMITY_MIN_WEIGHT))
    cached_w <- load_from_cache(wt_key, "weights")
    if (!is.null(cached_w)) {
      data_sf$proximity_weight <- cached_w
    } else {
      data_sf$proximity_weight <- calculate_proximity_weights(data_sf)
      save_to_cache(data_sf$proximity_weight, wt_key, "weights")
    }
    weight_col    <- "proximity_weight"
    output_prefix <- paste0(prefix, "_proximity")
  }
  
  do.call(performKDE, c(
    list(data_sf          = data_sf,
         species_col      = species_col,
         species_list     = species_list,
         weight_col       = weight_col,
         output_prefix    = output_prefix,
         sigma_val        = sigma,
         use_bw_diggle    = bw_diggle,
         data_source_label = group_label),
    base_args
  ))
}

# ................................
# SIMPLIFIED CONTOUR CREATION------
# ................................

create_contour_90_baseline <- function(kde_raster, output_path) {
  kde_raster[kde_raster <= 0] <- NA
  threshold_90 <- quantile(values(kde_raster, mat = FALSE, na.rm = TRUE),
                           probs = 0.90, type = 6, na.rm = TRUE)
  cat("    Threshold value:", threshold_90, "\n")
  kde_binary <- kde_raster
  kde_binary[kde_binary < threshold_90] <- NA
  kde_binary[!is.na(kde_binary)] <- 1
  kde_sf <- st_as_sf(as.polygons(kde_binary, dissolve = TRUE)) %>%
    st_make_valid() %>%
    mutate(contour = "0.90") %>%
    select(contour, geometry)
  cat("    Created", nrow(kde_sf), "polygon(s)\n")
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  return(kde_sf)
}

# ................................
# UNIFIED KDE FUNCTION------
# ................................

performKDE <- function(data_sf,
                       species_col,
                       species_list,
                       weight_col       = NULL,
                       buffer_percent   = 0.05,
                       sigma_val        = 10000,
                       window           = NULL,
                       output_prefix,
                       threshold_quantile       = NULL,
                       output_dir               = raster_dir,
                       plot_output_dir          = plot_dir,
                       shapefile_output_dir     = shapefile_dir,
                       use_bw_diggle            = FALSE,
                       data_source_label        = "KDE",
                       resolution               = RESOLUTION,
                       coastline_for_plot       = NULL,
                       owa_for_plot             = NULL,
                       study_area_for_plot      = NULL) {
  
  if (!dir.exists(output_dir))      dir.create(output_dir,      recursive = TRUE)
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
  
  shapefile_crs <- st_crs(data_sf)$wkt
  raster_list   <- list()
  sigma_store   <- list()
  global_min    <- Inf
  global_max    <- -Inf
  
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("'", "", name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    gsub(" ", "_", name)
  }
  
  cache_key  <- get_cache_hash(species_list, weight_col, sigma_val, use_bw_diggle,
                               output_prefix, resolution, buffer_percent,
                               if (!is.null(weight_col))
                                 list(PROXIMITY_K_NEIGHBORS, PROXIMITY_MAX_DISTANCE,
                                      PROXIMITY_MIN_WEIGHT) else NULL)
  cache_name <- paste0(output_prefix, "_", cache_key)
  
  #.............................
  # FIRST PASS: Create or Load KDE rasters-----
  #.............................
  
  if (FAST_MODE) {
    cat("\n=== FAST MODE: Loading existing rasters ===\n")
    for (species in species_list) {
      species_clean <- sanitize_filename(species)
      bw_label <- ifelse(use_bw_diggle, "diggle_bw*",
                         paste0("fixed_bw", round(sigma_val, 0)))
      tif_files <- list.files(output_dir,
                              pattern = glob2rx(paste0("KDE_", output_prefix, "_",
                                                       species_clean, "_", bw_label, ".tif")),
                              full.names = TRUE)
      if (length(tif_files) > 0) {
        raster_kd  <- rast(tif_files[1])
        bw_match   <- regmatches(basename(tif_files[1]),
                                 regexpr("bw[0-9]+", basename(tif_files[1])))
        sigma_used <- as.numeric(gsub("bw", "", bw_match))
        bw_method  <- ifelse(grepl("diggle", basename(tif_files[1])), "diggle", "fixed")
        raster_list[[species]] <- list(raster = raster_kd, sigma = sigma_used,
                                       method = bw_method,
                                       bw_label = paste0(bw_method, "_bw", sigma_used),
                                       raster_path = tif_files[1])
        global_min <- min(global_min, min(values(raster_kd), na.rm = TRUE))
        global_max <- max(global_max, max(values(raster_kd), na.rm = TRUE))
        cat("  Loaded raster for", species, "\n")
      } else {
        cat("  Warning: No raster found for", species, "\n")
      }
    }
    
  } else {
    cached_metadata <- load_from_cache(cache_name, "rasters")
    
    if (!is.null(cached_metadata)) {
      cat("\n=== Using cached KDE results ===\n")
      for (species in names(cached_metadata$raster_list)) {
        rp <- cached_metadata$raster_list[[species]]$raster_path
        if (file.exists(rp)) {
          raster_list[[species]] <- cached_metadata$raster_list[[species]]
          raster_list[[species]]$raster <- rast(rp)
          cat("  Loaded", species, "from", basename(rp), "\n")
        } else {
          cat("  Cached raster missing for", species, "- will recalculate\n")
          cached_metadata <- NULL; break
        }
      }
      if (!is.null(cached_metadata)) {
        sigma_store <- cached_metadata$sigma_store
        global_min  <- cached_metadata$global_min
        global_max  <- cached_metadata$global_max
      }
    }
    
    if (is.null(cached_metadata)) {
      cat("\n=== Calculating KDE rasters ===\n")
      for (species in species_list) {
        current_sf <- data_sf[data_sf[[species_col]] == species, ]
        if (nrow(current_sf) == 0) next
        
        cat("Processing", species, "(", nrow(current_sf), "records)...\n")
        tryCatch({
          coords <- st_coordinates(current_sf)
          
          sp_window <- if (is.null(window)) {
            xr <- range(coords[, "X"]); yr <- range(coords[, "Y"])
            xb <- diff(xr) * buffer_percent; yb <- diff(yr) * buffer_percent
            owin(xrange = c(xr[1] - xb, xr[2] + xb),
                 yrange = c(yr[1] - yb, yr[2] + yb))
          } else window
          
          if (!is.null(weight_col)) {
            pts <- ppp(coords[,"X"], coords[,"Y"], window = sp_window,
                       marks = data.frame(weights = current_sf[[weight_col]] * 10000))
          } else {
            pts <- ppp(coords[,"X"], coords[,"Y"], window = sp_window)
          }
          pts <- rjitter(pts, retry = TRUE, nsim = 1, drop = TRUE)
          
          sigma_used <- if (isTRUE(use_bw_diggle)) {
            bw_method <- "diggle"; bw.diggle(pts)
          } else {
            bw_method <- "fixed"
            if (exists("SPECIES_BANDWIDTH") && species %in% names(SPECIES_BANDWIDTH))
              SPECIES_BANDWIDTH[[species]] else sigma_val
          }
          sigma_store[[species]] <- list(sigma = sigma_used, method = bw_method)
          
          dimx <- round(diff(sp_window$xrange) / resolution)
          dimy <- round(diff(sp_window$yrange) / resolution)
          cat("  Resolution:", resolution, "m (grid:", dimx, "x", dimy, ")\n")
          
          kd <- density.ppp(pts, sigma = sigma_used, positive = TRUE,
                            kernel = "gaussian",
                            weights = if (!is.null(weight_col)) marks(pts) else NULL,
                            dimyx = c(dimy, dimx), diggle = TRUE)
          
          raster_kd <- rast(kd)
          crs(raster_kd) <- shapefile_crs
          
          max_val <- max(values(raster_kd), na.rm = TRUE)
          if (!is.na(max_val) && max_val > 0) {
            raster_kd <- raster_kd / max_val
            if (!is.null(threshold_quantile)) {
              thr <- quantile(values(raster_kd, mat = FALSE, na.rm = TRUE),
                              threshold_quantile, na.rm = TRUE)
              raster_kd[raster_kd < thr] <- NA
            }
          }
          cat("  Actual resolution:", round(mean(res(raster_kd)), 2), "m\n")
          
          global_min <- min(global_min, min(values(raster_kd), na.rm = TRUE))
          global_max <- max(global_max, max(values(raster_kd), na.rm = TRUE))
          
          bw_label      <- paste0(bw_method, "_bw", round(sigma_used, 0))
          raster_fname  <- paste0(output_dir, "KDE_", output_prefix, "_",
                                  sanitize_filename(species), "_", bw_label, ".tif")
          writeRaster(raster_kd, filename = raster_fname, overwrite = TRUE)
          
          raster_list[[species]] <- list(raster = raster_kd, sigma = sigma_used,
                                         method = bw_method, bw_label = bw_label,
                                         raster_path = raster_fname)
        }, error = function(e) cat("  Error:", species, "-", e$message, "\n"))
      }
      
      save_to_cache(list(
        raster_list = lapply(raster_list, function(x)
          x[c("sigma","method","bw_label","raster_path")]),
        sigma_store = sigma_store,
        global_min  = global_min,
        global_max  = global_max
      ), cache_name, "rasters")
    }
  }
  
  #.............................
  # GENERATE QUANTILE SHAPEFILES------
  #.............................
  
  cat("\n=== GENERATING Q90 CONTOUR SHAPEFILES ===\n")
  q90_dir <- file.path(shapefile_output_dir, "q90")
  if (!dir.exists(q90_dir)) dir.create(q90_dir, recursive = TRUE)
  
  for (species in names(raster_list)) {
    tryCatch({
      base_name     <- tools::file_path_sans_ext(basename(raster_list[[species]]$raster_path))
      species_clean <- str_extract(base_name, "(?<=KDE_)[^_]+_.*")
      contour_path  <- file.path(q90_dir, paste0(species_clean, "_contour90.shp"))
      cat("  Creating contour for", species, "\n")
      contour_sf    <- create_contour_90_baseline(rast(raster_list[[species]]$raster_path),
                                                  contour_path)
      raster_list[[species]]$shapefile_path <- contour_path
      raster_list[[species]]$contour_sf     <- contour_sf
    }, error = function(e) cat("  Error:", species, "-", conditionMessage(e), "\n"))
  }
  cat("\n+ Created", length(raster_list), "baseline contour shapefiles\n\n")
  
  #.............................
  # CREATE PLOTS------
  #.............................
  
  readable_names <- setNames(as.list(species_list), species_list)
  code_map <- c(Bb="Sei Whale", Bm="Blue Whale", Bp="Fin Whale",
                Mn="Humpback Whale", Ba="Minke Whale",
                Eg="North Atlantic Right Whale",
                Ha="Northern Bottlenose Whale",
                Mb="Sowerby's Beaked Whale", Zc="Goose Beaked Whale",
                MmMe="True's/Gervais' Beaked Whale")
  for (sp in species_list)
    if (sp %in% names(code_map)) readable_names[[sp]] <- code_map[[sp]]
  
  xlims <- ylims <- NULL
  if (!is.null(study_area_for_plot)) {
    plot_extent <- st_as_sfc(st_bbox(study_area_for_plot))
    if (nrow(data_sf) > 0) {
      plot_extent <- st_union(plot_extent, st_as_sfc(st_bbox(data_sf)))
    }
    bb <- st_bbox(st_buffer(plot_extent, dist = 5000))
    xlims <- c(bb["xmin"], bb["xmax"])
    ylims <- c(bb["ymin"], bb["ymax"])
  }
  
  plot_list <- list()
  cat("=== Generating plots ===\n")
  
  for (species in species_list) {
    rd <- raster_list[[species]]
    if (is.null(rd)) next
    cat("  Plotting", species, "...\n")
    
    kde_df <- as.data.frame(rd$raster, xy = TRUE) %>%
      filter(!is.na(lyr.1)) %>%
      mutate(lyr.1_sqrt = sqrt(lyr.1))
    
    sp_pts  <- data_sf[data_sf[[species_col]] == species, ]
    is_pam  <- any(c("site","deployment") %in% names(sp_pts))
    bw_text <- paste0(toupper(rd$method), " BW: ", round(rd$sigma, 0), "m")
    ct_text <- paste0(CONTOUR_QUANTILE * 100, "% contour (baseline)")
    
    p <- ggplot() +
      geom_tile(data = kde_df, aes(x = x, y = y, fill = lyr.1_sqrt)) +
      scale_fill_viridis_c(option = "viridis", name = "sqrt(Density)",
                           limits = c(0, max(kde_df$lyr.1_sqrt)),
                           begin = 0.15, end = 1.0)
    
    if (!is.null(rd$contour_sf) && nrow(rd$contour_sf) > 0)
      p <- p + map_layer_kde_contour(rd$contour_sf)
    if (!is.null(coastline_for_plot))
      p <- p + geom_sf(data = coastline_for_plot, fill = "grey90",
                       color = "grey70", linewidth = 0.3)
    if (!is.null(owa_for_plot))
      p <- p + map_layer_wea(owa_for_plot)
    if (nrow(sp_pts) > 0) {
      p <- p + geom_sf(data = sp_pts, color = "white",
                       size = ifelse(is_pam, 1.5, 0.8),
                       shape = 21, fill = "black",
                       stroke = ifelse(is_pam, 0.5, 0.3),
                       alpha = ifelse(is_pam, 0.8, 0.6))
    }
    if (!is.null(study_area_for_plot))
      p <- p + map_layer_study_area(study_area_for_plot)
    
    coord_args <- list(crs = st_crs(UTM20), expand = FALSE)
    if (!is.null(xlims)) { coord_args$xlim <- xlims; coord_args$ylim <- ylims }
    p <- p + do.call(coord_sf, coord_args) +
      map_base_theme() +
      labs(title = paste0(data_source_label, ": ", readable_names[[species]])) +
      theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
            panel.grid = element_line(color = "grey95"))
    
    if (!is.null(xlims))
      p <- p +
      map_annotation_label(xlims[1], ylims[2], bw_text, vjust = 1.3) +
      map_annotation_label(xlims[1], ylims[2], ct_text, vjust = 2.8)
    
    plot_list[[species]] <- p
  }
  
  # Save plots
  plot_combined <- NULL
  if (length(plot_list) > 0) {
    is_grouped <- grepl("grouped|combined_Mysticetes|combined_Odontocetes|guild",
                        output_prefix, ignore.case = TRUE)
    if (is_grouped && length(plot_list) > 1) {
      plot_combined <- wrap_plots(plot_list, guides = "collect")
      ggsave(paste0(plot_output_dir, "combined_KDE_", output_prefix, "_plots.png"),
             plot_combined, width = 16, height = 12, dpi = 300)
    }
    for (sp in names(plot_list)) {
      if (!is.null(plot_list[[sp]]))
        ggsave(paste0(plot_output_dir, "KDE_", output_prefix, "_",
                      sanitize_filename(sp), "_", raster_list[[sp]]$bw_label, ".png"),
               plot_list[[sp]], width = 8, height = 6, dpi = 300)
    }
  }
  
  return(list(plot = plot_combined, bandwidths = sigma_store,
              global_range = c(global_min, global_max)))
}

#.............................
# LOAD STUDY AREA (WITH CACHING)-------
#.............................

cat("\n=== LOADING STUDY AREA ===\n\n")

study_area    <- load_from_cache("study_area", "spatial")

if (is.null(study_area)) {
  study_area <- load_spatial_layer("study_area", crs = 32620)
  save_to_cache(study_area, "study_area", "spatial")
}

study_area_pam_buffer <- st_buffer(study_area, dist = PAM_BUFFER_DIST_M)

common_window_cache_name <- paste0("common_window_pam_buffer_", PAM_BUFFER_DIST_M, "m")
common_window <- load_from_cache(common_window_cache_name, "spatial")

if (is.null(common_window)) {
  bb <- st_bbox(study_area_pam_buffer)
  xb <- diff(c(bb["xmin"], bb["xmax"])) * buffer_percent
  yb <- diff(c(bb["ymin"], bb["ymax"])) * buffer_percent
  common_window <- owin(xrange = c(bb["xmin"] - xb, bb["xmax"] + xb),
                        yrange = c(bb["ymin"] - yb, bb["ymax"] + yb))
  save_to_cache(common_window, common_window_cache_name, "spatial")
  cat("Common window - X:", common_window$xrange, "Y:", common_window$yrange, "\n")
}
cat("Study area loaded | PAM inclusion buffer:", PAM_BUFFER_DIST_M / 1000,
    "km | Window buffer:", buffer_percent * 100, "%\n\n")

#.............................
# LOAD BASE LAYERS FOR PLOTTING (WITH CACHING)-------
#.............................

cat("=== LOADING BASE LAYERS FOR PLOTTING ===\n\n")

# Crop coastline to study area (matching CE_Collab)
coastline <- load_from_cache("coastline_cropped", "spatial")
if (is.null(coastline)) {
  coastline <- load_spatial_layer("land", crs = 32620, required = FALSE)
  if (!is.null(coastline)) {
    coastline <- coastline %>% st_crop(st_bbox(study_area_pam_buffer))
    save_to_cache(coastline, "coastline_cropped", "spatial")
    cat("Loaded and cropped coastline\n")
  } else { coastline <- NULL; cat("Coastline file not found\n") }
} else cat("Loaded coastline from cache\n")

# Crop OWA to study area (matching CE_Collab)
owa <- load_from_cache("owa_cropped", "spatial")
if (is.null(owa)) {
  owa <- load_spatial_layer("wea", crs = 32620, required = FALSE)
  if (!is.null(owa)) {
    owa <- owa %>% st_crop(st_bbox(study_area_pam_buffer))
    save_to_cache(owa, "owa_cropped", "spatial")
    cat("Loaded and cropped OWA\n")
  } else { owa <- NULL; cat("OWA file not found\n") }
} else cat("Loaded OWA from cache\n")

# Shared args passed to every performKDE() call via run_kde_group()
base_sightings_args <- list(
  window             = common_window,
  resolution         = RESOLUTION,
  buffer_percent     = 0.05,
  coastline_for_plot = coastline,
  owa_for_plot       = owa,
  study_area_for_plot = study_area
)
base_pam_args <- base_sightings_args  # identical for now; split if they diverge

cat("\n")

#.............................
# PROCESS COMBINED SIGHTINGS DATA (WITH CACHING)----
#.............................

if (RUN_SIGHTINGS) {
  
  cat("\n=== PROCESSING COMBINED SIGHTINGS DATA ===\n\n")
  
  cetacean_sf <- load_from_cache("cetacean_sf_date", "sightings")
  species_min <- load_from_cache("species_min",      "sightings")
  
  if (is.null(cetacean_sf) || is.null(species_min)) {
    combined_data <- read_csv("output/data/combined_dedup_1km_day.csv",
                              show_col_types = FALSE) %>%
      filter(as.Date(date_utc) >= DATE_FILTER_SIGHTINGS)
    cat("Records after date filter:", nrow(combined_data), "\n")
    
    cetacean_sf <- st_as_sf(combined_data, coords = c("lon","lat"), crs = 4326) %>%
      st_transform(UTM20) %>% st_make_valid() %>%
      filter(family %in% c("Odontocete","Mysticete"))
    cat("Cetacean records:", nrow(cetacean_sf), "\n")
    
    species_min <- cetacean_sf %>% st_drop_geometry() %>%
      count(common_name) %>%
      filter(n >= MIN_RECORDS_SPECIES) %>%
      pull(common_name)
    
    save_to_cache(cetacean_sf, "cetacean_sf_date", "sightings")
    save_to_cache(species_min, "species_min",      "sightings")
  } else cat("Loaded sightings from cache\n")
  
  cat(length(species_min), "species with >=", MIN_RECORDS_SPECIES, "records\n")
  
  #.............................
  # INDIVIDUAL SPECIES KDEs----
  #.............................
  
  if (RUN_INDIVIDUAL_SPECIES && length(species_min) > 0) {
    cat("\n=== INDIVIDUAL SPECIES KDEs ===\n\n")
    species_sf <- cetacean_sf %>% filter(common_name %in% species_min)
    run_kde_group(species_sf, "common_name", species_min,
                  use_prox     = USE_PROXIMITY_WEIGHTING_INDIVIDUAL,
                  cache_label  = "indiv_species",
                  prefix       = "sightings_species",
                  group_label  = "Sightings",
                  sigma        = DEFAULT_SIGHTINGS_SIGMA,
                  bw_diggle    = USE_BW_DIGGLE_INDIVIDUAL,
                  base_args    = base_sightings_args)
  }
  
  #.............................
  # FAMILY-LEVEL KDEs----
  #.............................
  
  if (RUN_FAMILY_GROUPS) {
    cat("\n=== FAMILY-LEVEL KDEs ===\n\n")
    
    family_groups <- list(
      list(filter_expr = quote(family == "Odontocete"),
           label       = "Odontocetes (All Toothed Whales)",
           col         = "family_group",
           cache_lbl   = "odonto",
           prefix      = "sightings_Odontocetes"),
      list(filter_expr = quote(family == "Mysticete"),
           label       = "Mysticetes (All Baleen Whales)",
           col         = "family_group",
           cache_lbl   = "mysti",
           prefix      = "sightings_Mysticetes")
    )
    
    for (fg in family_groups) {
      grp_sf <- cetacean_sf %>%
        filter(eval(fg$filter_expr)) %>%
        mutate(family_group = fg$label)
      if (nrow(grp_sf) == 0) next
      run_kde_group(grp_sf, fg$col, fg$label,
                    use_prox    = USE_PROXIMITY_WEIGHTING_FAMILY,
                    cache_label = fg$cache_lbl,
                    prefix      = fg$prefix,
                    group_label = "Sightings",
                    sigma       = DEFAULT_SIGHTINGS_SIGMA,
                    bw_diggle   = USE_BW_DIGGLE_FAMILY,
                    base_args   = base_sightings_args)
    }
  }
  
  #.............................
  # GROUPED BEAKED WHALE SIGHTINGS KDE----
  # Pools all beaked whale sightings (common_name %in% beaked_species_sightings)
  # into a single group and runs one KDE.
  # Edit beaked_species_sightings at the top of the script to match your data.
  #.............................
  
  if (RUN_GROUPED_BEAKED_SIGHTINGS) {
    cat("\n=== GROUPED BEAKED WHALE SIGHTINGS KDE ===\n\n")
    beaked_sf <- cetacean_sf %>%
      filter(common_name %in% beaked_species_sightings) %>%
      mutate(beaked_group = "All Beaked Whales (Sightings)")
    cat("Beaked records:", nrow(beaked_sf), "| Species:",
        paste(unique(beaked_sf$common_name), collapse = ", "), "\n\n")
    if (nrow(beaked_sf) > 0) {
      run_kde_group(beaked_sf, "beaked_group", "All Beaked Whales (Sightings)",
                    use_prox    = USE_PROXIMITY_WEIGHTING_FAMILY,
                    cache_label = "beaked_sightings",
                    prefix      = "sightings_beaked_grouped",
                    group_label = "Sightings",
                    sigma       = DEFAULT_SIGHTINGS_SIGMA,
                    bw_diggle   = USE_BW_DIGGLE_FAMILY,
                    base_args   = base_sightings_args)
    } else {
      cat("  Warning: No beaked sightings matched beaked_species_sightings.\n")
    }
  }
  
  #.............................
  # GROUPED DEEP DIVERS SIGHTINGS KDE----
  # Deep Divers = beaked whales (beaked_species_sightings) + Physeter sperm whales
  # (matched via PHYSETER_REGEX). Consistent with CE_Collab Deep_divers grouping.
  # Edit PHYSETER_REGEX and beaked_species_sightings at top of script as needed.
  #.............................
  
  if (RUN_GROUPED_DEEP_DIVERS_SIGHTINGS) {
    cat("\n=== GROUPED DEEP DIVERS SIGHTINGS KDE (Beaked + Sperm Whales) ===\n\n")
    dd_sf <- cetacean_sf %>%
      filter(grepl(PHYSETER_REGEX, common_name, perl = TRUE) |
               common_name %in% beaked_species_sightings) %>%
      mutate(deep_divers_group = "Deep Divers - Beaked and Sperm Whales (Sightings)")
    cat("Deep Divers records:", nrow(dd_sf), "| Species:",
        paste(sort(unique(dd_sf$common_name)), collapse = ", "), "\n\n")
    if (nrow(dd_sf) > 0) {
      run_kde_group(dd_sf, "deep_divers_group",
                    "Deep Divers - Beaked and Sperm Whales (Sightings)",
                    use_prox    = USE_PROXIMITY_WEIGHTING_FAMILY,
                    cache_label = "deep_divers_sightings",
                    prefix      = "sightings_deep_divers_grouped",
                    group_label = "Sightings",
                    sigma       = DEFAULT_SIGHTINGS_SIGMA,
                    bw_diggle   = USE_BW_DIGGLE_FAMILY,
                    base_args   = base_sightings_args)
    } else {
      cat("  Warning: No Deep Diver sightings found.\n",
          " Check beaked_species_sightings and PHYSETER_REGEX vs cetacean_sf common_name.\n")
    }
  }
  
} # end RUN_SIGHTINGS

#.............................
# PROCESS PAM DATA (WITH CACHING)---------
#.............................

if (RUN_PAM) {
  cat("\n=== PROCESSING PAM DATA ===\n\n")
  
  filepath_raw       <- "input/raw_data/PAM/baleen_presence_laura_2025.csv"
  filepath_processed <- "input/raw_data/PAM/CAM_baleen_presence_days_laura_2025.csv"
  
  if (file.exists(filepath_raw)) {
    cat("Using raw PAM file (matching CE_Collab format)\n")
    baleen_PA <- read.csv(filepath_raw) %>%
      mutate(rec_date = as.Date(rec_date)) %>%
      group_by(site, latitude, longitude, species) %>%
      summarise(detection_days  = sum(presence, na.rm = TRUE),
                effort_days     = n_distinct(rec_date),
                proportion_det  = detection_days / effort_days,
                .groups = "drop") %>%
      filter(effort_days > 0)
    cat("Station-species combinations:", nrow(baleen_PA), "\n")
  } else if (file.exists(filepath_processed)) {
    cat("Warning: Using processed PAM file (may differ from CE_Collab)\n")
    baleen_PA <- read.csv(filepath_processed)
  } else stop("No PAM file found!")
  
  baleen_sf <- baleen_PA %>%
    st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
    st_transform(UTM20) %>%
    st_filter(study_area_pam_buffer)
  
  cat("PAM records:", nrow(baleen_sf), "| Unique sites:",
      length(unique(baleen_sf$site)), "\n\n")
  
  #.............................
  # INDIVIDUAL SPECIES PAM KDEs----
  #.............................
  
  if (RUN_INDIVIDUAL_PAM) {
    cat("=== INDIVIDUAL SPECIES PAM KDEs ===\n\n")
    pam_names <- c(Bb="Sei_Whale", Bm="Blue_Whale", Bp="Fin_Whale",
                   Mn="Humpback_Whale", Eg="North_Atlantic_Right_Whale", Ba="Minke_Whale")
    for (sp in species_list_pam) {
      sp_data <- baleen_sf %>% filter(species == sp)
      if (nrow(sp_data) == 0) next
      cat("Processing", pam_names[sp], "(n =", nrow(sp_data), ")\n")
      do.call(performKDE, c(
        list(data_sf = sp_data, species_col = "species", species_list = sp,
             weight_col = "proportion_det",
             output_prefix = paste0("pam_baleen_", pam_names[sp]),
             sigma_val = DEFAULT_PAM_SIGMA, use_bw_diggle = FALSE,
             data_source_label = "Baleen PAM"),
        base_pam_args))
    }
  }
  
  #.............................
  # GROUPED PAM----
  #.............................
  
  if (RUN_GROUPED_PAM) {
    cat("\n=== GROUPED PAM KDE ===\n\n")
    baleen_grp <- baleen_sf %>%
      group_by(geometry) %>%
      summarise(detection_days = sum(detection_days, na.rm = TRUE),
                effort_days    = dplyr::first(effort_days), .groups = "drop") %>%
      mutate(proportion_det = detection_days / effort_days,
             baleen_group   = "All Baleen Whales (PAM)")
    do.call(performKDE, c(
      list(data_sf = baleen_grp, species_col = "baleen_group",
           species_list = "All Baleen Whales (PAM)",
           weight_col = "proportion_det",
           output_prefix = "pam_grouped_baleen",
           sigma_val = DEFAULT_PAM_SIGMA, use_bw_diggle = FALSE,
           data_source_label = "PAM"),
      base_pam_args))
  }
}

#.............................
# PROCESS BEAKED WHALE PAM DATA---------
#.............................

if (RUN_ANY_BEAKED_PAM) {
  cat("\n=== PROCESSING BEAKED WHALE PAM DATA ===\n\n")
  
  beaked_cache_name <- paste0("beaked_sf_pam_buffer_", PAM_BUFFER_DIST_M, "m")
  beaked_sf <- load_from_cache(beaked_cache_name, "pam")
  if (is.null(beaked_sf)) {
    beaked_names <- c(Ha="Northern_Bottlenose_Whale", Mb="Sowerbys_Beaked_Whale",
                      Zc="Goose_Beaked_Whale", MmMe="Trues_Gervais_Beaked_Whale")
    beaked_raw <- read_csv("input/raw_data/PAM/beaked_pam_results_2026-01-12.csv",
                           show_col_types = FALSE)
    beaked_station_effort <- beaked_raw %>%
      group_by(deployment, station, latitude, longitude) %>%
      summarise(effort_days = n_distinct(rec_date), .groups = "drop")
    beaked_observed <- beaked_raw %>%
      group_by(deployment, station, latitude, longitude, species) %>%
      summarise(detection_days = sum(presence, na.rm = TRUE), .groups = "drop")
    beaked_sf <- tidyr::crossing(
        beaked_station_effort,
        species = names(beaked_names)
      ) %>%
      left_join(
        beaked_observed,
        by = c("deployment", "station", "latitude", "longitude", "species")
      ) %>%
      mutate(
        detection_days = coalesce(detection_days, 0),
        proportion_det = detection_days / effort_days
      ) %>%
      st_as_sf(coords = c("longitude","latitude"), crs = 4326) %>%
      st_transform(UTM20) %>%
      st_filter(study_area_pam_buffer)
    save_to_cache(beaked_sf, beaked_cache_name, "pam")
    cat("Cached beaked PAM data within", PAM_BUFFER_DIST_M / 1000,
        "km of study area\n")
  } else cat("Loaded beaked PAM data from cache\n")
  
  cat("Beaked PAM records:", nrow(beaked_sf), "\n\n")
  
  beaked_names <- c(Ha="Northern_Bottlenose_Whale", Mb="Sowerbys_Beaked_Whale",
                    Zc="Goose_Beaked_Whale", MmMe="Trues_Gervais_Beaked_Whale")
  if (RUN_BEAKED_PAM) {
    cat("=== INDIVIDUAL BEAKED WHALE SPECIES PAM KDEs ===\n\n")
    for (sp in unique(beaked_sf$species)) {
      sp_data <- beaked_sf %>% filter(species == sp)
      if (nrow(sp_data) == 0) next
      cat("Processing", beaked_names[sp], "(n =", nrow(sp_data), ")\n")
      do.call(performKDE, c(
        list(data_sf = sp_data, species_col = "species", species_list = sp,
             weight_col = "proportion_det",
             output_prefix = paste0("pam_beaked_", beaked_names[sp]),
             sigma_val = DEFAULT_PAM_SIGMA, use_bw_diggle = FALSE,
             data_source_label = "Beaked Whale PAM"),
        base_pam_args))
    }
  }
  
  #.............................
  # GROUPED BEAKED WHALE PAM----
  #.............................
  
  if (RUN_GROUPED_BEAKED_PAM) {
    cat("\n=== GROUPED BEAKED WHALE PAM KDE ===\n\n")
    beaked_grp <- beaked_sf %>%
      group_by(deployment, station, geometry) %>%
      summarise(detection_days = sum(detection_days, na.rm = TRUE),
                effort_days    = dplyr::first(effort_days), .groups = "drop") %>%
      mutate(proportion_det = detection_days / effort_days,
             beaked_group   = "All Beaked Whales (PAM)")
    do.call(performKDE, c(
      list(data_sf = beaked_grp, species_col = "beaked_group",
           species_list = "All Beaked Whales (PAM)",
           weight_col = "proportion_det",
           output_prefix = "pam_grouped_beaked",
           sigma_val = DEFAULT_PAM_SIGMA, use_bw_diggle = FALSE,
           data_source_label = "Beaked Whale PAM"),
      base_pam_args))
  }
}

#.............................
# SUMMARY------
#.............................

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Resolution:", RESOLUTION, "m | Buffer:", buffer_percent * 100,
    "% | Mode:", ifelse(FAST_MODE, "FAST", "FULL"),
    "| Caching:", ifelse(ENABLE_CACHING, "ON", "OFF"), "\n")
if (RUN_SIGHTINGS)                                         cat("+ Sightings processed\n")
if (RUN_SIGHTINGS && RUN_GROUPED_BEAKED_SIGHTINGS)         cat("+ Grouped beaked sightings KDE\n")
if (RUN_SIGHTINGS && RUN_GROUPED_DEEP_DIVERS_SIGHTINGS)    cat("+ Grouped Deep Divers sightings KDE\n")
if (RUN_PAM)                                               cat("+ Baleen PAM processed\n")
if (RUN_ANY_BEAKED_PAM)                                    cat("+ Beaked whale PAM processed\n")
cat("Rasters:", raster_dir, "| Plots:", plot_dir, "| Shapefiles:", shapefile_dir, "\n")
cat("==============================================================================\n")

# ---- Export parameters for report ----
kde_params <- tibble::tribble(
  ~Parameter,                           ~Value,                                            ~Script,               ~Notes,
  "RESOLUTION",                         as.character(RESOLUTION),                          "08_run_kde_models.R", "KDE raster cell size (m)",
  "DEFAULT_SIGHTINGS_SIGMA",            as.character(DEFAULT_SIGHTINGS_SIGMA),             "08_run_kde_models.R", "Fixed bandwidth sightings KDEs (m)",
  "DEFAULT_PAM_SIGMA",                  as.character(DEFAULT_PAM_SIGMA),                   "08_run_kde_models.R", "Fixed bandwidth PAM KDEs (m)",
  "CONTOUR_QUANTILE",                   as.character(CONTOUR_QUANTILE),                    "08_run_kde_models.R", "Density quantile threshold for 90% contour",
  "buffer_percent",                     as.character(buffer_percent),                      "08_run_kde_models.R", "Buffer added to KDE computation window",
  "USE_PROXIMITY_WEIGHTING_INDIVIDUAL", as.character(USE_PROXIMITY_WEIGHTING_INDIVIDUAL),  "08_run_kde_models.R", "Proximity weights for individual sightings KDEs",
  "USE_PROXIMITY_WEIGHTING_FAMILY",     as.character(USE_PROXIMITY_WEIGHTING_FAMILY),      "08_run_kde_models.R", "Proximity weights for family-level KDEs",
  "PROXIMITY_K_NEIGHBORS",              as.character(PROXIMITY_K_NEIGHBORS),               "08_run_kde_models.R", "k nearest neighbours for proximity weight",
  "PROXIMITY_MAX_DISTANCE",             as.character(PROXIMITY_MAX_DISTANCE),              "08_run_kde_models.R", "Distance at which weight -> min (m)",
  "PROXIMITY_MIN_WEIGHT",               as.character(PROXIMITY_MIN_WEIGHT),                "08_run_kde_models.R", "Floor weight for isolated points",
  "USE_BW_DIGGLE_INDIVIDUAL",           as.character(USE_BW_DIGGLE_INDIVIDUAL),            "08_run_kde_models.R", "Fixed bandwidth used for individual species",
  "USE_BW_DIGGLE_FAMILY",               as.character(USE_BW_DIGGLE_FAMILY),                "08_run_kde_models.R", "Fixed bandwidth used for family groups",
  "DATE_FILTER_SIGHTINGS",              as.character(DATE_FILTER_SIGHTINGS),               "08_run_kde_models.R", "Records before this date excluded from KDE",
  "MIN_RECORDS_SPECIES",                as.character(MIN_RECORDS_SPECIES),                 "08_run_kde_models.R", "Minimum records per species for KDE",
  "species_list_pam",                   paste(species_list_pam, collapse = " "),           "08_run_kde_models.R", "Baleen PAM species codes",
  "RUN_GROUPED_BEAKED_SIGHTINGS",       as.character(RUN_GROUPED_BEAKED_SIGHTINGS),        "08_run_kde_models.R", "Whether grouped beaked whale sightings KDE was run",
  "beaked_species_sightings",           paste(beaked_species_sightings, collapse = "; "),  "08_run_kde_models.R", "Common names pooled for grouped beaked whale sightings KDE",
  "RUN_GROUPED_DEEP_DIVERS_SIGHTINGS",  as.character(RUN_GROUPED_DEEP_DIVERS_SIGHTINGS),   "08_run_kde_models.R", "Whether grouped Deep Divers sightings KDE was run",
  "PHYSETER_REGEX",                     PHYSETER_REGEX,                                    "08_run_kde_models.R", "Regex to identify Physeter sperm whales in cetacean_sf common_name",
  "PAM_BUFFER_DIST_M",                  as.character(PAM_BUFFER_DIST_M),                   "08_run_kde_models.R", "Buffer around study area for PAM station inclusion (m)"
)
readr::write_csv(kde_params, here::here("output/data/params_kde.csv"))
