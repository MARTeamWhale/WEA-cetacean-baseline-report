# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-05-27
# Script: 09_plot_kde_maps.R
# Description: Standalone script to re-generate KDE map PNGs from existing
#              .tif rasters. Does NOT recompute KDEs and does NOT regenerate
#              contour shapefiles. Run this after tweaking map aesthetics
#              (colours, themes, point styling) without touching the models.
#              Requires: output/tif/*.tif from a previous 08_run_kde_models.R run.
#
# Changes from previous version:
#   - Initial version: extracted from 08_run_kde_models.R plotting section.
# ==============================================================================


# ==============================================================================
# USER SETTINGS ----------------------------------------------------------------
# ==============================================================================

# ---- Toggle observation point overlay ----------------------------------------
# TRUE  = overlay sightings / PAM station points on each map.
# FALSE = plot raster and contour only (faster, no data loading required).
PLOT_OBSERVATION_POINTS <- TRUE

# ---- Contour quantile (must match the value used in 08_run_kde_models.R) -----
CONTOUR_QUANTILE <- 0.90

# ---- Output directory --------------------------------------------------------
# Overwrites PNGs in the same folder as 08_run_kde_models.R produces them.
PLOT_OUTPUT_DIR  <- "output/figs/KDE_maps/"


# ==============================================================================
# SETUP ------------------------------------------------------------------------
# ==============================================================================

# ---- Packages ----------------------------------------------------------------
library(here)
source(here("scripts", "00_load_helpers.R"))   # loads helper functions + spatial paths

pacman::p_load(terra, sf, ggplot2, dplyr, stringr, readr, patchwork)

# ---- CRS ---------------------------------------------------------------------
UTM20 <- SPATIAL_CRS_UTM20

# ---- Output folder -----------------------------------------------------------
if (!dir.exists(PLOT_OUTPUT_DIR)) dir.create(PLOT_OUTPUT_DIR, recursive = TRUE)


# ==============================================================================
# LOAD SPATIAL BASE LAYERS -----------------------------------------------------
# ==============================================================================

cat("\n=== Loading spatial base layers ===\n")

study_area <- load_spatial_layer("study_area", crs = 32620)
coastline  <- load_spatial_layer("land",       crs = 32620, required = FALSE)
owa        <- load_spatial_layer("wea",        crs = 32620, required = FALSE)

if (!is.null(coastline)) {
  buf <- st_buffer(study_area, dist = 10000)
  coastline <- suppressWarnings(st_crop(coastline, st_bbox(buf)))
}
if (!is.null(owa)) {
  buf <- st_buffer(study_area, dist = 10000)
  owa <- suppressWarnings(st_crop(owa, st_bbox(buf)))
}

# Plot extent: bounding box of study area with a small buffer
bb    <- st_bbox(st_buffer(study_area, dist = 5000))
xlims <- c(bb["xmin"], bb["xmax"])
ylims <- c(bb["ymin"], bb["ymax"])

cat("Base layers loaded.\n\n")


# ==============================================================================
# OPTIONALLY LOAD OBSERVATION POINTS -------------------------------------------
# ==============================================================================

cetacean_sf <- NULL
pam_sf      <- NULL

if (PLOT_OBSERVATION_POINTS) {
  cat("=== Loading observation points ===\n")

  # ---- Sightings --------------------------------------------------------------
  sightings_csv <- here("output", "data", "combined_dedup_1km_day.csv")
  if (file.exists(sightings_csv)) {
    cetacean_sf <- read_csv(sightings_csv, show_col_types = FALSE) %>%
      filter(family %in% c("Odontocete", "Mysticete")) %>%
      st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
      st_transform(UTM20) %>%
      st_make_valid()
    cat("Sightings loaded:", nrow(cetacean_sf), "records\n")
  } else {
    cat("Sightings CSV not found -- skipping point overlay for sightings maps.\n")
  }

  # ---- PAM --------------------------------------------------------------------
  pam_raw <- here("input", "raw_data", "PAM", "baleen_presence_laura_2025.csv")
  if (file.exists(pam_raw)) {
    pam_sf <- read.csv(pam_raw) %>%
      mutate(rec_date = as.Date(rec_date)) %>%
      group_by(site, latitude, longitude, species) %>%
      summarise(detection_days = sum(presence, na.rm = TRUE),
                effort_days    = n_distinct(rec_date),
                proportion_det = detection_days / effort_days,
                .groups        = "drop") %>%
      filter(effort_days > 0) %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_transform(UTM20)
    cat("PAM stations loaded:", nrow(pam_sf), "records\n")
  } else {
    cat("PAM CSV not found -- skipping point overlay for PAM maps.\n")
  }

  cat("\n")
}


# ==============================================================================
# FIND RASTERS -----------------------------------------------------------------
# ==============================================================================

tif_paths <- list.files(here("output", "tif"), pattern = "\\.tif$", full.names = TRUE)
if (length(tif_paths) == 0) {
  stop("No .tif rasters found in output/tif/.\n",
       "Run 08_run_kde_models.R first (set RUN_KDE_MODELS <- TRUE in RUN_PIPELINE.R).")
}
cat("=== Found", length(tif_paths), "rasters -- generating plots ===\n\n")


# ==============================================================================
# HELPER: parse a tif filename into a human-readable title ---------------------
# ==============================================================================

parse_tif_title <- function(path) {
  # filename format: KDE_{output_prefix}_{species_clean}_{bw_label}.tif
  base <- tools::file_path_sans_ext(basename(path))
  # strip leading "KDE_"
  rest <- sub("^KDE_", "", base)
  # strip trailing bandwidth label (e.g. _fixed_bw10000 or _diggle_bw*)
  rest <- sub("_(fixed|diggle)_bw[0-9*]+$", "", rest)
  # replace underscores with spaces for readability
  gsub("_", " ", rest)
}


# ==============================================================================
# HELPER: load matching q90 contour shapefile if it exists --------------------
# ==============================================================================

load_contour <- function(tif_path) {
  # Shapefile name mirrors the tif name: strip KDE_ prefix, append _contour90.shp
  # e.g. KDE_pam_baleen_Blue_Whale_Bm_fixed_bw10000.tif
  #   -> pam_baleen_Blue_Whale_Bm_fixed_bw10000_contour90.shp
  base     <- tools::file_path_sans_ext(basename(tif_path))
  shp_stem <- sub("^KDE_", "", base)
  shp_path <- here("output", "shapes", "q90", paste0(shp_stem, "_contour90.shp"))
  if (file.exists(shp_path)) {
    tryCatch(st_read(shp_path, quiet = TRUE), error = function(e) NULL)
  } else NULL
}


# ==============================================================================
# HELPER: get observation points for a given raster ---------------------------
# ==============================================================================

get_obs_points <- function(tif_path) {
  if (!PLOT_OBSERVATION_POINTS) return(NULL)
  base   <- tools::file_path_sans_ext(basename(tif_path))
  is_pam <- grepl("pam", base, ignore.case = TRUE)
  if (is_pam) pam_sf else cetacean_sf
}


# ==============================================================================
# PLOT EACH RASTER -------------------------------------------------------------
# ==============================================================================

for (tif_path in tif_paths) {

  base  <- tools::file_path_sans_ext(basename(tif_path))
  title <- parse_tif_title(tif_path)
  is_pam <- grepl("pam", base, ignore.case = TRUE)

  cat("  Plotting:", title, "\n")

  # ---- Load raster ------------------------------------------------------------
  rast_obj <- tryCatch(rast(tif_path), error = function(e) {
    cat("    ERROR loading raster:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(rast_obj)) next

  kde_df <- as.data.frame(rast_obj, xy = TRUE) %>%
    filter(!is.na(lyr.1)) %>%
    mutate(lyr.1_sqrt = sqrt(lyr.1))

  if (nrow(kde_df) == 0) { cat("    Raster is empty -- skipping.\n"); next }

  # ---- Load contour shapefile (no regeneration) -------------------------------
  contour_sf <- load_contour(tif_path)

  # ---- Observation points -----------------------------------------------------
  obs_pts <- get_obs_points(tif_path)

  # ---- Build plot -------------------------------------------------------------
  p <- ggplot() +
    geom_tile(data = kde_df, aes(x = x, y = y, fill = lyr.1_sqrt)) +
    scale_fill_viridis_c(option = "viridis", name = "sqrt(Density)",
                         limits = c(0, max(kde_df$lyr.1_sqrt)),
                         begin = 0.15, end = 1.0)

  if (!is.null(contour_sf) && nrow(contour_sf) > 0)
    p <- p + map_layer_kde_contour(contour_sf)

  if (!is.null(coastline))
    p <- p + geom_sf(data = coastline, fill = "grey90", color = "grey70", linewidth = 0.3)

  if (!is.null(owa))
    p <- p + map_layer_wea(owa)

  if (!is.null(obs_pts) && nrow(obs_pts) > 0)
    p <- p + geom_sf(data = obs_pts, color = "white",
                     size  = ifelse(is_pam, 1.5, 0.8),
                     shape = 21, fill = "black",
                     stroke = ifelse(is_pam, 0.5, 0.3),
                     alpha  = ifelse(is_pam, 0.8, 0.6))

  if (!is.null(study_area))
    p <- p + map_layer_study_area(study_area)

  p <- p +
    coord_sf(crs = st_crs(UTM20), xlim = xlims, ylim = ylims, expand = FALSE) +
    map_base_theme() +
    labs(title = title) +
    theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
          panel.grid    = element_line(color = "grey95"))

  ct_text <- paste0(CONTOUR_QUANTILE * 100, "% contour (baseline)")
  p <- p + map_annotation_label(xlims[1], ylims[2], ct_text, vjust = 1.3)

  # ---- Save -------------------------------------------------------------------
  out_path <- file.path(PLOT_OUTPUT_DIR, paste0(base, ".png"))
  ggsave(out_path, p, width = 8, height = 6, dpi = 300)
  cat("    Saved:", basename(out_path), "\n")
}

cat("\n=== Done -- PNGs written to", PLOT_OUTPUT_DIR, "===\n\n")
