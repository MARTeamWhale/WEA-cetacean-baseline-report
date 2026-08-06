# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-05-27
# Script: RUN_PIPELINE.R
# Description: Master controller for the WEA cetacean baseline report pipeline.
#              Set the stage toggles below to TRUE/FALSE to control which parts
#              of the pipeline run. Outputs from each stage are saved to disk,
#              so you can skip upstream stages if their outputs already exist.
#
# Changes from previous version:
#   - Initial version: controller script created from numbered pipeline scripts.
#   - Added RUN_KDE_MAPS_ONLY toggle: re-plots KDE PNGs from existing .tif rasters
#     without recomputing KDEs (sets FAST_MODE = TRUE in 08_run_kde_models.R).
#   - Added conflict guard: stops if both RUN_KDE_MODELS and RUN_KDE_MAPS_ONLY are TRUE.
#   - Removed dangling 11_compare_kde_legacy.R block (script deleted from repo).
#   - Updated "Tweak map layouts" preset to include RUN_KDE_MAPS_ONLY.
# ==============================================================================

# ==============================================================================
# SCRIPT TOGGLES ----------------------------------------------------------------
# Set to TRUE to run a stage, FALSE to skip it.
# Read the "When to re-run" notes before changing a toggle.
# ==============================================================================

# ---- GROUP A: One-time setup (run once, or when shapefiles change) ------------
# When to re-run: Only when the study area boundary or shapefiles are updated.
# Requires:  Raw shapefiles in the paths defined in helper_spatial.R
# Produces:  Nothing saved to disk -- builds spatial objects for downstream use.
RUN_SETUP_STUDY_AREA <- F

# ---- GROUP A: Aerial data cleaning (run once per new data delivery) ----------
# When to re-run: Only when a new aerial survey CSV is delivered.
# Requires:  input/raw_data/aerial/<survey_file>.csv
# Produces:  input/processed_data/aerial/deduplicated_sightings.csv
RUN_CLEAN_AERIAL <- F

# ---- GROUP B: Data preparation (run when inputs change) --------------------
# When to re-run: After GROUP A (new aerial data) or when WSDB is updated.
# Requires:  input/processed_data/aerial/deduplicated_sightings.csv
#            input/raw_data/wsdb/WSDB_DFO_Sep2025.csv
# Produces:  output/data/combined_dedup_1km_day.csv  (feeds ALL downstream)
#            output/data/deduplication_*.csv
#            output/data/sightings_*.csv
RUN_COMBINE_SIGHTINGS <- F

# ---- GROUP C: Summaries and plots and sightings plots (fast, re-run freely) -------
# When to re-run: Whenever you tweak plots, palettes, or table formatting.
# Requires:  output/data/combined_dedup_1km_day.csv  (GROUP B output)
# Produces:  output/figs/sights/*.png
#            output/data/params_effort.csv / params_kde.csv / etc.
RUN_PLOT_SIGHTINGS         <- FALSE  # 03_plot_combined_sightings.R
RUN_SUMMARIZE_EFFORT       <- F  # 04_summarize_wsdb_effort.R
RUN_SUMMARIZE_CONFIDENCE   <- FALSE  # 05_summarize_wsdb_confidence.R
RUN_SUMMARIZE_SEASONALITY  <- F  # 06_summarize_wsdb_seasonality.R
RUN_SUMMARIZE_PAM          <- FALSE  # 07_summarize_pam_detections.R

# ---- GROUP D: KDE models (SLOW -- 20+ min, skip if rasters exist) -----------
# When to re-run: Only when sightings data changes (GROUP B) or KDE parameters
#                 change (bandwidth, resolution, species list).
# Requires:  output/data/combined_dedup_1km_day.csv
#            input/raw_data/PAM/baleen_presence_laura_2025.csv
# Produces:  output/tif/*.tif   (KDE rasters -- fed into GROUP E)
#            output/figs/KDE_maps/*.png
RUN_KDE_MODELS <- FALSE

# ---- GROUP D (fast): Re-plot KDE maps only (skips computation) --------------
# When to re-run: After tweaking map aesthetics without
#                 changing the KDE parameters or sightings data.
# Requires:  output/tif/*.tif  (existing rasters from a prior GROUP D run)
# Produces:  output/figs/KDE_maps/*.png   (overwrites previous PNGs)
# NOTE: Cannot be TRUE at the same time as RUN_KDE_MODELS.
RUN_KDE_MAPS_ONLY <- F

# ---- GROUP E: Post-KDE outputs (fast, reads existing rasters) ---------------
# When to re-run: After GROUP D, or to tweak contour quantiles / map aesthetics.
# Requires:  output/tif/*.tif  (GROUP D output)
# Produces:  output/shapes/multi/*.shp
#            output/figs/compare_sightings_pam/*.png
RUN_KDE_CONTOURS          <- FALSE  # 10_create_kde_contours.R  
RUN_COMPARE_PAM_SIGHTINGS <- T  # 11_compare_kde_pam_vs_sightings.R  


# ==============================================================================
# QUICK PRESETS ----------------------------------------------------------------
# Uncomment ONE preset block below to configure a common workflow in one step.
# Comment it back out before customising individual toggles above.
# ==============================================================================

# -- PRESET: Full pipeline (run everything from scratch) -----------------------
# RUN_SETUP_STUDY_AREA <- TRUE; RUN_CLEAN_AERIAL <- TRUE
# RUN_COMBINE_SIGHTINGS <- TRUE
# RUN_PLOT_SIGHTINGS <- TRUE; RUN_SUMMARIZE_EFFORT <- TRUE
# RUN_SUMMARIZE_CONFIDENCE <- TRUE; RUN_SUMMARIZE_SEASONALITY <- TRUE
# RUN_SUMMARIZE_PAM <- TRUE
# RUN_KDE_MODELS <- TRUE
# RUN_KDE_CONTOURS <- TRUE; RUN_COMPARE_PAM_SIGHTINGS <- TRUE

# # -- PRESET: Tweak map layouts only (KDEs already done) -----------------------
# RUN_PLOT_SIGHTINGS <- TRUE
# RUN_KDE_MAPS_ONLY  <- TRUE
# RUN_COMPARE_PAM_SIGHTINGS <- TRUE

# -- PRESET: New data delivery (re-run from aerial clean through KDEs) ---------
# RUN_CLEAN_AERIAL <- TRUE; RUN_COMBINE_SIGHTINGS <- TRUE
# RUN_PLOT_SIGHTINGS <- TRUE; RUN_SUMMARIZE_EFFORT <- TRUE
# RUN_SUMMARIZE_CONFIDENCE <- TRUE; RUN_SUMMARIZE_SEASONALITY <- TRUE
# RUN_SUMMARIZE_PAM <- TRUE
# RUN_KDE_MODELS <- TRUE
# RUN_KDE_CONTOURS <- TRUE; RUN_COMPARE_PAM_SIGHTINGS <- TRUE

# -- PRESET: Re-run summaries and PAM only (no KDEs) --------------------------
# RUN_SUMMARIZE_EFFORT <- TRUE; RUN_SUMMARIZE_CONFIDENCE <- TRUE
# RUN_SUMMARIZE_SEASONALITY <- TRUE; RUN_SUMMARIZE_PAM <- TRUE


# ==============================================================================
# EXECUTE PIPELINE -------------------------------------------------------------
# Do not edit below this line unless you are modifying the pipeline structure.
# ==============================================================================

library(here)

prefer_tidyverse_conflicts <- function() {
  if (!requireNamespace("conflicted", quietly = TRUE)) {
    return(invisible(NULL))
  }
  
  conflicted::conflicts_prefer(
    dplyr::arrange,
    dplyr::between,
    dplyr::filter,
    dplyr::first,
    dplyr::lag,
    dplyr::last,
    dplyr::mutate,
    dplyr::rename,
    dplyr::select,
    dplyr::summarise,
    lubridate::month,
    lubridate::year,
    purrr::map,
    purrr::walk,
    .quiet = TRUE
  )
}

prefer_tidyverse_conflicts()

cat("\n")
cat("==============================================================================\n")
cat("WEA CETACEAN BASELINE REPORT -- PIPELINE CONTROLLER\n")
cat("==============================================================================\n\n")


# ---- GROUP A: One-time setup (~1-3 min each, run once per data delivery) -----

## 01_clean_aerial_sightings.R ----
if (RUN_CLEAN_AERIAL) {
  cat(">>> Running setup/01_clean_aerial_sightings.R ...\n")
  source(here("scripts", "setup", "01_clean_aerial_sightings.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped setup/01_clean_aerial_sightings.R (RUN_CLEAN_AERIAL = FALSE)\n\n")
  aerial_out <- here("input", "processed_data", "aerial", "deduplicated_sightings.csv")
  if (!file.exists(aerial_out)) {
    stop("01_clean_aerial_sightings.R was skipped but its output is missing:\n  ",
         aerial_out, "\nSet RUN_CLEAN_AERIAL <- TRUE and re-run, or place the file manually.")
  }
}

## 02_build_study_area.R ----
if (RUN_SETUP_STUDY_AREA) {
  cat(">>> Running setup/02_build_study_area.R ...\n")
  source(here("scripts", "setup", "02_build_study_area.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped setup/02_build_study_area.R (RUN_SETUP_STUDY_AREA = FALSE)\n\n")
}


# ---- GROUP B: Data preparation (~1-2 min, re-run when inputs change) ---------

## 02_prepare_combined_sightings.R ----
if (RUN_COMBINE_SIGHTINGS) {
  cat(">>> Running 02_prepare_combined_sightings.R ...\n")
  source(here("scripts", "02_prepare_combined_sightings.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 02_prepare_combined_sightings.R (RUN_COMBINE_SIGHTINGS = FALSE)\n\n")
  combined_out <- here("output", "data", "combined_dedup_1km_day.csv")
  if (!file.exists(combined_out)) {
    stop("02_prepare_combined_sightings.R was skipped but its output is missing:\n  ",
         combined_out, "\nSet RUN_COMBINE_SIGHTINGS <- TRUE and re-run.")
  }
}


# ---- GROUP C: Summaries and plots (~2-10 min each, re-run freely) ------------

## 03_plot_combined_sightings.R ----
if (RUN_PLOT_SIGHTINGS) {
  cat(">>> Running 03_plot_combined_sightings.R ...\n")
  source(here("scripts", "03_plot_combined_sightings.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 03_plot_combined_sightings.R (RUN_PLOT_SIGHTINGS = FALSE)\n\n")
}

## 04_summarize_wsdb_effort.R ----
if (RUN_SUMMARIZE_EFFORT) {
  cat(">>> Running 04_summarize_wsdb_effort.R ...\n")
  source(here("scripts", "04_summarize_wsdb_effort.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 04_summarize_wsdb_effort.R (RUN_SUMMARIZE_EFFORT = FALSE)\n\n")
}

## 05_summarize_wsdb_confidence.R  (needs 04 objects in memory) ----
if (RUN_SUMMARIZE_CONFIDENCE) {
  if (!RUN_SUMMARIZE_EFFORT && !exists("cell_all")) {
    stop("05_summarize_wsdb_confidence.R requires objects from 04_summarize_wsdb_effort.R.\n",
         "Either set RUN_SUMMARIZE_EFFORT <- TRUE, or run 04 first in the same R session.")
  }
  cat(">>> Running 05_summarize_wsdb_confidence.R ...\n")
  source(here("scripts", "05_summarize_wsdb_confidence.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 05_summarize_wsdb_confidence.R (RUN_SUMMARIZE_CONFIDENCE = FALSE)\n\n")
}

## 06_summarize_wsdb_seasonality.R ----
if (RUN_SUMMARIZE_SEASONALITY) {
  if (!RUN_SUMMARIZE_EFFORT && !exists("cell_m")) {
    stop("06_summarize_wsdb_seasonality.R requires objects from 04_summarize_wsdb_effort.R.\n",
         "Either set RUN_SUMMARIZE_EFFORT <- TRUE, or run 04 first in the same R session.")
  }
  cat(">>> Running 06_summarize_wsdb_seasonality.R ...\n")
  source(here("scripts", "06_summarize_wsdb_seasonality.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 06_summarize_wsdb_seasonality.R (RUN_SUMMARIZE_SEASONALITY = FALSE)\n\n")
}

## 07_summarize_pam_detections.R ----
if (RUN_SUMMARIZE_PAM) {
  cat(">>> Running 07_summarize_pam_detections.R ...\n")
  source(here("scripts", "07_summarize_pam_detections.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 07_summarize_pam_detections.R (RUN_SUMMARIZE_PAM = FALSE)\n\n")
}


# ---- GROUP D: KDE models (SLOW -- 20-40 min, skip if rasters exist) ----------

## 08_run_kde_models.R ----
if (RUN_KDE_MODELS) {
  cat(">>> Running 08_run_kde_models.R (20-40 min) ...\n")
  source(here("scripts", "08_run_kde_models.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 08_run_kde_models.R (RUN_KDE_MODELS = FALSE)\n\n")
  tif_files <- list.files(here("output", "tif"), pattern = "\\.tif$")
  if (length(tif_files) == 0 && any(RUN_KDE_CONTOURS, RUN_COMPARE_PAM_SIGHTINGS)) {
    warning("08_run_kde_models.R was skipped but no .tif rasters found in output/tif/.\n",
            "Scripts 09-10 need these rasters. Set RUN_KDE_MODELS <- TRUE first.")
  }
}

## 09_plot_kde_maps.R ----
if (RUN_KDE_MAPS_ONLY) {
  cat(">>> Running 09_plot_kde_maps.R ...\n")
  source(here("scripts", "09_plot_kde_maps.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 09_plot_kde_maps.R (RUN_KDE_MAPS_ONLY = FALSE)\n\n")
}


# ---- GROUP E: Post-KDE outputs (~2-5 min each, reads existing rasters) -------

## 10_create_kde_contours.R ----
if (RUN_KDE_CONTOURS) {
  cat(">>> Running 10_create_kde_contours.R ...\n")
  source(here("scripts", "10_create_kde_contours.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 10_create_kde_contours.R (RUN_KDE_CONTOURS = FALSE)\n\n")
}

## 11_compare_kde_pam_vs_sightings.R ----
if (RUN_COMPARE_PAM_SIGHTINGS) {
  cat(">>> Running 11_compare_kde_pam_vs_sightings.R ...\n")
  source(here("scripts", "11_compare_kde_pam_vs_sightings.R"))
  cat("    Done.\n\n")
} else {
  cat("--- Skipped 11_compare_kde_pam_vs_sightings.R (RUN_COMPARE_PAM_SIGHTINGS = FALSE)\n\n")
}


# ---- Pipeline complete -------------------------------------------------------
cat("==============================================================================\n")
cat("PIPELINE COMPLETE\n")
cat("==============================================================================\n\n")
