# ==============================================================================
# Author: Laura Joan Feyrer
# Date Updated: 2026-03-15
# Script: 10_compare_kde_pam_vs_sightings.R
# Description: Creates side-by-side  SIGHTINGS vs PAM KDE comparison maps for
#              cetacean species. Outputs PNG maps per species to
#              output/figs/compare_sightings_pam/. Also exports key params
#              to output/data/params_compare.csv for the project report.
#
# Changes from previous version (2026-03-13):
#   - GROUPED_BEAKED entry updated: sightings_base now points to the grouped
#     beaked sightings KDE from v6 KDE script (sightings_beaked_grouped),
#     and pam_only changed from TRUE to FALSE so a  SIGHTINGS vs PAM comparison
#     map is produced instead of a PAM-only map.
#   - Added BEAKED_SPECIES_SIGHTINGS vector at top of script (mirrors
#     beaked_species_sightings in v6 KDE script) defining which common_name
#     values are pooled for the grouped beaked sightings panel.
#   - Added "BeakedWhale" group_type filter branch in create_comparison_map()
#     that pools sightings_data records matching BEAKED_SPECIES_SIGHTINGS and
#     passes all beaked_pam_data records to the PAM panel.
#   - Grouped beaked run block updated: pam_only = FALSE, passes beaked_pam_data.
#   - GROUPED_OdontoceteS_SIGHTINGS renamed to GROUPED_ODONTOCETES_SIGHTINGS
#     for consistency.
#   - params export updated to include BEAKED_SPECIES_SIGHTINGS.
# ==============================================================================

options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, ggplot2, viridis, ggspatial, patchwork, rlang)
suppressWarnings(source(here::here("scripts/00_load_helpers.R")))


# ==============================================================================
# CONFIGURATION ----
# ==============================================================================

CREATE_SIGHTINGS_ONLY_MAPS <- F

# ---- Projections ----
UTM20 <- st_crs(32620)

# ---- Input paths ----
SHAPEFILE_DIR    <- "output/shapes/multi"
KDE_SHAPEFILE_DIRS <- c("output/shapes/multi", "output/shapes/q90")
SIGHTINGS_DATA   <- "output/data/combined_dedup_1km_day.csv"
BALEEN_PAM_DATA  <- "input/raw_data/PAM/baleen_presence_laura_2025.csv"
BEAKED_PAM_DATA  <- "input/raw_data/PAM/beaked_pam_results_2026-01-12.csv"

# ---- Output directory ----
COMPARE_OUTPUT_DIR <- "output/figs/compare_sightings_pam/"
if (!dir.exists(COMPARE_OUTPUT_DIR)) dir.create(COMPARE_OUTPUT_DIR, recursive = TRUE)
old_compare_pngs <- list.files(COMPARE_OUTPUT_DIR, pattern = "\\.png$",
                               full.names = TRUE, ignore.case = TRUE)
if (length(old_compare_pngs) > 0) {
  unlink(old_compare_pngs)
  cat("Cleared", length(old_compare_pngs), "old comparison PNGs from",
      COMPARE_OUTPUT_DIR, "\n")
}

# ---- Colour palette settings ----
QUANTILE_PALETTE  <- "mako"
QUANTILE_HIDE     <- c("85%")   # contour levels to drop before plotting (too light to see)
POINT_COLOR      <- "gray40"
STATION_COLOR    <- "orange"
OSW_COLOR        <- "red"
OSW_ALPHA        <- 0.3

# ---- PAM inclusion buffer ----
# Stations within this distance of the study area boundary are included.
# Must match run_all_kdes and 07_summarize_pam_detections.R.
PAM_BUFFER_DIST_M <- 10000  # 10 km

# ---- Beaked whale sightings species list ----
# Common names pooled for the grouped beaked whale sightings panel.
# Must match beaked_species_sightings in the v6 KDE script exactly.
BEAKED_SPECIES_SIGHTINGS <- c(
  "Northern Bottlenose Whale",
  "Sowerby's Beaked Whale",
  "True's Beaked Whale",
  "Gervais' Beaked Whale",
  "Blainville's Beaked Whale",
  "Cuvier's Beaked Whale"
)


# ==============================================================================
# SPECIES CONFIG: BALEEN WHALES (PAM + SIGHTINGS) ----
# ==============================================================================

BALEEN_SPECIES_MAPPING <- list(
  "Sei_Whale" = list(
    pam_code       = "Bb",
    common_name    = "Sei Whale",
    display_name   = "Sei Whale",
    sightings_base = "Sei_Whale",
    pam_base       = "pam_baleen_Sei_Whale_Bb_fixed_bw10000"
  ),
  "Blue_Whale" = list(
    pam_code       = "Bm",
    common_name    = "Blue Whale",
    display_name   = "Blue Whale",
    sightings_base = "Blue_Whale",
    pam_base       = "pam_baleen_Blue_Whale_Bm_fixed_bw10000"
  ),
  "Fin_Whale" = list(
    pam_code       = "Bp",
    common_name    = "Fin Whale",
    display_name   = "Fin Whale",
    sightings_base = "Fin_Whale",
    pam_base       = "pam_baleen_Fin_Whale_Bp_fixed_bw10000"
  ),
  "Humpback_Whale" = list(
    pam_code       = "Mn",
    common_name    = "Humpback Whale",
    display_name   = "Humpback Whale",
    sightings_base = "Humpback_Whale",
    pam_base       = "pam_baleen_Humpback_Whale_Mn_fixed_bw10000"
  ),
  "Minke_Whale" = list(
    pam_code       = "Ba",
    common_name    = "Minke Whale",
    display_name   = "Minke Whale",
    sightings_base = "Minke_Whale",
    pam_base       = "pam_baleen_Minke_Whale_Ba_fixed_bw10000"
  ),
  "North_Atlantic_Right_Whale" = list(
    pam_code       = "Eg",
    common_name    = "North Atlantic Right Whale",
    display_name   = "North Atlantic Right Whale",
    sightings_base = "North_Atlantic_Right_Whale",
    pam_base       = "pam_baleen_North_Atlantic_Right_Whale_Eg_fixed_bw10000"
  )
)


# ==============================================================================
# SPECIES CONFIG: BEAKED WHALES ----
# ==============================================================================

BEAKED_SPECIES_MAPPING <- list(
  "Northern_Bottlenose_Whale" = list(
    pam_code       = "Ha",
    common_name    = "Northern Bottlenose Whale",
    display_name   = "Northern Bottlenose Whale",
    sightings_base = "Northern_Bottlenose_Whale",
    pam_base       = "pam_beaked_Northern_Bottlenose_Whale_Ha_fixed_bw10000",
    has_sightings  = TRUE
  ),
  "Sowerbys_Beaked_Whale" = list(
    pam_code       = "Mb",
    common_name    = "Sowerby's Beaked Whale",
    display_name   = "Sowerby's Beaked Whale",
    sightings_base = "Sowerbys_Beaked_Whale",
    pam_base       = "pam_beaked_Sowerbys_Beaked_Whale_Mb_fixed_bw10000",
    has_sightings  = TRUE
  ),
  "Goose_Beaked_Whale" = list(
    pam_code       = "Zc",
    common_name    = "Goose Beaked Whale",
    display_name   = "Goose Beaked Whale",
    sightings_base = "Goose_Beaked_Whale",
    pam_base       = "pam_beaked_Goose_Beaked_Whale_Zc_fixed_bw10000",
    has_sightings  = FALSE,
    alternate_sightings_names = c("WHALE- CUVIER'S BEAKED", "WHALE-CUVIER'S BEAKED")
  ),
  "Trues_Gervais_Beaked_Whale" = list(
    pam_code       = "MmMe",
    common_name    = "True's/Gervais' Beaked Whale",
    display_name   = "True's/Gervais' Beaked Whale",
    sightings_base = "Trues_Gervais_Beaked_Whale",
    pam_base       = "pam_beaked_Trues_Gervais_Beaked_Whale_MmMe_fixed_bw10000",
    has_sightings  = FALSE
  )
)


# ==============================================================================
# SPECIES CONFIG: GROUPED ----
# ==============================================================================

GROUPED_BALEEN <- list(
  "All_Baleen_Whales" = list(
    common_name    = "ALL_BALEEN_MysticeteS",
    display_name   = "All Baleen Whales",
    sightings_base = "KDE_sightings_mysticetes",
    pam_base       = "pam_grouped_baleen_All_Baleen_Whales_(PAM)_fixed_bw10000",
    is_grouped     = TRUE,
    group_type     = "Mysticete"
  )
)

# Updated from v1: sightings_base now points to the grouped beaked sightings
# KDE produced by the v6 KDE script. pam_only = FALSE produces a PAM vs
# sightings comparison. group_type "BeakedWhale" triggers the dedicated filter
# branch in create_comparison_map() which pools BEAKED_SPECIES_SIGHTINGS.
GROUPED_BEAKED <- list(
  "All_Beaked_Whales" = list(
    common_name    = "ALL_BEAKED_WHALES",
    display_name   = "All Beaked Whales",
    sightings_base = "sightings_beaked_grouped",
    pam_base       = "pam_grouped_beaked_All_Beaked_Whales_(PAM)_fixed_bw10000",
    is_grouped     = TRUE,
    group_type     = "BeakedWhale",
    pam_only       = FALSE
  )
)

GROUPED_ODONTOCETES_SIGHTINGS <- list(
  "All_Odontocetes" = list(
    common_name    = "OdontoceteS",
    display_name   = "All Odontocetes (Toothed Whales)",
    sightings_base = "KDE_sightings_Odontocetes",
    pam_base       = NA_character_,
    is_grouped     = TRUE,
    group_type     = "Odontocete",
    pam_only       = FALSE
  )
)


# ==============================================================================
# LOAD BASE LAYERS ----
# ==============================================================================

cat("Loading base layers...\n")

land <- load_spatial_layer("land", crs = 32620)
cont <- load_bathy_contours(crs = UTM20, required = FALSE)

study_area <- load_spatial_layer("study_area", crs = 32620)
study_area_pam_buffer <- st_buffer(study_area, dist = PAM_BUFFER_DIST_M)

bbox  <- st_bbox(study_area_pam_buffer)
xlims <- c(bbox["xmin"], bbox["xmax"])
ylims <- c(bbox["ymin"], bbox["ymax"])

cat("Loading OSW wind areas...\n")
osw_wind <- load_spatial_layer("wea", crs = 32620)

cat("+ Base layers loaded:", nrow(osw_wind), "OSW wind lease areas\n\n")


# ==============================================================================
# LOAD SOURCE DATA ----
# ==============================================================================

cat("Loading source data...\n")

# ---- Sightings ----
sightings_data <- read_csv(SIGHTINGS_DATA, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  filter(as.Date(date_utc) >= as.Date("2010-01-01")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM20)

cat("+ Loaded", nrow(sightings_data), "sightings (2010+)\n")

# ---- Baleen PAM ----
baleen_pam_raw <- read_csv(BALEEN_PAM_DATA, show_col_types = FALSE) %>%
  mutate(rec_date = as.Date(rec_date), month = month(rec_date)) %>%
  group_by(site, latitude, longitude, species) %>%
  summarise(
    total_days     = n(),
    detection_days = sum(presence, na.rm = TRUE),
    proportion_det = detection_days / total_days,
    .groups        = "drop"
  ) %>%
  filter(total_days > 0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

baleen_pam_data <- st_filter(baleen_pam_raw, study_area_pam_buffer)
cat("+ Loaded", nrow(baleen_pam_data), "baleen PAM records within",
    PAM_BUFFER_DIST_M / 1000, "km of study area\n")

# ---- Beaked whale PAM ----
cat("Loading beaked whale PAM data...\n")
beaked_pam_raw <- read_csv(BEAKED_PAM_DATA, show_col_types = FALSE)%>%rename("site" = "deployment")

beaked_species_codes <- unique(vapply(BEAKED_SPECIES_MAPPING, `[[`, character(1), "pam_code"))

beaked_station_effort <- beaked_pam_raw %>%
  group_by(station) %>%
  summarise(effort_days = n_distinct(rec_date), .groups = "drop")

beaked_pam_observed <- beaked_pam_raw %>%
  group_by(station, latitude, longitude, species) %>%
  summarise(
    detection_days = sum(presence, na.rm = TRUE),
    .groups        = "drop"
  )

beaked_pam_processed <- tidyr::crossing(
  beaked_station_effort
) %>%
  left_join(
    beaked_pam_observed,
    by = c( "station")
  ) %>%group_by(station, species)%>%
  mutate(
    detection_days = sum(detection_days, na.rm = TRUE)
  ) %>% ungroup() %>% 
  mutate(
proportion_det = detection_days / effort_days)%>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

tidyr::crossing(beaked_station_effort, species = beaked_species_codes)

beaked_pam_data <- st_filter(beaked_pam_processed, study_area_pam_buffer)
cat("+ Loaded", nrow(beaked_pam_data), "beaked whale PAM records within",
    PAM_BUFFER_DIST_M / 1000, "km of study area\n")

beaked_species_summary <- beaked_pam_data %>%
  st_drop_geometry() %>%
  group_by(species) %>%
  summarise(
    n_stations   = n(),
    total_det_days  = sum(detection_days),
    mean_proportion = mean(proportion_det),
    .groups = "drop"
  ) %>%
  arrange(desc(total_det_days))

cat("\nBeaked whale species in PAM data:\n")
print(beaked_species_summary, n = Inf)
cat("\n")

# Keep for backward compatibility
pam_data <- baleen_pam_data


# ==============================================================================
# HELPER FUNCTIONS ----
# ==============================================================================

# ---- extract_kde_params: parse bandwidth / method / proximity from filename ----
extract_kde_params <- function(filename) {
  fname <- basename(filename)
  
  has_proximity <- grepl("proximity", fname, ignore.case = TRUE)
  
  bw  <- NA_real_
  m   <- regexec("bw[_-]?(\\d+)", fname, ignore.case = TRUE)
  hit <- regmatches(fname, m)[[1]]
  if (length(hit) > 1) bw <- as.numeric(hit[2])
  
  method <- dplyr::case_when(
    grepl("diggle",   fname, ignore.case = TRUE) ~ "diggle",
    grepl("adaptive", fname, ignore.case = TRUE) ~ "adaptive",
    grepl("fixed",    fname, ignore.case = TRUE) ~ "fixed",
    TRUE ~ "unknown"
  )
  
  list(has_proximity = has_proximity, bandwidth = bw, method = method)
}

# ---- find_kde_shapefiles: locate matching KDE shapefiles for a base name ----
find_kde_shapefiles <- function(base_name,
                                shapefile_dir    = KDE_SHAPEFILE_DIRS,
                                prefer_proximity = FALSE,
                                pam_only         = FALSE) {
  if (is.null(base_name) || is.na(base_name)) return(NULL)
  
  files <- unlist(lapply(
    shapefile_dir,
    list.files,
    pattern = "\\.shp$",
    full.names = TRUE,
    ignore.case = TRUE
  ), use.names = FALSE)
  
  base_esc  <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", base_name)
  token_pat <- paste0("(?<![A-Za-z0-9])", base_esc, "(?![A-Za-z0-9])")
  
  hits <- files[
    grepl("^(KDE_)?", basename(files), ignore.case = TRUE) &
      grepl(token_pat, basename(files), ignore.case = TRUE, perl = TRUE)
  ]
  
  if (length(hits) == 0) return(NULL)
  
  hits <- hits[!grepl("old|backup|test", basename(hits), ignore.case = TRUE)]
  
  # Separate PAM and sightings files explicitly so a sightings search never
  # picks up a PAM shapefile and vice versa.
  is_pam_file <- grepl("^(KDE_)?pam_", basename(hits), ignore.case = TRUE)
  if (pam_only) {
    hits <- hits[is_pam_file]
  } else {
    hits <- hits[!is_pam_file]
  }
  
  if (length(hits) == 0) return(NULL)
  
  hits <- hits[!grepl("old|backup|test", basename(hits), ignore.case = TRUE)]
  
  df <- tibble::tibble(
    path   = hits,
    file   = basename(hits),
    mtime  = file.info(hits)$mtime,
    params = purrr::map(hits, extract_kde_params)
  ) |>
    dplyr::mutate(
      has_proximity = purrr::map_lgl(params, ~ .x$has_proximity),
      bandwidth     = purrr::map_dbl(params, ~ .x$bandwidth %||% NA_real_),
      method        = purrr::map_chr(params, ~ .x$method    %||% "unknown")
    )

  if (any(grepl("_multi\\.shp$", df$file, ignore.case = TRUE))) {
    df <- df |> dplyr::filter(grepl("_multi\\.shp$", file, ignore.case = TRUE))
  }
  
  if (prefer_proximity) {
    df <- df |>
      dplyr::mutate(priority = dplyr::if_else(has_proximity, 0L, 1L)) |>
      dplyr::arrange(priority, dplyr::desc(mtime))
  } else {
    df <- df |> dplyr::arrange(dplyr::desc(mtime))
  }
  
  df
}

# ---- load_kde_shapefile: read + standardise a KDE shapefile ----
load_kde_shapefile <- function(pattern_or_file, shapefile_dir = SHAPEFILE_DIR) {
  
  if (file.exists(pattern_or_file) && grepl("\\.shp$", pattern_or_file)) {
    shp_file <- pattern_or_file
  } else {
    shp_files <- list.files(shapefile_dir,
                            pattern    = paste0(pattern_or_file, "\\.shp$"),
                            full.names = TRUE, ignore.case = TRUE)
    
    if (length(shp_files) == 0) return(NULL)
    
    if (length(shp_files) > 1) {
      shp_file <- shp_files[which.max(file.info(shp_files)$mtime)]
      cat("  Note: multiple matches - using most recent:", basename(shp_file), "\n")
    } else {
      shp_file <- shp_files[1]
    }
  }
  
  kde_shp <- st_read(shp_file, quiet = TRUE) %>% st_transform(UTM20)
  
  # Standardise quantile column name
  quantile_col <- NULL
  for (col_name in c("Quantile", "Quantil", "quant_level", "quantil", "QUANTILE", "contour")) {
    if (col_name %in% names(kde_shp)) { quantile_col <- col_name; break }
  }
  
  if (is.null(quantile_col)) {
    cat("  WARNING: no quantile column found in", basename(shp_file), "\n")
    return(NULL)
  }
  
  kde_shp <- kde_shp %>% rename(Quantile = !!quantile_col)
  
  # Convert numeric quantiles to "XX%" labels
  unique_q <- unique(kde_shp$Quantile)
  unique_q <- unique_q[!is.na(unique_q)]
  if (is.numeric(unique_q)) {
    kde_shp <- kde_shp %>% mutate(Quantile = paste0(Quantile * 100, "%"))
  } else if (all(grepl("^0?\\.\\d+$", as.character(unique_q)))) {
    kde_shp <- kde_shp %>%
      mutate(Quantile = paste0(as.numeric(as.character(Quantile)) * 100, "%"))
  }
  
  # Set ordered factor
  q_levels <- unique(kde_shp$Quantile)
  q_nums   <- as.numeric(gsub("%", "", q_levels))
  kde_shp$Quantile <- factor(kde_shp$Quantile,
                             levels  = q_levels[order(q_nums)],
                             ordered = TRUE)
  kde_shp
}


# ==============================================================================
# FUNCTION: CREATE COMPARISON MAP ----
# ==============================================================================

create_comparison_map <- function(species_key, species_info,
                                  sightings_data, pam_data,
                                  is_grouped = FALSE,
                                  pam_only   = FALSE) {
  
  species_name <- species_info$common_name
  cat("Creating", ifelse(pam_only, "PAM-only", "comparison"),
      "map for:", species_info$display_name, "\n")
  
  # ---- Locate KDE shapefiles ----
  has_sightings_kde <- FALSE
  sightings_df      <- NULL
  
  if (!pam_only && !is.na(species_info$sightings_base)) {
    cat("  Searching sightings KDE base name:", species_info$sightings_base, "\n")
    sightings_df      <- find_kde_shapefiles(species_info$sightings_base,
                                             prefer_proximity = TRUE,
                                             pam_only = FALSE)
    has_sightings_kde <- !is.null(sightings_df) && nrow(sightings_df) > 0
  }
  
  cat("  Searching PAM KDE base name:", species_info$pam_base, "\n")
  pam_df      <- find_kde_shapefiles(species_info$pam_base,
                                     prefer_proximity = FALSE,
                                     pam_only         = TRUE)
  has_pam_kde <- !is.null(pam_df) && nrow(pam_df) > 0
  
  if (!has_sightings_kde && !has_pam_kde) {
    cat("  ERROR: no KDE shapefiles found for", species_info$display_name, "\n\n")
    return(NULL)
  }
  
  if (!has_pam_kde && !CREATE_SIGHTINGS_ONLY_MAPS) {
    cat("  SKIPPING: no PAM data and CREATE_SIGHTINGS_ONLY_MAPS is FALSE\n\n")
    return(NULL)
  }
  
  # ---- Load PAM reference KDE ----
  pam_ref_file <- NA_character_
  pam_ref_kde  <- NULL
  if (has_pam_kde) {
    pam_ref_file <- pam_df$path[1]
    cat("  Using PAM reference KDE:", basename(pam_ref_file), "\n")
    pam_ref_kde  <- load_kde_shapefile(pam_ref_file)
  } else {
    cat("  No PAM KDE found",
        ifelse(CREATE_SIGHTINGS_ONLY_MAPS, " (will create sightings-only map)\n", "\n"))
  }
  
  # Cap sightings variants to avoid runaway loops
  MAX_SIGHTINGS_VARIANTS <- 12
  if (has_sightings_kde && nrow(sightings_df) > MAX_SIGHTINGS_VARIANTS) {
    sightings_df <- sightings_df[seq_len(MAX_SIGHTINGS_VARIANTS), , drop = FALSE]
    cat("  Note: capped sightings variants to", MAX_SIGHTINGS_VARIANTS, "\n")
  }
  
  # ---- Filter source data to species ----
  # group_type controls which records are pulled for each panel.
  group_type <- species_info$group_type %||% NA_character_
  
  if (is_grouped && !pam_only && identical(group_type, "Mysticete")) {
    # All baleen whales
    cat("  Filtering all baleen (Mysticete) records...\n")
    species_sightings <- sightings_data %>% filter(family == "Mysticete")
    species_pam       <- pam_data
    
  } else if (is_grouped && !pam_only && identical(group_type, "Odontocete")) {
    # All toothed whales - sightings only, no PAM
    cat("  Filtering all toothed (Odontocete) records...\n")
    species_sightings <- sightings_data %>% filter(family == "Odontocete")
    species_pam       <- data.frame()
    
  } else if (is_grouped && !pam_only && identical(group_type, "BeakedWhale")) {
    # All beaked whales grouped: sightings pooled from BEAKED_SPECIES_SIGHTINGS,
    # PAM uses all records passed in as pam_data (beaked_pam_data upstream).
    cat("  Filtering grouped beaked whale records...\n")
    species_sightings <- sightings_data %>%
      filter(common_name %in% BEAKED_SPECIES_SIGHTINGS)
    species_pam <- pam_data
    cat("  Beaked sightings records:", nrow(species_sightings), "\n")
    cat("  Species in sightings panel:",
        paste(sort(unique(species_sightings$common_name)), collapse = ", "), "\n")
    
  } else if (pam_only) {
    cat("  PAM-only species - no sightings data used\n")
    species_sightings <- data.frame()
    species_pam <- if (isTRUE(is_grouped) || isTRUE(species_info$is_grouped)) {
      pam_data
    } else if (!is.null(species_info$pam_code)) {
      pam_data %>% dplyr::filter(species == species_info$pam_code)
    } else {
      data.frame()
    }
    
  } else {
    # Individual species - match on common_name
    cat("  Looking for common_name:", species_name, "\n")
    species_sightings <- sightings_data %>% filter(common_name == species_name)
    
    # Try alternate names (e.g. Cuvier's for Goose Beaked)
    if (nrow(species_sightings) == 0 && !is.null(species_info$alternate_sightings_names)) {
      for (alt_name in species_info$alternate_sightings_names) {
        cat("  Trying alternate name:", alt_name, "\n")
        species_sightings <- sightings_data %>% filter(common_name == alt_name)
        if (nrow(species_sightings) > 0) { cat("  Matched with alternate name\n"); break }
      }
    }
    
    if (nrow(species_sightings) == 0) {
      alt1 <- gsub("_", " ", species_name)
      species_sightings <- sightings_data %>% filter(common_name == alt1)
      if (nrow(species_sightings) > 0) cat("  Matched with underscores->spaces:", alt1, "\n")
    }
    
    if (nrow(species_sightings) == 0) {
      alt2 <- gsub("[-_]", " ", species_name)
      species_sightings <- sightings_data %>% filter(common_name == alt2)
      if (nrow(species_sightings) > 0) cat("  Matched with all->spaces:", alt2, "\n")
    }
    
    if (nrow(species_sightings) == 0)
      cat("  WARNING: could not match species in sightings data\n\n")
    
    species_pam <- if (!is.null(species_info$pam_code)) {
      pam_data %>% filter(species == species_info$pam_code)
    } else {
      data.frame()
    }
  }
  
  # ---- Count PAM stations ----
  stations <- if (nrow(species_pam) > 0 && "site" %in% names(species_pam)) {
    nrow(species_pam %>% group_by(site) %>% summarise(n()))
  } else if (nrow(species_pam) > 0 && "station" %in% names(species_pam)) {
    nrow(species_pam %>% group_by(station) %>% summarise(n()))
  } else { 0 }
  
  cat("  Sightings:", nrow(species_sightings), "  PAM Stations:", stations, "\n")
  
  if (nrow(species_sightings) == 0 && stations == 0) {
    cat("  ERROR: no data for sightings or PAM - skipping\n\n")
    return(NULL)
  }
  
  # ---- Helper: build method tag string ----
  make_tag <- function(file) {
    if (is.na(file) || !file.exists(file)) return("none")
    p    <- extract_kde_params(file)
    prox <- ifelse(isTRUE(p$has_proximity), "prox", "unw")
    meth <- ifelse(is.null(p$method) || is.na(p$method), "unknown", p$method)
    bw   <- ifelse(is.na(p$bandwidth), "bwNA", paste0("bw", p$bandwidth))
    paste(prox, meth, bw, sep = "_")
  }
  
  pam_tag         <- make_tag(pam_ref_file)
  sightings_files <- if (!has_sightings_kde || pam_only) character(0) else sightings_df$path
  
  plot_list  <- list()
  n_variants <- max(1, length(sightings_files))
  
  # ---- Loop over sightings KDE variants ----
  for (i in seq_len(n_variants)) {
    
    if (length(sightings_files) >= i) {
      sightings_file <- sightings_files[i]
      s_tag          <- make_tag(sightings_file)
      cat("  Variant", i, "sightings KDE:", basename(sightings_file), "\n")
      sightings_kde  <- load_kde_shapefile(sightings_file)
    } else {
      sightings_file <- NA_character_
      s_tag          <- "none"
      sightings_kde  <- NULL
    }
    
    pam_kde <- pam_ref_kde
    
    if (is.null(sightings_kde) && is.null(pam_kde)) {
      cat("  WARNING: could not load KDEs for variant", i, "\n"); next
    }
    
    # Build named palette from ALL levels BEFORE filtering, so colours are
    # stable -- removing a level just drops its entry, not shifts every colour.
    all_levels <- unique(c(
      if (!is.null(sightings_kde)) as.character(levels(sightings_kde$Quantile)) else character(0),
      if (!is.null(pam_kde))       as.character(levels(pam_kde$Quantile))       else character(0)
    ))
    pal_full <- setNames(
      rev(get(QUANTILE_PALETTE)(length(all_levels) + 1))[seq_along(all_levels)],
      all_levels
    )

    # Drop unwanted contour levels (e.g. outermost 85% which renders near-white)
    if (length(QUANTILE_HIDE) > 0) {
      if (!is.null(sightings_kde))
        sightings_kde <- sightings_kde %>%
          filter(!as.character(Quantile) %in% QUANTILE_HIDE) %>%
          mutate(Quantile = droplevels(Quantile))
      if (!is.null(pam_kde))
        pam_kde <- pam_kde %>%
          filter(!as.character(Quantile) %in% QUANTILE_HIDE) %>%
          mutate(Quantile = droplevels(Quantile))
    }

    # Palette for remaining levels only -- colours unchanged from full set
    pal <- pal_full[!names(pal_full) %in% QUANTILE_HIDE]
    n_q <- length(pal)
    
    base_theme <- theme_bw() +
      theme(
        panel.grid    = element_blank(),
        axis.text     = element_blank(),
        axis.ticks    = element_blank(),
        axis.title    = element_blank(),
        plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position      = c(0.99, 0.02),
        legend.justification = c(1, 0),
        legend.background    = element_rect(fill = "white", color = NA),
        legend.key.height    = unit(0.55, "cm"),
        legend.key.width     = unit(0.30, "cm"),
        legend.title         = element_text(size = 10),
        legend.text          = element_text(size = 9),
        legend.box           = "horizontal"
      )
    
    # ---- Sightings panel ----
    if (!is.null(sightings_kde) && nrow(species_sightings) > 0) {
      
      s_params    <- extract_kde_params(sightings_file)
      bw_text     <- if (!is.na(s_params$bandwidth))
        paste0("BW: ", round(s_params$bandwidth / 1000, 2), "km") else "BW: unknown"
      weight_text <- ifelse(isTRUE(s_params$has_proximity),
                            "Proximity weighted", "Unweighted")
      param_text  <- paste(bw_text, "|", s_params$method, "|", weight_text)
      
      p_sightings <- ggplot() +
        map_layer_bathy_contours(cont) +
        geom_sf(data = sightings_kde, aes(fill = Quantile), color = NA, alpha = 0.6) +
        map_layer_wea(osw_wind, alpha = 0.18, linewidth = 0.8) +
        geom_sf(data = species_sightings,
                color = if (is_grouped) "grey40" else "grey20",
                fill  = if (is_grouped) "grey40" else "grey30",
                alpha = if (is_grouped) 0.5  else 0.7,
                size  = if (is_grouped) 1.0  else 1.5,
                shape = 21) +
        geom_point(data = data.frame(x = NA_real_, y = NA_real_),
                   aes(x = x, y = y, shape = "Sightings"), inherit.aes = FALSE,
                   color = if (is_grouped) "grey40" else "grey20",
                   fill  = if (is_grouped) "grey40" else "grey30",
                   size  = if (is_grouped) 1.0  else 1.5,
                   alpha = if (is_grouped) 0.5  else 0.7) +
        scale_shape_manual(name = "", values = c("Sightings" = 21),
                           guide = guide_legend(override.aes = list(
                             size  = if (is_grouped) 1.0  else 1.5,
                             color = if (is_grouped) "grey40" else "grey20",
                             fill  = if (is_grouped) "grey40" else "grey30",
                             alpha = if (is_grouped) 0.5  else 0.7))) +
        map_layer_land(land) +
        map_layer_study_area(study_area) +
        scale_fill_manual(values = pal, name = "Relative\nOccurrence",
                          drop = FALSE, guide = "none") +
        map_coord(xlim = xlims, ylim = ylims, crs = UTM20) +
        labs(title = "Sightings", subtitle = paste0(scales::comma(nrow(species_sightings)), " sightings records")) +
        base_theme +
        theme(legend.position = c(.97, 0.05)) +
        annotation_scale(location = "br", width_hint = 0.25)
      
    } else {
      p_sightings <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No sightings data", size = 6) +
        theme_void()
    }
    
    # ---- PAM panel ----
    if (!is.null(pam_kde) && stations > 0) {
      
      pam_zero       <- species_pam %>% filter(proportion_det == 0)
      pam_detections <- species_pam %>% filter(proportion_det  > 0)
      
      p_params <- if (!is.na(pam_ref_file)) extract_kde_params(pam_ref_file) else
        list(bandwidth = NA_real_, method = "unknown")
      bw_text  <- if (!is.na(p_params$bandwidth))
        paste0("BW: ", round(p_params$bandwidth / 1000, 2), "km") else "BW: unknown"
      sub_text <- if (is_grouped) {
        paste0(stations, " stations with detections")
      } else {
        paste0(stations, " stations (", nrow(pam_zero), " with 0 detections)")
      }
      
      p_pam <- ggplot() +
        map_layer_bathy_contours(cont) +
        map_layer_wea(osw_wind, alpha = 0.18, linewidth = 0.8) +
        geom_sf(data = pam_kde, aes(fill = Quantile), color = NA, alpha = 0.6) +
        {if (nrow(pam_detections) > 0)
          geom_sf(data = pam_detections,
                  aes(size = proportion_det, color = proportion_det), alpha = 0.8)}
      
      
      if (!is_grouped && nrow(pam_zero) > 0) {
        dummy_zero <- data.frame(x = NA, y = NA, label = "No Validated Detections")
        p_pam <- p_pam +
          geom_sf(data = pam_zero, color = "gray20", fill = "gray80",
                  size = 1.5, shape = 21) +
          geom_point(data = dummy_zero, aes(x = x, y = y, shape = label),
                     color = "gray20", fill = "gray80", size = 1.5) +
          scale_shape_manual(name = "", values = 21,
                             labels = "No Validated\nDetections",
                             guide  = guide_legend(order = 1,
                                                   override.aes = list(size = 3)))
      }
      
      p_pam <- p_pam +
        map_layer_land(land) +
        map_layer_study_area(study_area) +
        scale_fill_manual(values = pal, name = "Relative\nOccurrence",
                          drop = FALSE, guide = guide_legend(order = 2)) +
        scale_size_continuous(name   = "Detection\nProportion",
                              range  = c(1, 4), limits = c(0, 1),
                              breaks = c(0.2, 0.5, 0.9),
                              guide  = guide_legend(order = 3)) +
        scale_color_gradient(low = "orange", high = "darkorange",
                             name   = "Detection\nProportion",
                             limits = c(0, 1), breaks = c(0.2, 0.5, 0.9),
                             guide  = guide_legend(order = 3)) +
        map_coord(xlim = xlims, ylim = ylims, crs = UTM20) +
        labs(title = "PAM", subtitle = sub_text) +
        base_theme
      
    } else {
      p_pam <- NULL
    }
    
    # ---- Combine and save ----
    if (!is.null(pam_kde) && stations > 0 && !pam_only) {
      
      p_combined <- (p_sightings + p_pam) +
        plot_annotation(
          title = paste0(species_info$display_name, " Sightings vs PAM"),
          theme = theme(plot.title = element_text(hjust = 0.5))
        )
      
      out_name <- paste0("compare_",
                         gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
                         "__S_", s_tag, "__P_", pam_tag, ".png")
      
      suppressWarnings(ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 16, height = 8, dpi = 300, bg = "white"))
      
    } else if (pam_only && !is.null(pam_kde) && stations > 0) {
      
      p_combined <- p_pam +
        labs(title = paste0(species_info$display_name, "  (PAM Only)")) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      
      out_name <- paste0("pam_only_",
                         gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
                         "__P_", pam_tag, ".png")
      
      suppressWarnings(ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 10, height = 8, dpi = 300, bg = "white"))
      
    } else {

      p_combined <- p_sightings +
        # Re-enable fill legend (suppressed on sightings panel in paired plots)
        guides(fill = guide_legend(title = "Relative\nOccurrence", order = 1)) +
        # Add a sightings point entry to the legend
        geom_point(data = data.frame(x = NA_real_, y = NA_real_),
                   aes(x = x, y = y, shape = "Sightings"), inherit.aes = FALSE,
                   color = "grey20", fill = "grey30", size = 1.5, alpha = 0.7) +
        scale_shape_manual(name = "", values = c("Sightings" = 21),
                           guide = guide_legend(order = 2,
                                                override.aes = list(size = 1.5, color = "grey20",
                                                                     fill = "grey30", alpha = 0.7))) +
        labs(title = paste0(species_info$display_name, " Sightings")) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              legend.position = c(0.99, 0.05))
      
      out_name <- paste0("sightings_only_",
                         gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
                         "__S_", s_tag, ".png")
      
      suppressWarnings(ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 10, height = 8, dpi = 300, bg = "white"))
    }
    
    cat("  Saved:", out_name, "\n\n")
    plot_list[[length(plot_list) + 1]] <- p_combined
  }
  
  if (length(plot_list) == 0) {
    cat("  WARNING: no plots generated for", species_info$display_name, "\n\n")
    return(NULL)
  }
  
  if (length(plot_list) > 1) patchwork::wrap_plots(plot_list, ncol = 1) else plot_list[[1]]
}


# ==============================================================================
# HELPER: DISCOVER SIGHTINGS-ONLY SPECIES FROM KDE FILES ----
# ==============================================================================

discover_sightings_only_species <- function(shapefile_dir   = SHAPEFILE_DIR,
                                            exclude_species = character(),
                                            sightings_data  = NULL) {
  
  cat("Discovering sightings-only species from KDE files...\n")
  
  all_kde_files <- list.files(shapefile_dir, pattern = "^KDE_.*\\.shp$",
                              full.names = FALSE, ignore.case = TRUE)
  
  # Only auto-discover individual proximity-weighted sightings KDEs
  sightings_kdes <- all_kde_files[
    grepl("^KDE_sightings_species_proximity_", all_kde_files, ignore.case = TRUE)
  ]
  
  if (length(sightings_kdes) == 0) {
    cat("  No sightings KDE files found\n")
    return(list())
  }
  
  species_list <- list()
  
  for (file in sightings_kdes) {
    
    # Strip prefix and method/bandwidth suffixes to isolate the species token
    temp         <- sub("^KDE_sightings_species_proximity_", "", file, ignore.case = TRUE)
    temp         <- sub("\\.shp$", "", temp, ignore.case = TRUE)
    temp         <- sub("_(fixed|adaptive|diggle).*", "", temp, ignore.case = TRUE)
    temp         <- sub("_bw\\d+.*", "", temp, ignore.case = TRUE)
    species_base <- sub("_+$", "", temp)
    
    if (species_base %in% exclude_species) next
    if (species_base %in% names(species_list)) next
    
    # Look up correct common_name directly from sightings data, stripping
    # apostrophes so e.g. "RISSOS" correctly matches "Risso's Dolphin".
    epithet <- sub("^(WHALE|DOLPHINS|PORPOISE)-", "", species_base, ignore.case = TRUE)
    
    matched_name <- NULL
    if (!is.null(sightings_data)) {
      epithet_clean <- gsub("[-_']", "", tolower(epithet))
      matched_name <- sightings_data %>%
        st_drop_geometry() %>%
        filter(grepl(epithet_clean,
                     gsub("[-_' ]", "", tolower(common_name)),
                     ignore.case = TRUE)) %>%
        pull(common_name) %>%
        unique()
      if (length(matched_name) != 1) matched_name <- NULL
    }
    
    display_name <- if (!is.null(matched_name)) {
      matched_name
    } else {
      fallback <- tools::toTitleCase(tolower(gsub("[-_]", " ", epithet)))
      cat("  WARNING: could not uniquely match", species_base,
          "in sightings data - using fallback:", fallback, "\n")
      fallback
    }
    
    species_list[[species_base]] <- list(
      common_name       = display_name,
      display_name      = display_name,
      sightings_base    = species_base,
      pam_base          = NA_character_,
      is_sightings_only = TRUE,
      is_grouped        = FALSE
    )
    
    cat("  Found:", display_name, "(base token:", species_base, ")\n")
  }
  
  cat("  Total sightings-only species found:", length(species_list), "\n\n")
  species_list
}


# ==============================================================================
# RUN: BALEEN WHALES ----
# ==============================================================================

cat("=== CREATING  SIGHTINGS vs PAM COMPARISON MAPS ===\n")
cat("CREATE_SIGHTINGS_ONLY_MAPS =", CREATE_SIGHTINGS_ONLY_MAPS, "\n\n")

cat("--- Processing baleen whales (with PAM data) ---\n\n")
for (species_key in names(BALEEN_SPECIES_MAPPING)) {
  create_comparison_map(species_key, BALEEN_SPECIES_MAPPING[[species_key]],
                        sightings_data, baleen_pam_data, is_grouped = FALSE)
}


# ==============================================================================
# RUN: BEAKED WHALES ----
# ==============================================================================

cat("\n--- Processing beaked whales with sightings ---\n\n")
for (species_key in names(BEAKED_SPECIES_MAPPING)) {
  sp <- BEAKED_SPECIES_MAPPING[[species_key]]
  if (isTRUE(sp$has_sightings))
    create_comparison_map(species_key, sp, sightings_data, beaked_pam_data,
                          is_grouped = FALSE, pam_only = FALSE)
}

cat("\n--- Processing PAM-only beaked whales ---\n\n")
for (species_key in names(BEAKED_SPECIES_MAPPING)) {
  sp <- BEAKED_SPECIES_MAPPING[[species_key]]
  if (!isTRUE(sp$has_sightings))
    create_comparison_map(species_key, sp, sightings_data, beaked_pam_data,
                          is_grouped = FALSE, pam_only = TRUE)
}


# ==============================================================================
# RUN: SIGHTINGS-ONLY SPECIES ----
# ==============================================================================

if (CREATE_SIGHTINGS_ONLY_MAPS) {
  cat("\n--- Processing sightings-only species ---\n\n")
  
  exclude_species <- c(
    sapply(BALEEN_SPECIES_MAPPING, `[[`, "sightings_base"),
    sapply(BEAKED_SPECIES_MAPPING[
      sapply(BEAKED_SPECIES_MAPPING, `[[`, "has_sightings")],
      `[[`, "sightings_base")
  )
  
  sightings_only_species <- discover_sightings_only_species(
    exclude_species = exclude_species,
    sightings_data  = sightings_data
  )
  
  if (length(sightings_only_species) > 0) {
    cat("Processing", length(sightings_only_species), "sightings-only species...\n\n")
    for (species_key in names(sightings_only_species)) {
      create_comparison_map(species_key, sightings_only_species[[species_key]],
                            sightings_data, pam_data, is_grouped = FALSE)
    }
  }
}


# ==============================================================================
# RUN: GROUPED MAPS ----
# ==============================================================================

cat("\n--- Processing grouped Odontocetes (sightings) ---\n\n")
for (species_key in names(GROUPED_ODONTOCETES_SIGHTINGS)) {
  create_comparison_map(species_key, GROUPED_ODONTOCETES_SIGHTINGS[[species_key]],
                        sightings_data, pam_data, is_grouped = TRUE, pam_only = FALSE)
}

cat("\n--- Processing grouped all baleen whales ---\n\n")
for (species_key in names(GROUPED_BALEEN)) {
  create_comparison_map(species_key, GROUPED_BALEEN[[species_key]],
                        sightings_data, baleen_pam_data, is_grouped = TRUE)
}

# Updated from v1: now a  SIGHTINGS vs PAM comparison map.
# Sightings panel: records matching BEAKED_SPECIES_SIGHTINGS from sightings_data.
# PAM panel: all beaked_pam_data records (already filtered to beaked species).
cat("\n--- Processing grouped all beaked whales ( Sightings vs PAM) ---\n\n")
for (species_key in names(GROUPED_BEAKED)) {
  create_comparison_map(species_key, GROUPED_BEAKED[[species_key]],
                        sightings_data, beaked_pam_data,
                        is_grouped = TRUE, pam_only = FALSE)
}

cat("\n=== COMPARISON MAPPING COMPLETE ===\n")
cat("Maps saved to:", COMPARE_OUTPUT_DIR, "\n")


# ==============================================================================
# EXPORT PARAMETERS FOR REPORT ----
# ==============================================================================

compare_params <- tibble::tribble(
  ~Parameter,                 ~Value,                                            ~Script,                     ~Notes,
  "OSW_ALPHA",                as.character(OSW_ALPHA),                           "10_compare_kde_pam_vs_sightings.R", "OSW polygon overlay transparency",
  "POINT_COLOR",              as.character(POINT_COLOR),                         "10_compare_kde_pam_vs_sightings.R", "Sightings point colour on comparison maps",
  "PAM_BUFFER_DIST_M",        as.character(PAM_BUFFER_DIST_M),                   "10_compare_kde_pam_vs_sightings.R", "PAM station inclusion buffer (m)",
  "BEAKED_SPECIES_SIGHTINGS", paste(BEAKED_SPECIES_SIGHTINGS, collapse = "; "),  "10_compare_kde_pam_vs_sightings.R", "Common names pooled for grouped beaked sightings panel"
)
readr::write_csv(compare_params, here::here("output/data/params_compare.csv"))
cat("+ Params exported to output/data/params_compare.csv\n")

