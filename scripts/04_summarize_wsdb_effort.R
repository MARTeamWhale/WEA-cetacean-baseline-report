# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-03-12
# Script: 04_summarize_wsdb_effort.R
# Description: Generates gridded maps of cetacean sightings from opportunistic
#              vessel survey data. Produces: (1) an all-cetacean coverage map,
#              (2) a data confidence/consistency map, (3) raw-count and
#              effort-normalised maps for each target species group, and
#              (4) individual species maps. Also outputs a species-level
#              appendix table (XLSX) and a seasonal summary table comparing
#              WEA cells vs. the broader study region for priority species and
#              guilds.
#
# Changes from previous version:
#     This resolves near-boundary species showing artificially
#     low WEA record counts in downstream confidence summaries.
#   - Added seasonal summary table (Section 15): records by season x region
#     for target species and guilds, saved as CSV and printed to console.
# ==============================================================================

suppressWarnings(source(here::here("scripts/00_load_helpers.R")))

# ---- 1. Configuration -----------------------------------------------------------

# Performance options
use_datatable <- TRUE     # Faster aggregation with data.table
use_parallel  <- TRUE     # Parallel processing for multiple targets

# Toggle survey area overlay (5 km WEA buffer shapefile)
# Set TRUE only if the WEA_5km_buffer shapefile is relevant to your project.
# When FALSE, in_survey tagging uses the WEA polygons directly (WEA_wind).
use_survey_area <- FALSE

# Input / output paths
infile  <- "output/data/combined_dedup_1km_day.csv"
out_dir <- "output/figs/Effort_maps/"

# Grid and filtering parameters
grid_km           <- 25      # Grid cell size in kilometres
min_records_all   <- 5       # Minimum records to show normalised maps
year_min          <- 2010    # Start year
year_max          <- 2025    # End year
bbox_lonlat       <- NULL    # Optional spatial crop: c(xmin, xmax, ymin, ymax)
log_transform_all <- TRUE    # Log-transform coverage map
weight_col        <- NULL    # Record weighting column (NULL = 1 per row)

# Spatial reference files
WEA_polygons_path    <- spatial_path("wea")
survey_polygons_path <- spatial_path("wea_5km_buffer")
study_area_path      <- spatial_path("study_area")
land_path            <- spatial_path("land")

# Map styling
WEA_color    <- "#FF6B35"
WEA_alpha    <- 0.3
WEA_linewidth <- 1
map_width    <- 10
map_height   <- 8
map_dpi      <- 300

# Seasonal definitions
spring_months <- c(3, 4, 5)    # March, April, May
summer_months <- c(6, 7, 8)    # June, July, August
fall_months   <- c(9, 10, 11)  # September, October, November
winter_months <- c(12, 1, 2)   # December, January, February

# Target species groups (regex patterns)
targets <- list(
  NARW              = "Eubalaena\\s+glacialis|RIGHT\\s+WHALE|N\\s*ATLANTIC\\s+RIGHT",
  Blue_whale        = "Balaenoptera\\s+musculus|BLUE\\s+WHALE",
  Shelf_baleen      = "Megaptera|HUMPBACK|acutorostrata|MINKE|physalus|FIN\\s+WHALE|borealis|SEI|musculus|BLUE|glacialis|RIGHT|RORQUAL|BALEEN",
  Harbour_porpoise  = "Phocoena\\s+phocoena|HARBOU?R\\s+PORPOISE",
  Small_odontocetes = "DOLPHIN|Delphin|Stenella|Tursiops|Lagenorhynchus|Grampus|WHITE[-\\s]?SIDED|WHITE[-\\s]?|PILOT|Phocoena\\s+phocoena|HARBOU?R\\s+PORPOISE",
  Deep_divers       = "Mesoplodon|Hyperoodon|Ziphius|BEAKED\\s+WHALE|SOWERBY|CUVIER|BOTTLENOSE"
)

# Human-readable titles
target_titles <- list(
  NARW              = "North Atlantic Right Whale",
  Blue_whale        = "Blue Whale",
  Shelf_baleen      = "Baleen Whales",
  Harbour_porpoise  = "Harbour Porpoise",
  Small_odontocetes = "Small Odontocetes",
  Deep_divers       = "Deep Divers"
)

# Viridis palette per target
target_palettes <- list(
  NARW              = "inferno",
  Blue_whale        = "cividis",
  Shelf_baleen      = "magma",
  Harbour_porpoise  = "turbo",
  Small_odontocetes = "plasma",
  Deep_divers       = "rocket"
)

# Species-level map settings
species_field       <- "common_name"
appendix_xlsx       <- file.path(out_dir, "APPENDIX_species_summary.xlsx")
make_species_maps   <- TRUE
top_n_species_maps  <- 30
species_map_dir     <- file.path(out_dir, "species_maps")
species_map_bins    <- c(0, 5, 25, 50, 100, Inf)
species_map_labels  <- c("1-5", "6-25", "26-50", "51-100", "100+")
exclude_zero_WEA    <- TRUE
exclude_beaked_whales <- TRUE


# ---- 2. Packages ----------------------------------------------------------------

setup_packages <- function() {
  pkgs <- c("dplyr", "readr", "lubridate", "ggplot2", "stringr",
            "tidyr", "sf", "ggspatial", "ggpattern", "writexl", "png", "here")
  if (exists("use_datatable") && use_datatable) pkgs <- c(pkgs, "data.table")
  if (exists("use_parallel")  && use_parallel)  pkgs <- c(pkgs, "future", "furrr")
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

setup_packages()

# Resolve conflicts explicitly -- parallel workers do not inherit session
# conflict preferences so these must be declared before furrr::future_walk.
conflicted::conflicts_prefer(lubridate::year,    .quiet = TRUE)
conflicted::conflicts_prefer(lubridate::month,   .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::filter,      .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::select,      .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::mutate,      .quiet = TRUE)
conflicted::conflicts_prefer(dplyr::between,     .quiet = TRUE)


# ---- 3. Helper functions --------------------------------------------------------

# Safe filename string
nm_safe <- function(x) stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")

# Longest consecutive year streak (temporal consistency metric)
longest_streak <- function(years_vec) {
  y <- sort(unique(years_vec[!is.na(years_vec)]))
  if (length(y) == 0) return(0L)
  grp <- cumsum(c(TRUE, diff(y) != 1))
  max(as.integer(table(grp)))
}

# Standardised map theme
theme_map <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text    = element_blank(),
      axis.ticks   = element_blank(),
      axis.title   = element_blank(),
      plot.title    = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption  = element_text(size = 8, hjust = 0)
    )
}

# Add standard spatial overlays (land, WEA boundary, study area, scale bar)
# use_survey_overlay: if TRUE and survey_area is not NULL, draws the 5 km buffer
#                     boundary instead of the WEA polygons (species maps only)
add_spatial_layers <- function(p, land, WEA_wind, study_area,
                               survey_area = NULL, use_survey_overlay = FALSE) {
  if (!is.null(land)) {
    p <- p + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
  }
  if (use_survey_overlay && !is.null(survey_area)) {
    p <- p + geom_sf(data = survey_area, fill = NA, color = "blue",
                     linetype = "solid", alpha = 0.7,
                     linewidth = 0.6, inherit.aes = FALSE)
  } else if (!is.null(WEA_wind)) {
    p <- p + geom_sf(data = WEA_wind, fill = NA, color = WEA_color,
                     linetype = "dashed", alpha = WEA_alpha,
                     linewidth = WEA_linewidth, inherit.aes = FALSE)
  }
  if (!is.null(study_area)) {
    p <- p + geom_sf(data = study_area, fill = NA, color = "black",
                     linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
  }
  p <- p + ggspatial::annotation_scale(location = "br", width_hint = 0.25)
  return(p)
}


# ---- 4. Load and clean data -----------------------------------------------------

message("\n=== LOADING DATA ===")
message("Reading: ", infile)

required_cols <- c("lon", "lat", "date_utc", "scientific_name", "common_name")
optional_cols <- c("is_duplicate", if (!is.null(weight_col)) weight_col else NULL)

df <- readr::read_csv(
  infile,
  col_select    = all_of(c(required_cols, optional_cols)),
  show_col_types = FALSE
)

message("Initial rows: ", nrow(df))

df <- df %>%
  mutate(
    lon            = as.numeric(lon),
    lat            = as.numeric(lat),
    date_utc       = lubridate::as_date(date_utc),
    year           = lubridate::year(date_utc),
    month          = lubridate::month(date_utc),
    common_name    = as.character(common_name),
    scientific_name = as.character(scientific_name),
    is_duplicate   = if ("is_duplicate" %in% names(.)) as.logical(is_duplicate) else FALSE
  ) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(date_utc),
         year >= year_min, year <= year_max)

if (!is.null(bbox_lonlat)) {
  df <- df %>%
    filter(!is.na(lon), !is.na(lat),
           between(lon, bbox_lonlat["xmin"], bbox_lonlat["xmax"]),
           between(lat, bbox_lonlat["ymin"], bbox_lonlat["ymax"]))
}

message("Rows after filtering: ", nrow(df))

# Weighting
df <- df %>%
  mutate(
    w = if (!is.null(weight_col) && weight_col %in% names(.)) {
      as.numeric(.data[[weight_col]])
    } else { 1 },
    w = ifelse(is.na(w) | w <= 0, 1, w)
  )

# Name blob for target regex matching
df <- df %>%
  mutate(
    name_blob = paste0(scientific_name, " | ", common_name),
    # Season assignment
    season = case_when(
      month %in% spring_months ~ "Spring",
      month %in% summer_months ~ "Summer",
      month %in% fall_months   ~ "Fall",
      month %in% winter_months ~ "Winter",
      TRUE ~ NA_character_
    ),
    is_spring = month %in% spring_months
  )

# Aggregate ambiguous species groups
df <- df %>%
  mutate(
    common_name = case_when(
      grepl("SEAL", common_name, ignore.case = TRUE) &
        !grepl("SEA LION|FUR SEAL", common_name, ignore.case = TRUE) ~ "All Seals",
      grepl("WHALE[-\\s]*FIN[-/\\s]*SEI|WHALE[-\\s]*SEI(?!.*FIN)",
            common_name, ignore.case = TRUE, perl = TRUE)             ~ "Fin/Sei Whale",
      grepl(paste0("WHALE.*\\(NS\\)|\\(NS\\).*WHALE|CETACEAN.*\\(NS\\)|",
                   "\\(NS\\).*CETACEAN|UNIDENTIFIED.*WHALE|WHALE.*UNIDENTIFIED|",
                   "UNIDENTIFIED.*CETACEAN|CETACEAN.*UNIDENTIFIED|",
                   "WHALE.*NOT\\s+SPECIFIED|CETACEAN.*NOT\\s+SPECIFIED"),
            common_name, ignore.case = TRUE)                          ~ "Unidentified Whale",
      TRUE ~ common_name
    ),
    scientific_name = case_when(
      common_name == "All Seals"         ~ "Phocidae spp.",
      common_name == "Fin/Sei Whale"     ~ "Balaenoptera physalus/borealis",
      common_name == "Unidentified Whale" ~ "Cetacea spp.",
      TRUE ~ scientific_name
    ),
    name_blob = paste0(scientific_name, " | ", common_name)
  )

message("Species aggregated: Seals, Fin/Sei Whale, Unidentified Whale")


# ---- 5. Grid calculation --------------------------------------------------------

message("\n=== CALCULATING GRID CELLS ===")
cell_m <- grid_km * 1000

# Project lon/lat to UTM zone 20N for metre-based grid calculation
df_sf <- df %>%
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  sf::st_transform(crs = 32620)

coords_m <- sf::st_coordinates(df_sf)

df <- df %>%
  mutate(
    x_m     = coords_m[, 1],
    y_m     = coords_m[, 2],
    gx      = floor(x_m / cell_m),
    gy      = floor(y_m / cell_m),
    cell_id = paste0(gx, "_", gy),
    x0      = gx * cell_m,
    y0      = gy * cell_m,
    xc      = x0 + cell_m / 2,
    yc      = y0 + cell_m / 2
  )

rm(df_sf, coords_m)

# ---- 6. Load spatial layers -----------------------------------------------------

message("\n=== LOADING SPATIAL LAYERS ===")

message("  Loading WEA polygons (Designated_WEAs) ...")
WEA_wind    <- load_spatial_layer("wea")
survey_area <- if (use_survey_area) load_spatial_layer("wea_5km_buffer", required = FALSE) else NULL
message("  Loading study area ...")
study_area  <- load_spatial_layer("study_area")
message("  Loading land polygons ...")
land        <- load_spatial_layer("land", required = FALSE)

# Total study area in km²
study_area_km2 <- as.numeric(sf::st_area(sf::st_union(study_area))) / 1e6
WEA_wind_km2 <- as.numeric(sf::st_area(sf::st_union(WEA_wind))) / 1e6

WEA_wind_km2/study_area_km2

# ---- 7. Spatial tagging ---------------------------------------------------------
# in_survey: TRUE for any record whose grid cell overlaps the WEA polygon by at
# least wea_overlap_threshold (25%). This matches the cell-selection logic used
# in 05_summarize_wsdb_confidence.R, ensuring species appendix counts, seasonal
# summaries, and confidence tables are all consistent.
#
# Cell-level tagging is more appropriate than point-in-polygon at 25 km grid
# resolution: a sighting 2 km outside the WEA boundary is ecologically within
# the WEA footprint at this scale. The 25% threshold excludes cells that only
# clip the WEA corner (confirmed: edge cell 24_199 = 0.18% overlap, excluded;
# genuine WEA cell 27_198 = 25.1% overlap, retained).

wea_overlap_threshold <- 0.25   # Must match 05_summarize_wsdb_confidence.R

message("\n=== TAGGING RECORDS BY WEA CELL (>=", wea_overlap_threshold * 100, "% overlap) ===")

# Build temporary cell polygons for overlap calculation
# (cell_all exists after Section 9; tagging is deferred to after aggregation)
# NOTE: df is converted to dt after tagging below -- keep this section before
#       the data.table conversion.

if (!is.null(WEA_wind)) {
  # Build cell polygons from df grid coordinates
  cells_for_tagging <- df %>%
    distinct(cell_id, xc, yc) %>%
    mutate(
      geometry = purrr::map2(xc, yc, function(x, y) {
        half <- cell_m / 2
        sf::st_polygon(list(matrix(
          c(x - half, y - half,
            x + half, y - half,
            x + half, y + half,
            x - half, y + half,
            x - half, y - half),
          ncol = 2, byrow = TRUE
        )))
      })
    ) %>%
    sf::st_as_sf(crs = 32620)
  
  wea_union  <- sf::st_union(WEA_wind)
  cell_area  <- as.numeric(cell_m^2)
  
  inter <- sf::st_intersection(cells_for_tagging, wea_union) %>%
    mutate(overlap_pct = as.numeric(sf::st_area(.)) / cell_area) %>%
    sf::st_drop_geometry() %>%
    select(cell_id, overlap_pct)
  
  wea_cell_ids <- inter %>%
    filter(overlap_pct >= wea_overlap_threshold) %>%
    pull(cell_id)
  
  df$in_survey <- df$cell_id %in% wea_cell_ids
  
  message(sprintf("  Cells intersecting WEA:         %d", nrow(inter)))
  message(sprintf("  Cells meeting >=%.0f%% threshold: %d",
                  wea_overlap_threshold * 100, length(wea_cell_ids)))
  message(sprintf("  Records tagged in WEA cells:    %d", sum(df$in_survey)))
  
} else {
  df$in_survey <- FALSE
  message("  WARNING: WEA_wind is NULL -- all in_survey set FALSE.")
}

# NARW diagnostic
narw_n <- sum(grepl("RIGHT|glacialis", df$name_blob, ignore.case = TRUE) & df$in_survey)
message("  NARW records in WEA cells: ", narw_n)

if (use_datatable) {
  message("\n=== Converting to data.table ===")
  dt <- data.table::as.data.table(df)
  data.table::setkey(dt, cell_id)
} else {
  dt <- df
}


# ---- 8. Map extent --------------------------------------------------------------

message("\n=== Calculating map extent ===")

xlims <- range(dt$xc, na.rm = TRUE) + c(-cell_m / 2, cell_m / 2)
ylims <- range(dt$yc, na.rm = TRUE) + c(-cell_m / 2, cell_m / 2)
message(sprintf("  X: [%.0f, %.0f] | Y: [%.0f, %.0f]",
                xlims[1], xlims[2], ylims[1], ylims[2]))

if (!is.null(land)) {
  bbox_crop <- sf::st_bbox(
    c(xmin = xlims[1], xmax = xlims[2], ymin = ylims[1], ymax = ylims[2]),
    crs = sf::st_crs(32620)
  )
  land <- sf::st_crop(land, bbox_crop)
}


# ---- 9. Aggregate all cetacean records ------------------------------------------

message("\n=== AGGREGATING ALL CETACEAN RECORDS ===")

if (use_datatable && data.table::is.data.table(dt)) {
  cell_all <- dt[, .(
    n_all   = sum(w, na.rm = TRUE),
    n_years = data.table::uniqueN(year, na.rm = TRUE)
  ), by = .(cell_id, gx, gy, xc, yc)] %>% as_tibble()
} else {
  cell_all <- dt %>%
    group_by(cell_id, gx, gy, xc, yc) %>%
    summarise(n_all   = sum(w, na.rm = TRUE),
              n_years = n_distinct(year, na.rm = TRUE),
              .groups = "drop")
}


# ---- 10. Calculate confidence index ---------------------------------------------
# Confidence combines two dimensions:
#   Coverage:  total cetacean records per cell (low <5, medium 5-12, high >12)
#   Temporal:  distinct years recorded      (low <3, medium 3-5, high >5)
# Thresholds are derived from the 33rd and 67th percentiles of cells with data.

message("\n=== CALCULATING DATA CONFIDENCE ===")

cells_with_data <- cell_all %>% filter(n_all > 0)

n_all_q33   <- quantile(cells_with_data$n_all,   0.33, na.rm = TRUE)
n_all_q67   <- quantile(cells_with_data$n_all,   0.67, na.rm = TRUE)
n_years_q33 <- quantile(cells_with_data$n_years, 0.33, na.rm = TRUE)
n_years_q67 <- quantile(cells_with_data$n_years, 0.67, na.rm = TRUE)

message(sprintf("  Coverage: Low < %.0f | Med %.0f-%.0f | High > %.0f",
                n_all_q33, n_all_q33, n_all_q67, n_all_q67))
message(sprintf("  Temporal: Low < %.0f yrs | Med %.0f-%.0f yrs | High > %.0f yrs",
                n_years_q33, n_years_q33, n_years_q67, n_years_q67))

cell_all <- cell_all %>%
  mutate(
    coverage_cat = case_when(
      n_all == 0            ~ "None",
      n_all < n_all_q33     ~ "Low",
      n_all < n_all_q67     ~ "Medium",
      TRUE                  ~ "High"
    ),
    temporal_cat = case_when(
      n_years == 0          ~ "None",
      n_years < n_years_q33 ~ "Low",
      n_years < n_years_q67 ~ "Medium",
      TRUE                  ~ "High"
    ),
    confidence = case_when(
      n_all == 0 | n_years == 0                              ~ "No Data",
      coverage_cat == "High"   & temporal_cat == "High"     ~ "High",
      coverage_cat == "High"   | temporal_cat == "High"     ~ "Medium-High",
      coverage_cat == "Medium" & temporal_cat == "Medium"   ~ "Medium",
      coverage_cat == "Medium" | temporal_cat == "Medium"   ~ "Low-Medium",
      TRUE                                                   ~ "Low"
    ),
    confidence = factor(confidence,
                        levels = c("High", "Medium-High", "Medium", "Low-Medium", "Low"))
  )


# ---- 11. Coverage map -----------------------------------------------------------

message("\n=== CREATING COVERAGE MAP ===")

if (log_transform_all) {
  cell_all   <- cell_all %>% mutate(n_all_display = log10(n_all + 1))
  log_breaks <- seq(0, ceiling(max(cell_all$n_all_display, na.rm = TRUE)), by = 1)
  actual_values        <- 10^log_breaks
  actual_values[1]     <- 0
} else {
  cell_all      <- cell_all %>% mutate(n_all_display = n_all)
  log_breaks    <- NULL
  actual_values <- NULL
}

p_all <- ggplot(cell_all) +
  geom_tile(aes(x = xc, y = yc, fill = n_all_display)) +
  {if (log_transform_all) {
    scale_fill_viridis_c(option = "viridis", na.value = "grey90", name = "Count",
                         breaks = log_breaks, labels = actual_values,
                         begin = 0.15, end = 0.95)
  } else {
    scale_fill_viridis_c(option = "viridis", na.value = "grey90", name = "Count",
                         begin = 0.15, end = 0.95)
  }} +
  labs(title = "Density of all cetacean sightings") +
  annotate("text", x = -Inf, y = Inf,
           label = paste0(year_min, "-", year_max, " | ", grid_km, " km grid",
                          if (log_transform_all) " | log scale" else ""),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_map()

p_all <- suppressMessages(add_spatial_layers(p_all, land, WEA_wind, study_area, survey_area))
suppressMessages(ggsave(file.path(out_dir, "00_all_cetacean_records.png"),
       p_all, width = map_width, height = map_height, dpi = map_dpi))
message("  Saved: 00_all_cetacean_records.png")


# ---- 12. Confidence map ---------------------------------------------------------

message("\n=== CREATING CONFIDENCE MAP ===")

confidence_colors <- c(
  "Low"         = "#fee5d9",
  "Low-Medium"  = "#fdd49e",
  "Medium"      = "#fdbb84",
  "Medium-High" = "#a1d99b",
  "High"        = "#31a354"
)

cell_all_pattern <- cell_all %>%
  mutate(low_coverage = confidence %in% c("No Data", "Low", "Low-Medium"))

p_conf <- ggplot(cell_all_pattern) +
  ggpattern::geom_tile_pattern(
    aes(x = xc, y = yc, fill = confidence,
        pattern = ifelse(low_coverage, "stripe", "none")),
    pattern_fill = "white", pattern_color = "white",
    width = cell_m, height = cell_m,
    pattern_density = 0.02, pattern_spacing = 0.015,
    pattern_angle = 45, pattern_alpha = 0.5
  ) +
  scale_fill_manual(values = confidence_colors, name = "Data coverage", drop = FALSE) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(
    pattern = c("none", "none", "none", "stripe", "stripe"),
    pattern_fill = "white", pattern_color = "white"
  ))) +
  labs(
    title   = "Cetacean sighting coverage and temporal consistency",
    caption = sprintf(
      "Coverage: Low < %.0f, Med %.0f-%.0f, High > %.0f records | Temporal: Low < %.0f, Med %.0f-%.0f, High > %.0f years",
      n_all_q33, n_all_q33, n_all_q67, n_all_q67,
      n_years_q33, n_years_q33, n_years_q67, n_years_q67
    )
  ) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0(year_min, "-", year_max, " | ", grid_km, " km grid"),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_map()

p_conf <- suppressMessages(add_spatial_layers(p_conf, land, WEA_wind, study_area, survey_area))
suppressMessages(ggsave(file.path(out_dir, "00_data_consistency.png"),
       p_conf, width = map_width, height = map_height, dpi = map_dpi))
message("  Saved: 00_data_consistency.png")


# ---- 13. Target species map function --------------------------------------------
# Produces two maps per target group:
#   Map 1: Raw sightings count per cell (reflects effort + distribution)
#   Map 2: Effort-normalised share = target / all cetaceans per cell
#          Controls for uneven observer effort across the study area.
#          Cells with < min_records_all total cetacean sightings are masked.

create_target_maps <- function(nm, target_regex, dt_input, cell_all_conf, use_dt = FALSE) {
  
  message("\n=== Processing: ", nm, " ===")
  
  if (use_dt && data.table::is.data.table(dt_input)) {
    cell_sum <- dt_input[, .(
      n_all    = sum(w, na.rm = TRUE),
      n_target = sum(w * stringr::str_detect(name_blob, target_regex), na.rm = TRUE)
    ), by = .(cell_id, gx, gy, x0, y0, xc, yc)] %>% as_tibble()
    
    target_years <- dt_input[
      stringr::str_detect(name_blob, target_regex),
      .(n_target_years = data.table::uniqueN(year, na.rm = TRUE)),
      by = .(cell_id)
    ] %>% as_tibble()
    
  } else {
    cell_sum <- dt_input %>%
      mutate(is_target = stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id, gx, gy, xc, yc) %>%
      summarise(n_all    = sum(w, na.rm = TRUE),
                n_target = sum(w[is_target], na.rm = TRUE),
                .groups  = "drop")
    
    target_years <- dt_input %>%
      filter(stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(n_target_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
  }
  
  cell_sum <- cell_sum %>%
    left_join(target_years, by = "cell_id") %>%
    mutate(
      n_target_years = tidyr::replace_na(n_target_years, 0),
      share_target   = ifelse(
        n_all >= min_records_all & n_all > 0,
        n_target / n_all,
        NA_real_
      )
    )
  
  # Evidence classification for confidence overlay
  cell_sum <- cell_sum %>%
    mutate(
      evidence = case_when(
        n_target == 0                         ~ "None",
        n_target >= 5 & n_target_years >= 3  ~ "High",
        n_target >= 2 & n_target_years >= 2  ~ "Medium",
        TRUE                                  ~ "Low"
      ),
      evidence = factor(evidence, levels = c("High", "Medium", "Low", "None"))
    ) %>%
    left_join(cell_all_conf %>% select(cell_id, confidence), by = "cell_id") %>%
    mutate(
      target_conf = case_when(
        confidence %in% c("No Data", "Low") & evidence %in% c("None", "Low") ~ "Low (both)",
        evidence == "None"                                                     ~ "No evidence",
        TRUE                                                                   ~ "Acceptable"
      ),
      target_conf  = factor(target_conf, levels = c("Acceptable", "No evidence", "Low (both)")),
      low_coverage = target_conf %in% c("No evidence", "Low (both)")
    )
  
  nm_file <- nm_safe(nm)
  
  # Map 1: raw counts
  p1 <- ggplot(cell_sum) +
    geom_tile(aes(x = xc, y = yc, fill = n_target)) +
    scale_fill_viridis_c(option = target_palettes[[nm]], na.value = "grey90",
                         name = "Count", begin = 0.10, end = 0.95) +
    labs(title   = paste0("Opportunistic sightings: ", target_titles[[nm]]),
         caption = "Raw counts reflect both species distribution and observer effort.") +
    annotate("text", x = -Inf, y = Inf,
             label = paste0(year_min, "-", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    theme_map()
  p1 <- suppressMessages(add_spatial_layers(p1, land, WEA_wind, study_area, survey_area))
  suppressMessages(ggsave(file.path(out_dir, paste0(nm_file, "_1_target_records.png")),
         p1, width = map_width, height = map_height, dpi = map_dpi))
  
  # Map 2: effort-normalised share with hatching for low-confidence cells
  low_conf_cells <- cell_sum %>% filter(low_coverage)
  
  p2 <- ggplot(cell_sum) +
    ggpattern::geom_tile_pattern(
      aes(x = xc, y = yc, fill = share_target,
          pattern = ifelse(low_coverage, "stripe", "none")),
      pattern_fill = "white", pattern_color = "white",
      width = cell_m, height = cell_m,
      pattern_density = 0.03, pattern_spacing = 0.015,
      pattern_angle = 45, pattern_alpha = 1.0, pattern_linewidth = 0.5
    ) +
    scale_fill_viridis_c(option = "viridis", na.value = "grey90",
                         limits = c(0, 1), name = "Proportion of\nall sightings",
                         direction = 1, begin = 0.1, end = 0.95) +
    scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
    {if (nrow(low_conf_cells) > 0)
      geom_tile(data = low_conf_cells, aes(x = xc, y = yc),
                fill = NA, color = "white", linewidth = 0.6)} +
    labs(title   = paste0("Effort-normalised sightings: ", target_titles[[nm]]),
         caption = paste0(
           "Each cell: (target sightings) / (all cetacean sightings). Controls for observer effort.\n",
           "White hatching: unreliable proportion (poor coverage or zero target records).\n",
           "Grey cells: total cetacean records < ", min_records_all, "."
         )) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0(year_min, "-", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    theme_map()
  p2 <- suppressMessages(add_spatial_layers(p2, land, WEA_wind, study_area, survey_area))
  suppressMessages(ggsave(file.path(out_dir, paste0(nm_file, "_2_share.png")),
         p2, width = map_width, height = map_height, dpi = map_dpi))
  
  message("  Saved 2 maps for ", nm)
  return(invisible(NULL))
}


# ---- 14. Process all target groups ----------------------------------------------

message("\n=== CREATING TARGET MAPS ===")

target_patterns <- lapply(targets, function(x) stringr::regex(x, ignore_case = TRUE))

if (use_parallel && length(targets) >= 3) {
  message("Using parallel processing (", future::availableCores() - 1, " workers)")
  future::plan(future::multisession,
               workers = min(length(targets), future::availableCores() - 1))

  # Capture PROJ_LIB in the main session so workers can inherit it.
  # Multisession workers are fresh R processes and may not find proj.db otherwise.
  .proj_lib <- Sys.getenv("PROJ_LIB")

  furrr::future_walk(
    names(targets),
    function(nm) {
      # Restore PROJ_LIB in the worker session
      if (nzchar(.proj_lib)) Sys.setenv(PROJ_LIB = .proj_lib)
      # Re-declare conflict preferences in each worker session
      conflicted::conflicts_prefer(lubridate::year,  .quiet = TRUE)
      conflicted::conflicts_prefer(lubridate::month, .quiet = TRUE)
      conflicted::conflicts_prefer(dplyr::filter,    .quiet = TRUE)
      conflicted::conflicts_prefer(dplyr::select,    .quiet = TRUE)
      create_target_maps(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
    },
    .options = furrr::furrr_options(seed = TRUE)
  )
  future::plan(future::sequential)
} else {
  for (nm in names(targets)) {
    create_target_maps(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
  }
}


# ---- 15. Seasonal summary table (WEA vs study region) ---------------------------
# Counts records per season x region for each target group and priority species.
# Priority species: NARW, Blue Whale, Harbour Porpoise.
# Guilds: Baleen Whales (Shelf_baleen), Small Odontocetes.
# Output: CSV saved to out_dir and printed to console.

message("\n=== CREATING SEASONAL SUMMARY TABLE ===")

# Pre-compile patterns for priority targets only
priority_targets <- targets[c("NARW", "Blue_whale", "Shelf_baleen",
                              "Harbour_porpoise", "Small_odontocetes")]
priority_patterns <- lapply(priority_targets, function(x) stringr::regex(x, ignore_case = TRUE))

# Tag each record with matching target group(s)
# A record can match multiple groups (e.g. NARW matches both NARW and Shelf_baleen)
df_tagged <- df %>%
  mutate(
    target_group = purrr::map_chr(name_blob, function(nb) {
      matches <- names(priority_patterns)[sapply(priority_patterns,
                                                 function(rx) stringr::str_detect(nb, rx))]
      if (length(matches) == 0) return(NA_character_)
      paste(matches, collapse = "|")
    })
  ) %>%
  filter(!is.na(target_group)) %>%
  tidyr::separate_rows(target_group, sep = "\\|")

# Seasonal summary by region
seasonal_summary <- df_tagged %>%
  mutate(region = ifelse(in_survey, "WEA", "Study Region")) %>%
  group_by(target_group, region, season) %>%
  summarise(
    n_records  = sum(w, na.rm = TRUE),
    n_years    = n_distinct(year, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  # Add both regions for all combinations (fill 0 for missing)
  tidyr::complete(target_group, region, season,
                  fill = list(n_records = 0, n_years = 0)) %>%
  mutate(
    target_label = dplyr::recode(target_group,
                                 NARW              = "North Atlantic Right Whale",
                                 Blue_whale        = "Blue Whale",
                                 Shelf_baleen      = "Baleen Whales (guild)",
                                 Harbour_porpoise  = "Harbour Porpoise",
                                 Small_odontocetes = "Small Odontocetes (guild)"),
    season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
  ) %>%
  arrange(target_label, region, season) %>%
  select(target_label, region, season, n_records, n_years)

# Wide format: seasons as columns for compact reporting
seasonal_wide <- seasonal_summary %>%
  tidyr::pivot_wider(
    names_from  = season,
    values_from = c(n_records, n_years),
    names_glue  = "{season}_{.value}"
  ) %>%
  select(target_label, region,
         Spring_n_records, Spring_n_years,
         Summer_n_records, Summer_n_years,
         Fall_n_records,   Fall_n_years,
         Winter_n_records, Winter_n_years)

readr::write_csv(seasonal_wide,
                 file.path(out_dir, "seasonal_summary_WEA_vs_study.csv"))
readr::write_csv(seasonal_wide,
                 here::here("output/data/seasonal_summary_WEA_vs_study.csv"))

message("  Saved: seasonal_summary_WEA_vs_study.csv")
message("\nSEASONAL SUMMARY (WEA vs Study Region):")
print(seasonal_wide, n = Inf, width = Inf)

# Export all-cetacean seasonal totals for Rmd denominator
# Uses the full df (all species, no guild grouping) -- avoids double-counting
# that occurs when summing across guild rows in seasonal_wide.
all_cet_seasonal_export <- df %>%
  mutate(region = ifelse(in_survey, "WEA", "Study Region")) %>%
  group_by(region, season) %>%
  summarise(
    total_all_cet = sum(w, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(season))

readr::write_csv(all_cet_seasonal_export,
                 here::here("output/data/all_cet_seasonal_totals.csv"))
message("  Saved: all_cet_seasonal_totals.csv")

# ---- 16. Species appendix table -------------------------------------------------

message("\n=== CREATING SPECIES SUMMARY TABLE ===")

yr_min_data   <- min(df$year, na.rm = TRUE)
yr_max_data   <- max(df$year, na.rm = TRUE)
n_years_total <- yr_max_data - yr_min_data + 1

if (use_datatable && data.table::is.data.table(dt)) {
  sp <- dt[, .(
    n_study_area        = sum(w, na.rm = TRUE),
    n_survey_any        = sum(w[in_survey == TRUE], na.rm = TRUE),
    n_spring_study      = sum(w[is_spring == TRUE], na.rm = TRUE),
    n_spring_survey_any = sum(w[is_spring == TRUE & in_survey == TRUE], na.rm = TRUE),
    years_study_area    = data.table::uniqueN(year, na.rm = TRUE),
    years_survey_any    = data.table::uniqueN(year[in_survey == TRUE], na.rm = TRUE),
    first_year_wea      = suppressWarnings(min(year[in_survey == TRUE], na.rm = TRUE)),
    last_year_wea       = suppressWarnings(max(year[in_survey == TRUE], na.rm = TRUE)),
    streak_wea          = longest_streak(year[in_survey == TRUE])
  ), by = .(species = get(species_field))] %>% as_tibble()
} else {
  sp <- df %>%
    group_by(species = .data[[species_field]]) %>%
    summarise(
      n_study_area        = sum(w, na.rm = TRUE),
      n_survey_any        = sum(w[in_survey], na.rm = TRUE),
      n_spring_study      = sum(w[is_spring], na.rm = TRUE),
      n_spring_survey_any = sum(w[is_spring & in_survey], na.rm = TRUE),
      years_study_area    = n_distinct(year),
      years_survey_any    = n_distinct(year[in_survey]),
      first_year_wea      = suppressWarnings(min(year[in_survey], na.rm = TRUE)),
      last_year_wea       = suppressWarnings(max(year[in_survey], na.rm = TRUE)),
      streak_wea          = longest_streak(year[in_survey]),
      .groups             = "drop"
    )
}

sp <- sp %>%
  mutate(
    first_year_wea        = ifelse(is.infinite(first_year_wea), NA, first_year_wea),
    last_year_wea         = ifelse(is.infinite(last_year_wea),  NA, last_year_wea),
    prop_years_survey     = years_survey_any / n_years_total,
    pct_in_wea            = ifelse(n_study_area > 0,
                                   round(100 * n_survey_any / n_study_area, 1), 0),
    consistency_wea       = case_when(
      years_survey_any >= 5 ~ "frequent (5+ years in WEA)",
      years_survey_any >= 2 ~ "intermittent (2-4 years in WEA)",
      years_survey_any == 1 ~ "rare (1 year in WEA)",
      TRUE                  ~ "none in WEA"
    ),
    present_in_wea = years_survey_any > 0
  ) %>%
  arrange(desc(n_survey_any), desc(n_study_area), species)

sp_renamed <- sp %>%
  rename(
    n_wea              = n_survey_any,
    n_spring_wea       = n_spring_survey_any,
    years_wea_count    = years_survey_any,
    prop_years_wea     = prop_years_survey
  )

metadata <- tibble::tribble(
  ~Column,            ~Description,                                                     ~Data_Type, ~Notes,
  "species",          "Species common name",                                            "text",     "Based on species_field config",
  "n_study_area",     "Total sightings across entire study area",                       "integer",  "All locations within study boundary",
  "n_wea",            "Total sightings within WEA polygons",                            "integer",  "Subset of n_study_area",
  "n_spring_study",   "Spring sightings in study area (March-May)",                     "integer",  "Seasonal subset",
  "n_spring_wea",     "Spring sightings within WEA polygons (March-May)",               "integer",  "Seasonal subset",
  "years_study_area", "Distinct years seen in study area",                              "integer",  paste0("Out of ", n_years_total, " years (", yr_min_data, "-", yr_max_data, ")"),
  "years_wea_count",  "Distinct years seen in WEA",                                    "integer",  "Subset of years_study_area",
  "first_year_wea",   "First year observed in WEA",                                    "integer",  "NA if never seen in WEA",
  "last_year_wea",    "Most recent year observed in WEA",                              "integer",  "NA if never seen in WEA",
  "streak_wea",       "Longest consecutive year streak in WEA",                        "integer",  "Max consecutive years with sightings",
  "prop_years_wea",   "Proportion of total years with WEA sightings",                  "decimal",  "years_wea_count / n_years_total",
  "pct_in_wea",       "Percentage of total sightings in WEA",                          "decimal",  "(n_wea / n_study_area) x 100",
  "consistency_wea",  "Temporal consistency category for WEA sightings",               "text",     "frequent / intermittent / rare / none",
  "present_in_wea",   "Presence flag for WEA",                                         "logical",  "TRUE if years_wea_count > 0"
)

analysis_info <- tibble::tribble(
  ~Parameter,            ~Value,
  "Date generated",      as.character(Sys.Date()),
  "Study area",          basename(study_area_path),
  "WEA polygons",        basename(WEA_polygons_path),
  "Survey area used",    as.character(use_survey_area),
  "Year range",          paste0(year_min, "-", year_max),
  "Grid size (km)",      as.character(grid_km),
  "Spring months",       paste(month.name[spring_months], collapse = ", "),
  "Species aggregation", "Seals, Fin/Sei Whale, Unidentified Whale",
  "Total species",       as.character(nrow(sp_renamed)),
  "Species in WEA",      as.character(sum(sp_renamed$present_in_wea)),
  "Data years",          paste0(yr_min_data, "-", yr_max_data, " (", n_years_total, " years)")
)

writexl::write_xlsx(
  list("Species_Data"    = sp_renamed,
       "Column_Metadata" = metadata,
       "Analysis_Info"   = analysis_info),
  path = appendix_xlsx
)
message("  Saved: ", appendix_xlsx)
message("  Species richness: ", nrow(sp_renamed))


# ---- 17. Individual species maps ------------------------------------------------

if (make_species_maps) {
  
  message("\n=== CREATING INDIVIDUAL SPECIES MAPS ===")
  
  create_species_map <- function(sp_name, dt_input, use_dt = FALSE) {
    sp_label <- as.character(sp_name)
    sp_file  <- nm_safe(sp_label)
    
    if (use_dt && data.table::is.data.table(dt_input)) {
      cell_sp <- dt_input[get(species_field) == sp_label, .(
        n_sp = sum(w, na.rm = TRUE)
      ), by = .(cell_id, gx, gy, xc, yc)] %>% as_tibble()
      pts_sp  <- as.data.frame(dt_input[get(species_field) == sp_label, .(x_m, y_m)])
    } else {
      cell_sp <- dt_input %>%
        filter(.data[[species_field]] == sp_label) %>%
        group_by(cell_id, gx, gy, xc, yc) %>%
        summarise(n_sp = sum(w, na.rm = TRUE), .groups = "drop")
      pts_sp  <- dt_input %>%
        filter(.data[[species_field]] == sp_label) %>%
        select(x_m, y_m)
    }
    
    pts_sp <- pts_sp %>% filter(!is.na(x_m), !is.na(y_m))
    if (nrow(cell_sp) == 0) return(invisible(NULL))
    
    category_colors <- c("1-5" = "#FDE724", "6-25" = "#5DC863",
                         "26-50" = "#21908C", "51-100" = "#3B528B", "100+" = "#440154")
    names(category_colors) <- species_map_labels
    
    cell_sp <- cell_sp %>%
      mutate(count_category = cut(n_sp, breaks = species_map_bins,
                                  labels = species_map_labels,
                                  right = TRUE, include.lowest = FALSE),
             count_category = factor(count_category, levels = species_map_labels))
    
    # Dummy rows to force full legend
    all_categories <- data.frame(
      cell_id = paste0("dummy_", species_map_labels),
      gx = NA_real_, gy = NA_real_, xc = NA_real_, yc = NA_real_,
      n_sp = NA_real_,
      count_category = factor(species_map_labels, levels = species_map_labels)
    )
    cell_sp_plot <- dplyr::bind_rows(cell_sp, all_categories)
    
    p <- ggplot(cell_sp_plot) +
      geom_tile(aes(x = xc, y = yc, fill = count_category),
                width = cell_m, height = cell_m) +
      scale_fill_manual(values = category_colors, name = "Sightings",
                        na.value = "grey90", drop = FALSE) +
      geom_point(data = pts_sp, aes(x = x_m, y = y_m),
                 inherit.aes = FALSE, shape = 21, fill = NA,
                 color = "grey40", alpha = 0.6, size = 1, stroke = 0.35) +
      labs(title   = paste0("Opportunistic sightings: ", sp_label),
           caption = "Standardised bins allow cross-species comparison. Raw counts reflect distribution and effort.") +
      annotate("text", x = -Inf, y = Inf,
               label = paste0(year_min, "-", year_max, " | ", grid_km, " km grid"),
               hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
      coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
      theme_map()
    
    # Use survey overlay for species maps only if use_survey_area is TRUE
    p <- suppressMessages(add_spatial_layers(p, land, WEA_wind, study_area, survey_area,
                            use_survey_overlay = use_survey_area))

    suppressMessages(ggsave(file.path(species_map_dir, paste0("species_", sp_file, "_raw.png")),
           p, width = map_width, height = map_height, dpi = map_dpi))
    invisible(NULL)
  }
  
  sp_filtered <- sp %>%
    filter(!is.na(species), species != "")
  
  if (exclude_zero_WEA) {
    sp_filtered <- sp_filtered %>% filter(n_survey_any > 0)
  }
  if (exclude_beaked_whales) {
    sp_filtered <- sp_filtered %>%
      filter(!grepl("BEAKED|Mesoplodon|Hyperoodon|Ziphius", species, ignore.case = TRUE))
  }
  sp_filtered <- sp_filtered %>%
    filter(!grepl("NORTHERN\\s+BOTTLENOSE|N\\s+BOTTLENOSE|BELUGA",
                  species, ignore.case = TRUE))
  
  sp_to_map <- sp_filtered %>% slice_head(n = top_n_species_maps) %>% pull(species)
  
  message("  Creating maps for ", length(sp_to_map), " species")
  
  if (use_parallel && length(sp_to_map) >= 8) {
    future::plan(future::multisession,
                 workers = min(length(sp_to_map), future::availableCores() - 1))
    furrr::future_walk(
      sp_to_map,
      function(spp) {
        conflicted::conflicts_prefer(lubridate::year,  .quiet = TRUE)
        conflicted::conflicts_prefer(lubridate::month, .quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::filter,    .quiet = TRUE)
        conflicted::conflicts_prefer(dplyr::select,    .quiet = TRUE)
        create_species_map(spp, dt, use_dt = use_datatable)
      },
      .options = furrr::furrr_options(seed = TRUE)
    )
    future::plan(future::sequential)
  } else {
    for (spp in sp_to_map) create_species_map(spp, dt, use_dt = use_datatable)
  }
  
  message("  Species maps saved to: ", species_map_dir)
  
  # Compile species maps into a single PDF
  species_pngs <- list.files(species_map_dir, pattern = "\\.png$", full.names = TRUE)
  pdf(file.path(out_dir, "SPECIES_MAPS.pdf"), width = 11, height = 8.5)
  for (i in seq_along(species_pngs)) {
    if (i > 1) grid::grid.newpage()
    img <- png::readPNG(species_pngs[i])
    grid::grid.raster(img)
  }
  dev.off()
  message("  Saved species maps PDF: SPECIES_MAPS.pdf")
}


# ---- 18. Export parameters for report -------------------------------------------

effort_params <- tibble::tribble(
  ~Parameter,        ~Value,                         ~Script,         ~Notes,
  "grid_km",         as.character(grid_km),          "04_summarize_wsdb_effort.R", "Grid cell size (km)",
  "year_min",        as.character(year_min),         "04_summarize_wsdb_effort.R", "Earliest year in effort maps",
  "year_max",        as.character(year_max),         "04_summarize_wsdb_effort.R", "Latest year in effort maps",
  "min_records_all", as.character(min_records_all),  "04_summarize_wsdb_effort.R", "Min records for effort-normalised maps",
  "use_survey_area", as.character(use_survey_area),  "04_summarize_wsdb_effort.R", "Survey area overlay active"
)
readr::write_csv(effort_params, here::here("output/data/params_effort.csv"))

# ---- 18. Source confidence summary ----------------------------------------------
# Run 05_summarize_wsdb_confidence.R immediately after effort maps so all downstream
# CSVs (overall_conf_summary.csv, target_conf_summary.csv) reflect the current
# df and cell_all objects. This ensures the Rmd always reads fresh data and
# never silently uses stale counts from a previous run.
#
# Required objects passed through from this session (no re-loading needed):
#   cell_all, dt, targets, target_patterns, target_titles,
#   WEA_wind, out_dir, cell_m, use_datatable

# NOTE: 05_summarize_wsdb_confidence.R is run separately via RUN_PIPELINE.R.
# The objects it needs (cell_all, dt, targets, target_patterns, target_titles,
# WEA_wind, out_dir, cell_m, use_datatable) remain available in the session
# as long as 04 and 05 are run in the same R session via the controller.
# ---- 19. Completion summary -----------------------------------------------------

message("\n", paste(rep("=", 70), collapse = ""))
message("04_summarize_wsdb_effort.R COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message("Output directory: ", out_dir)
message("Maps:   coverage, confidence, ", length(targets), " groups x 2 = ",
        length(targets) * 2, " target maps")
if (make_species_maps) message("        ", length(sp_to_map), " individual species maps")
message("Tables: APPENDIX_species_summary.xlsx")
message("        seasonal_summary_WEA_vs_study.csv")
message("        params_effort.csv")
message("Objects available for 05_summarize_wsdb_confidence.R:")
message("  cell_all, dt, targets, target_patterns, target_titles,")
message("  WEA_wind, out_dir, cell_m, use_datatable")
message(paste(rep("=", 70), collapse = ""))
