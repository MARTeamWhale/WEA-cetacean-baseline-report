# HEADER --------------------------------------------
#
# Author: Laura Joan Feyrer
# Email:  ljfeyrer@dal.ca
# Most Recent Date Updated: 2026-03-10
#
# Script Name: 07_summarize_pam_detections.R
#
# Description:
## Summarises baleen whale PAM detection data by season and species.
## Produces plots of effort, detection days, geographic range, and site maps.
## Exports summary CSVs for use in RMarkdown reports.
#
# Changes from previous version:
# - Added PAM_BUFFER_DIST_M (10 km) to match KDE inclusion buffer
# - Station clipping now uses buffered study area (was: exact boundary)
# - Added WHALE_PALETTE: named colour palette by common name (reusable across scripts)
# - Species labels on plots now show common names (not Latin initials)
# - Combined layout trimmed to top two panels (p1, p2) only; legend in one row
# - Added species x season detection days summary table (Section 10)
# - Exports: pam_species_season_det_days.csv
#

# SET OPTIONS -------------------------------------------------------
options(scipen = 999)

# LIBRARIES --------------------------------------------------------
pacman::p_load(terra, spatstat, sf, dplyr, patchwork, ggplot2, lubridate, grafify, tidyr, readr)
suppressWarnings(source(here::here("scripts/00_load_helpers.R")))

filepath = "input/raw_data/PAM/baleen_presence_laura_2025.csv"
dir.create("output/figs/Effort_maps/PAM", showWarnings = FALSE, recursive = TRUE)
dir.create("output/data", showWarnings = FALSE, recursive = TRUE)

# PROJECTIONS ------------------------------------------------------
UTM20 <- SPATIAL_CRS_UTM20 # UTM Zone 20N

# PAM BUFFER DISTANCE ----------------------------------------------
# Match the buffer used in run_all_kdes: stations within this distance
# of the study area boundary are included (was: exact intersection only)
PAM_BUFFER_DIST_M <- 10000  # 10 km — must match KDE script

# WHALE SPECIES PALETTE --------------------------------------------
# Inspired by Wes Anderson "Life Aquatic" + "Grand Budapest Hotel" —
# desaturated, dusty, film-stock feel. No light yellows.
#
# Design logic:
#   Baleen whales  → blue/teal/slate family + terracotta accent for Humpback
#                    burgundy for Right Whale (rare/endangered visual weight)
#   Beaked whales  → warm sandy-earth family, clearly distinct from baleen group
#
# Copy SPECIES_COMMON_NAMES + WHALE_PALETTE into any script for consistent colours.
# Use with:  scale_fill_manual(values = WHALE_PALETTE)  (fill = common_name)
#            scale_colour_manual(values = WHALE_PALETTE) (colour = common_name)

# Latin code → common name lookup (baleen + beaked species)
SPECIES_COMMON_NAMES <- c(
  "Bb"   = "Sei Whale",
  "Bm"   = "Blue Whale",
  "Bp"   = "Fin Whale",
  "Mn"   = "Humpback Whale",
  "Ba"   = "Minke Whale",
  "Eg"   = "North Atlantic Right Whale",
  "Ha"   = "Northern Bottlenose Whale",
  "Mb"   = "Sowerby's Beaked Whale",
  "Zc"   = "Cuvier's Beaked Whale",
  "MmMe" = "True's/Gervais' Beaked Whale"
)

# Palette keyed by common name — keep in sync with 03_plot_combined_sightings.R
WHALE_PALETTE <- c(
  # ── Baleen: blue-teal family ──────────────────────────────────────────────
  "Blue Whale"                   = "#3288BD",
  "Fin Whale"                    = "#276B95",
  "Sei Whale"                    = "#2A5857",
  "Fin/Sei Whale"                = "#7AB0C0",
  "Humpback Whale"               = "#6DAFB1",
  "Minke Whale"                  = "#9bc4f8",
  "North Atlantic Right Whale"   = "#7C6FB3",
  # ── Beaked: warm sandy-earth family ──────────────────────────────────────
  "Northern Bottlenose Whale"    = "#B07D62",
  "Sowerby's Beaked Whale"       = "#8B6B4E",
  "Cuvier's Beaked Whale"        = "#C9A97A",
  "True's/Gervais' Beaked Whale" = "#7A6651",
  # ── Dolphins & other odontocetes ─────────────────────────────────────────
  "Common Dolphin"               = "#F1B2A1",
  "Atlantic Bottlenose Dolphin"  = "#D53E4F",
  "Atlantic White-Sided Dolphin" = "#E57A7D",
  "White-Beaked Dolphin"         = "#ABDDA4",
  "Striped Dolphin"              = "#66C2A5",
  "Risso's Dolphin"              = "#A8627A",
  "Long-Finned Pilot Whale"      = "#f7c6d5",
  "Sperm Whale"                  = "#6B9527",
  "Harbour Porpoise"             = "#FDAE61"
)

# Read baleen whale PAM season data------
baleen_DOY <- read.csv(filepath)
baleen_DOY %>% group_by(species) %>% summarise(count = n())

whale_data <- baleen_DOY %>%
  mutate(UTC = as.POSIXct(rec_date, format = "%Y-%m-%d")) %>%
  mutate(Date = format(as_date(rec_date), "%Y-%m-%d")) %>%
  mutate(
    month = lubridate::month(UTC),
    Season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(site, Season) %>%
  mutate(
    season_days = n_distinct(Date)
  ) %>%
  ungroup()

# ── CLIP whale_data TO BUFFERED STUDY AREA ───────────────────────────────────
# All summaries, plots, and exports downstream use whale_data, so filtering
# here ensures only stations within PAM_BUFFER_DIST_M of the study area are
# included — consistent with the KDE inclusion logic in run_all_kdes.
# Done here (not later) so baleen_PA_season and all section 1-8 summaries
# are already scoped to the study area.
# Load study_area if not already in environment
if (!exists("study_area") || is.null(study_area)) {
  study_area <- load_spatial_layer("study_area", crs = UTM20)
  cat("Loaded study_area from", spatial_path("study_area"), "\n")
}

study_area_pam_buffer <- st_buffer(
  st_transform(study_area, UTM20),
  dist = PAM_BUFFER_DIST_M
)
whale_data_sf <- whale_data %>%
  distinct(site, latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

sites_in_buffer <- whale_data_sf %>%
  st_filter(study_area_pam_buffer) %>%
  pull(site)

n_before <- n_distinct(whale_data$site)
whale_data <- whale_data %>% filter(site %in% sites_in_buffer)
n_after  <- n_distinct(whale_data$site)

cat(sprintf(
  "\nStudy area filter: %d stations total → %d retained within %g km buffer\n",
  n_before, n_after, PAM_BUFFER_DIST_M / 1000
))

baleen_PA_season <- whale_data %>%
  group_by(site, Season, species) %>%
  summarise(
    efort_days     = dplyr::first(season_days),
    detection_days = sum(presence > 0),
    proportion_det = detection_days / efort_days,
    latitude       = dplyr::first(latitude),
    longitude      = dplyr::first(longitude),
    .groups = "drop"
  )

# check data
season <- whale_data %>% group_by(Season, site) %>% summarise(days = n())

# Date range
whale_data %>% summarise(min(Date), max(Date))

# Seasonal Site Distribution and Species Analysis Plots
# Required: baleen_PA_season and whale_data from your existing code

# 1. SUMMARY STATISTICS ----

# Number of sites per year
sites_per_year <- whale_data %>%mutate(year = lubridate::year(UTC))%>%
  group_by(year) %>%
  summarise(
    n_sites = n_distinct(site),
    lat_min = min(latitude),
    lat_max = max(latitude),
    lon_min = min(longitude),
    lon_max = max(longitude),
    lat_range = lat_max - lat_min,
    lon_range = lon_max - lon_min
  )

# Number of sites per season
sites_per_season <- whale_data %>%
  group_by(Season) %>%
  summarise(
    n_sites = n_distinct(site),
    lat_min = min(latitude),
    lat_max = max(latitude),
    lon_min = min(longitude),
    lon_max = max(longitude),
    lat_range = lat_max - lat_min,
    lon_range = lon_max - lon_min
  )

# Number of years per season
years_per_season <- whale_data %>%
  mutate(year = lubridate::year(UTC)) %>%
  group_by(Season) %>%
  summarise(
    n_years = n_distinct(year),
    year_min = min(year),
    year_max = max(year)
  )

# Combine seasonal summaries
seasonal_summary <- left_join(sites_per_season, years_per_season, by = "Season")

print("Seasonal Summary:")
print(seasonal_summary)

# 2. PLOT: Number of Sites Sampled per Year by Season ----
sites_year_season <- whale_data %>%
  mutate(year = lubridate::year(UTC)) %>%
  group_by(Season, year) %>%
  summarise(n_sites = n_distinct(site), .groups = "drop")

p1 <- ggplot(
  sites_year_season,
  aes(
    x = factor(year), y = n_sites,
    fill = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))
  )
) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "Winter" = "#2E86AB", "Spring" = "#06A77D",
      "Summer" = "#F18F01", "Fall" = "#A23B72"
    )
  ) +
  labs(
    title = "",
    x = "",
    y = "Number of Sites",
    fill = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 3. PLOT: Geographic Range per Season ----
range_data <- seasonal_summary %>%
  tidyr::pivot_longer(
    cols = c(lat_range, lon_range),
    names_to = "degrees",
    values_to = "range"
  ) %>%
  mutate(degrees = ifelse(degrees == "lat_range", "Latitude", "Longitude"))

p2 <- ggplot(range_data, aes(
  x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
  y = range, fill = degrees
)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Latitude" = "#06A77D", "Longitude" = "#D4A373")) +
  labs(
    title = "Distribution of Sites by Season",
    x = "",
    y = "Range *",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# 4. SEASONAL SUMMARY BY SPECIES ----

species_season_summary <- whale_data %>%
  mutate(year = lubridate::year(UTC)) %>%
  group_by(Season, species) %>%
  summarise(
    n_sites        = n_distinct(site),
    n_years        = n_distinct(year),
    detection_days = sum(presence > 0),  # days with actual detections
    effort_days    = n(),                 # total days sampled
    .groups = "drop"
  )

print("\nSpecies by Season Summary:")
print(species_season_summary)

# ── MAP LATIN CODES TO COMMON NAMES for plotting ──────────────────────────────
# Adds a common_name column to any data frame that has a 'species' code column.
# Uses WHALE_PALETTE (keyed by common name) for consistent colours across plots.

species_season_summary <- species_season_summary %>%
  mutate(common_name = recode(species, !!!SPECIES_COMMON_NAMES))

baleen_PA_season_named <- baleen_PA_season %>%
  mutate(common_name = recode(species, !!!SPECIES_COMMON_NAMES))

# Subset palette to only species present in this dataset
present_names  <- unique(species_season_summary$common_name)
palette_subset <- WHALE_PALETTE[names(WHALE_PALETTE) %in% present_names]

# 5. PLOT: Sites with Detections per Season by Species ----
# Check variability in site counts by species
site_species_check <- species_season_summary %>%
  select(Season, common_name, n_sites) %>%
  tidyr::pivot_wider(names_from = common_name, values_from = n_sites)

print("\nSites per season by species (check for variability):")
print(site_species_check)

# Proportion of sites with detections for each species
species_detection_prop <- baleen_PA_season_named %>%
  group_by(Season, common_name) %>%
  summarise(
    sites_with_detections = sum(detection_days > 0),
    total_sites           = n(),
    prop_sites_detected   = sites_with_detections / total_sites,
    .groups = "drop"
  )

p3 <- ggplot(
  species_detection_prop,
  aes(
    x    = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
    y    = prop_sites_detected,
    fill = common_name
  )
) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = palette_subset, name = "") +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(title = "% Sites with Detections", x = "", y = "% Sites") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical",
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 2))

# 6. PLOT: Detection Days by Season and Species ----
p4 <- ggplot(
  species_season_summary,
  aes(
    x    = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
    y    = detection_days,
    fill = common_name
  )
) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = palette_subset, name = "") +
  labs(title = "Detection Days", x = "", y = "Detection Days") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.direction = "vertical",
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 2))

# # 7. PLOT: Years Sampled by Season and Species (kept for reference, not in layout) ----
# p5 <- ggplot(
#   species_season_summary,
#   aes(
#     x    = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
#     y    = n_years,
#     fill = common_name
#   )
# ) +
#   geom_col(position = "dodge", alpha = 0.85) +
#   scale_fill_manual(values = palette_subset, name = "") +
#   labs(title = "Sampling Years", x = "", y = "Number of Years") +
#   theme_minimal(base_size = 12) +
#   theme(legend.position = "none", panel.grid.minor = element_blank())

# 8. COMBINED LAYOUT — top two panels (p3 % sites detected | p4 detection days) ----
# Correct patchwork pattern: legend.position for the collected guide goes in
# plot_annotation(theme = ...), NOT via & theme() or on individual plots.

p3_bottom <- p3 + theme(legend.position = "bottom", legend.direction = "vertical", )
p4_bottom <- p4 + theme(legend.position = "bottom", legend.direction = "vertical", )

combined_seasonal <- p3_bottom / p4_bottom +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Seasonal Summary of Baleen Whale Detections",
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  legend.position = "bottom", legend.direction = "vertical")
  )

print(combined_seasonal)
ggsave("output/figs/Effort_maps/PAM/combined_seasonal.png", combined_seasonal, h = 12, w = 6, units = "in")

# 9. Site Distribution Map by Season & Year ----
# FIXED VERSION WITH PROPER STUDY AREA AND OSW OVERLAY

# Load OSW wind areas if not already loaded
if (!exists("osw_wind")) {
  osw_wind <- load_spatial_layer("wea", crs = UTM20, required = FALSE)
}

# Load contours if not already loaded
if (!exists("cont")) {
  cont <- load_bathy_contours(crs = UTM20, required = FALSE)
}

# Load land if not already loaded
if (!exists("land")) {
  land <- load_spatial_layer("land", crs = UTM20, required = FALSE)
} else if (exists("land") && !is.null(land)) {
  # Ensure land is in UTM20
  if (st_crs(land) != st_crs(UTM20)) {
    land <- st_transform(land, UTM20)
  }
}

# Prepare site locations and transform to UTM20
site_locations <- whale_data %>%
  mutate(year = lubridate::year(UTC)) %>%
  group_by(site, Season) %>%
  summarise(
    latitude = dplyr::first(latitude),
    longitude = dplyr::first(longitude),
    n_species = n_distinct(species),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

# Get bounding box from study area if available, otherwise from data
if (exists("study_area") && !is.null(study_area)) {
  bbox <- st_bbox(study_area)
  xlims <- c(bbox["xmin"], bbox["xmax"])
  ylims <- c(bbox["ymin"], bbox["ymax"])
} else {
  # Fallback to data extent with buffer
  bbox <- st_bbox(site_locations)
  buffer <- 50000 # 50km buffer in meters (UTM)
  xlims <- c(bbox["xmin"] - buffer, bbox["xmax"] + buffer)
  ylims <- c(bbox["ymin"] - buffer, bbox["ymax"] + buffer)
}

p_map <- ggplot() +
  # Add land first
  {
    if (exists("land") && !is.null(land)) {
      geom_sf(data = land, fill = "grey60", color = NA)
    }
  } +
  # Add points
  geom_sf(
    data = site_locations,
    aes(color = Season, alpha = n_years),
    size = 3
  ) +
  # Add OSW areas
  {
    if (exists("osw_wind") && !is.null(osw_wind)) {
      geom_sf(
        data = osw_wind, fill = NA, color = "#FF6B35",
        linewidth = 1, inherit.aes = FALSE
      )
    }
  } +
  # Add study area outline
  {
    if (exists("study_area") && !is.null(study_area)) {
      geom_sf(
        data = study_area, fill = NA, color = "black",
        linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
      )
    }
  } +
  scale_alpha_continuous(range = c(0.3, 1), name = "Years Sampled") +
  scale_color_manual(
    values = c(
      "Winter" = "#2E86AB", "Spring" = "#06A77D",
      "Summer" = "#F18F01", "Fall" = "#A23B72"
    )
  ) +
  coord_sf(xlim = xlims, ylim = ylims, crs = st_crs(UTM20), expand = FALSE) +
  facet_wrap(~Season) +
  labs(
    title = "Site Level Annual Sample Effort",
    x = "",
    y = ""
  ) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_map)
ggsave("output/figs/Effort_maps/PAM/Site_effort_map.png", p_map, width = 12, height = 10, dpi = 300)

# ============================
# PAM COVERAGE: CONCENTRIC RANGE BUFFERS (single map)--------
# ============================

# 1) Build a station layer
pam_tbl <- whale_data

has_site <- any(names(pam_tbl) %in% c("site", "station", "Station", "Site"))
site_col <- base::intersect(names(pam_tbl), c("site", "station", "Station", "Site"))[1]

stations_sf <- if (has_site) {
  pam_tbl %>%
    mutate(site_id = .data[[site_col]]) %>%
    group_by(site_id) %>%
    summarise(
      longitude = dplyr::first(longitude),
      latitude = dplyr::first(latitude),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM20)
} else {
  pam_tbl %>%
    mutate(
      site_id = paste0(round(longitude, 4), "_", round(latitude, 4))
    ) %>%
    group_by(site_id) %>%
    summarise(
      longitude = dplyr::first(longitude),
      latitude = dplyr::first(latitude),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM20)
}

# Store total stations before clipping
total_stations <- nrow(stations_sf)

# Clip stations_sf to buffered study area for coverage map.
# whale_data is already filtered to these stations (clipped at load time above),
# so counts will match. study_area must exist — loaded in load_basemap_shapes.R
# or the PAM data load block. If this errors, check study_area is in your env.
if (!exists("study_area") || is.null(study_area)) {
  stop("study_area not found — cannot clip stations_sf for p_pam_range. ",
       "Load it before running this section.")
}

study_area_pam_buffer <- st_buffer(
  st_transform(study_area, UTM20), dist = PAM_BUFFER_DIST_M
)
stations_sf <- st_intersection(stations_sf, study_area_pam_buffer)

cat(sprintf("stations_sf clipped: %d total → %d within %g km buffer\n",
            total_stations, nrow(stations_sf), PAM_BUFFER_DIST_M / 1000))

study_area_stations <- nrow(stations_sf)

# SUMMARY STATISTICS FOR STUDY AREA STATIONS ----
cat("\n=== PAM STATION SUMMARY ===\n")
cat("Total PAM stations (all data):", total_stations, "\n")
cat("Stations within study area:", study_area_stations, "\n")
cat("Stations excluded:", total_stations - study_area_stations, "\n")
cat("Percentage in study area:", round(100 * study_area_stations / total_stations, 1), "%\n\n")

# Get date range and species from whale_data
date_summary <- whale_data %>%
  summarise(
    min_date = min(UTC, na.rm = TRUE),
    max_date = max(UTC, na.rm = TRUE),
    total_days = n_distinct(Date),
    n_species = n_distinct(species)
  )

species_list <- whale_data %>%
  distinct(species) %>%
  pull(species) %>%
  sort()

cat("Date range:", format(date_summary$min_date, "%Y-%m-%d"), "to", 
    format(date_summary$max_date, "%Y-%m-%d"), "\n")
cat("Total unique detection days:", date_summary$total_days, "\n")
cat("Species detected:", paste(species_list, collapse = ", "), "\n")

# Seasonal coverage summary
seasonal_station_summary <- whale_data %>%
  filter(site %in% stations_sf$site_id) %>%
  group_by(Season) %>%
  summarise(
    n_stations = n_distinct(site),
    n_detection_days = n_distinct(Date),
    n_years = n_distinct(lubridate::year(UTC)),
    .groups = "drop"
  ) %>%
  arrange(match(Season, c("Winter", "Spring", "Summer", "Fall")))

cat("\n=== SEASONAL COVERAGE (Study Area Stations Only) ===\n")
print(seasonal_station_summary)

# Species-specific summary
species_station_summary <- whale_data %>%
  filter(site %in% stations_sf$site_id) %>%
  group_by(species) %>%
  summarise(
    n_stations = n_distinct(site),
    n_detection_days = sum(presence > 0),
    pct_stations = round(100 * n_distinct(site) / study_area_stations, 1),
    .groups = "drop"
  )

cat("\n=== SPECIES DETECTION SUMMARY (Study Area Stations Only) ===\n")
print(species_station_summary)
cat("\n")

# 2) Define range scenarios (km)-----
ranges_km <- c(50, 25, 10, 5) # plot big-to-small so circles look concentric
ranges_m <- ranges_km * 1000

# 3) Create buffers (after clipping stations)
buffers_sf <- purrr::map2_dfr(
  ranges_km, ranges_m,
  ~ {
    st_buffer(stations_sf, dist = .y) %>%
      mutate(range_km = .x)
  }
)

# 4) Style: blue gradient by range (darker near station)
# Need to match the order we'll use in the legend
range_fill <- c(
  "50" = scales::alpha("#6BAED6", 0.15),
  "25" = scales::alpha("#2171B5", 0.20),
  "10" = scales::alpha("#08519C", 0.30),
  "5" = scales::alpha("#08306B", 0.35)
)

buffers_sf$range_km <- factor(buffers_sf$range_km, levels = c(50, 25, 10, 5))

# 5) Get proper extent for the map -----
# Use study area bounds if available, otherwise use buffer extent
if (exists("study_area") && !is.null(study_area)) {
  map_bbox <- st_bbox(study_area)
  map_xlims <- c(map_bbox["xmin"], map_bbox["xmax"])
  map_ylims <- c(map_bbox["ymin"], map_bbox["ymax"])
  
  # Crop land to bbox (not study area polygon since study area is marine)
  if (exists("land") && !is.null(land)) {
    # Add buffer to bbox to ensure land edges are included
    bbox_buffered <- st_as_sfc(map_bbox) %>% 
      st_buffer(10000) %>%  # 10km buffer
      st_bbox()
    land_cropped <- suppressWarnings(st_crop(land, bbox_buffered))
  }
} else {
  # Use the outermost buffer extent
  map_bbox <- st_bbox(buffers_sf)
  map_xlims <- c(map_bbox["xmin"], map_bbox["xmax"])
  map_ylims <- c(map_bbox["ymin"], map_bbox["ymax"])
  
  if (exists("land") && !is.null(land)) {
    bbox_buffered <- st_as_sfc(map_bbox) %>% 
      st_buffer(10000) %>%
      st_bbox()
    land_cropped <- suppressWarnings(st_crop(land, bbox_buffered))
  }
}

# 6) Plot---
p_pam_range <- ggplot() +
  # Land first
  {
    if (exists("land_cropped") && !is.null(land_cropped)) {
      geom_sf(data = land_cropped, fill = "grey75", color = "grey50", 
              linewidth = 0.3, inherit.aes = FALSE)
    }
  } +
  # # Range buffers
  # geom_sf(
  #   data = buffers_sf,
  #   aes(fill = as.character(range_km)),
  #   color = NA,
  #   inherit.aes = FALSE
  # ) +
  # Stations on top
  geom_sf(data = stations_sf, shape = 21, fill = "darkblue", color = "darkblue", alpha = .5,
          size = 2, stroke = 0.8, inherit.aes = FALSE) +
  # OWA outline
  {
    if (exists("osw_wind") && !is.null(osw_wind)) {
      geom_sf(data = osw_wind, fill = NA, color = "#FF6B35", 
              linewidth = 1, inherit.aes = FALSE)
    }
  } +
# contours  
    { if (exists("cont") && !is.null(cont)) {
      geom_sf(data = cont, fill = NA, color = "grey", 
              linewidth = 0.2, inherit.aes = FALSE)
    }}+
  # Study area outline
  {
    if (exists("study_area") && !is.null(study_area)) {
      geom_sf(data = study_area, fill = NA, color = "black", 
              linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE)
    }
  } +
  # scale_fill_manual(
  #   values = range_fill,
  #   name = "Range from \nRecorder (Km)",
  #   breaks = c("10", "25", "50"),
  #   labels = c("10", "25", "50"))
  
  guides(fill = guide_legend(
    override.aes = list(alpha = 0.5), 
    nrow = 1,
    label.position = "bottom",
    title.position = "top"
  )) +
  coord_sf(xlim = map_xlims, ylim = map_ylims, crs = st_crs(UTM20), expand = T) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
  labs(
    title    = "PAM Station Coverage",
    subtitle = paste0("n = ", nrow(stations_sf), " stations within ",
                      PAM_BUFFER_DIST_M / 1000, " km of study area")
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    # axis.text = element_blank(),
    # axis.ticks = element_blank(),
    legend.position = c(0.98, 0.058),
    legend.key.spacing.x = unit(.02, "cm"),
    legend.key.justification = "center",
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.box.background = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(0.58, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b = 10))
  )

ggsave(
  filename = "output/figs/Effort_maps/PAM/PAM_station_coverage.png",
  plot = p_pam_range,
  width = 10, height = 8, dpi = 300
)

print(p_pam_range)

# pam_summary for rmd:
pam_summary <- whale_data %>%
  group_by(species) %>%
  summarise(
    n_sites         = n_distinct(site),
    total_days      = n(),
    detection_days  = sum(presence > 0, na.rm = TRUE),
    detection_rate  = round(detection_days / total_days, 3)
  )
write_csv(pam_summary, "output/data/pam_summary.csv")

# Seasonal station summary for rmd
write_csv(seasonal_station_summary, "output/data/pam_seasonal_summary.csv")

# ── 10. SPECIES x SEASON DETECTION DAYS SUMMARY ────────────────────────────

# Detection days per species per season, using buffered study area stations
# Filters to stations included in analysis (matching KDE buffer logic)

species_season_det <- whale_data %>%
  # Keep only stations that pass the buffer filter
  filter(site %in% stations_sf$site_id) %>%
  group_by(species, Season) %>%
  summarise(
    detection_days  = sum(presence > 0, na.rm = TRUE),
    effort_days     = n_distinct(Date),
    detection_rate  = round(detection_days / effort_days, 3),
    n_sites         = n_distinct(site),
    .groups         = "drop"
  ) %>%
  # Enforce season order for display
  mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))) %>%
  arrange(species, Season)

cat("\n=== SPECIES x SEASON DETECTION DAYS (buffered study area stations) ===\n")
print(species_season_det, n = Inf)

# Wide version — detection days as a species x season matrix (handy for tables)
species_season_wide <- species_season_det %>%
  select(species, Season, detection_days) %>%
  tidyr::pivot_wider(
    names_from  = Season,
    values_from = detection_days,
    values_fill = 0
  )

cat("\n=== DETECTION DAYS MATRIX (species x season) ===\n")
print(species_season_wide)

# Save both formats
write_csv(species_season_det,  "output/data/pam_species_season_det_days.csv")
write_csv(species_season_wide, "output/data/pam_species_season_det_days_wide.csv")

cat("\nSaved: output/data/pam_species_season_det_days.csv (long)\n")
cat("Saved: output/data/pam_species_season_det_days_wide.csv (wide matrix)\n")

# ── 11. PLOT: Detection Days Heatmap — Species x Season ─────────────────────

p_heatmap <- ggplot(
  species_season_det,
  aes(
    x    = Season,
    y    = species,
    fill = detection_days
  )
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = detection_days), size = 3.5, color = "white", fontface = "bold") +
  scale_fill_viridis_c(option = "mako", begin = 0.1, end = 0.9,
                       name = "Detection\nDays") +
  labs(
    title    = "Baleen Whale PAM — Detection Days by Species and Season",
    subtitle = paste0("Stations within ", PAM_BUFFER_DIST_M / 1000,
                      " km of study area boundary"),
    x = "",
    y = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid    = element_blank(),
    axis.text     = element_text(size = 11),
    plot.subtitle = element_text(size = 10, color = "grey40")
  )

print(p_heatmap)
ggsave(
  "output/figs/Effort_maps/PAM/pam_species_season_heatmap.png",
  p_heatmap, width = 8, height = 5, dpi = 300
)
