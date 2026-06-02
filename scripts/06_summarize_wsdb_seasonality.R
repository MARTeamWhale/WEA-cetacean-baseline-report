# Seasonal Analysis of Cetacean Sightings--------
# Creates faceted maps showing all cetacean records by season
# Run after 04_summarize_wsdb_effort.R has completed successfully

# ---------------------------
# SETUP
# ---------------------------
message("\n=== Seasonal Coverage Analysis ===\n")

# Check that main script has been run
if (!exists("dt")) {
  stop("dt object not found. Please run 04_summarize_wsdb_effort.R first or source it.")
}

if (!exists("cell_m")) {
  stop("cell_m object not found. Please run 04_summarize_wsdb_effort.R first or source it.")
}

# Load required packages
required_pkgs <- c("dplyr", "ggplot2", "sf", "lubridate")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---------------------------
# HELPERS
# ---------------------------
# Define seasons (Northern Hemisphere)
# Winter: Dec, Jan, Feb
# Spring: Mar, Apr, May
# Summer: Jun, Jul, Aug
# Fall: Sep, Oct, Nov

add_season <- function(data, date_col = date_utc) {
  data %>%
    dplyr::mutate(
      month = lubridate::month({{ date_col }}),
      season = dplyr::case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        month %in% c(9, 10, 11) ~ "Fall",
        TRUE ~ NA_character_
      ),
      season = factor(season, levels = SEASON_LEVELS)
    )
}

# ---------------------------
# ADD SEASONAL CLASSIFICATION
# ---------------------------

message("Classifying records by season...")

if (use_datatable && "data.table" %in% class(dt)) {
  dt <- as_tibble(dt)
}

dt <- add_season(dt, date_utc)

seasonal_summary <- dt %>%
  dplyr::group_by(season) %>%
  dplyr::summarise(n_records = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    pct = round(n_records / sum(n_records, na.rm = TRUE) * 100, 1)
  )
message("\nSeasonal distribution of records:")
print(seasonal_summary)

# Calculate percentages
if (use_datatable && "data.table" %in% class(seasonal_summary)) {
  total_records <- sum(seasonal_summary$n_records, na.rm = TRUE)
  seasonal_summary[, pct := round(n_records / total_records * 100, 1)]
} else {
  total_records <- sum(seasonal_summary$n_records, na.rm = TRUE)
  seasonal_summary <- seasonal_summary %>%
    mutate(pct = round(n_records / total_records * 100, 1))
}

message("\nPercentage by season:")
print(seasonal_summary)

# ---------------------------
# CALCULATE SEASONAL GRID SUMMARIES
# ---------------------------
message("\nCalculating seasonal grid summaries...")

if (use_datatable && "data.table" %in% class(dt)) {
  
  cell_seasonal <- dt[, .(
    n_all = sum(w, na.rm = TRUE)
  ), by = .(cell_id, gx, gy, x0, y0, season)]
  
  cell_seasonal <- tibble::as_tibble(cell_seasonal)
  
} else {
  
  cell_seasonal <- dt %>%
    dplyr::group_by(cell_id, gx, gy, x0, y0, season) %>%
    dplyr::summarise(
      n_all = sum(w, na.rm = TRUE),
      .groups = "drop"
    )
}

# Remove NA seasons and add centre coordinates for grid cells
cell_seasonal <- cell_seasonal %>%
  dplyr::filter(!is.na(season)) %>%
  dplyr::mutate(
    xc = x0 + cell_m / 2,
    yc = y0 + cell_m / 2
  )

message(sprintf(
  "  Created seasonal summaries for %d unique cell-season combinations",
  nrow(cell_seasonal)
))

# ---------------------------
# CREATE SEASONAL FACETED MAP
# ---------------------------
message("\nCreating seasonal faceted map...")

# Apply log transform to match main coverage map
cell_seasonal <- cell_seasonal %>%
  mutate(n_all_display = log10(n_all + 1))

# Create breaks at powers of 10 for intuitive reading
log_breaks <- seq(0, ceiling(max(cell_seasonal$n_all_display, na.rm = TRUE)), by = 1)
actual_values <- 10^log_breaks
actual_values[1] <- 0  # First break is 0, not 1

# Create the plot
p_seasonal <- ggplot(cell_seasonal) +
  geom_tile(aes(x = xc, y = yc, fill = n_all_display), width = cell_m, height = cell_m)+
  scale_fill_viridis_c(
    option = "viridis", 
    na.value = "grey90", 
    name = "Count",
    breaks = log_breaks,
    labels = actual_values,
    begin = 0.15, 
    end = 0.95
  ) +
  facet_wrap(~season, ncol = 2) +
  labs(
    title = "Seasonal Distribution of All Cetacean Sightings",
    subtitle = paste0(year_min, "–", year_max, " | ", grid_km, " km grid | log scale"),
    caption = paste0("Seasons: Spring (Mar-May), Summer (Jun-Aug), Fall (Sep-Nov), Winter (Dec-Feb)\n",
                     "Record counts: Spring=", seasonal_summary$n_records[seasonal_summary$season=="Spring"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Spring"], "%), ",
                     "Summer=", seasonal_summary$n_records[seasonal_summary$season=="Summer"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Summer"], "%), ",
                     "Fall=", seasonal_summary$n_records[seasonal_summary$season=="Fall"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Fall"], "%), ",
                     "Winter=", seasonal_summary$n_records[seasonal_summary$season=="Winter"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Winter"], "%)")
  ) +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

# Add spatial overlays to each facet
if (!is.null(land)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = land, fill = MAP_LAYER_STYLE$land_fill,
            color = MAP_LAYER_STYLE$land_color, inherit.aes = FALSE)
}

if (!is.null(WEA_wind)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = WEA_wind, fill = NA, color = WEA_color, 
            alpha = WEA_alpha, linewidth = WEA_linewidth, inherit.aes = FALSE)
}

if (!is.null(study_area)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = study_area, fill = NA, color = MAP_LAYER_STYLE$study_area_color, 
            linewidth = 0.6, linetype = "dashed", inherit.aes = FALSE)
}

# Save the plot
ggsave(
  file.path(out_dir, "00_seasonal_coverage.png"),
  p_seasonal,
  width = 12,
  height = 10,
  dpi = map_dpi
)

message("  ✓ Saved seasonal coverage map")

# ---------------------------
# CALCULATE SEASONAL STATISTICS FOR OSW AREAS
# ---------------------------
if (!is.null(WEA_wind)) {
  message("\nCalculating seasonal statistics for OSW areas...")
  
  # Create spatial object from seasonal grid
  cells_seasonal_sf <- make_cell_polygons(
    cells = cell_seasonal,
    cell_m = cell_m,
    x_col = "xc",
    y_col = "yc",
    crs = 32620
  ) 
  
  # Find cells that intersect OSW polygons
  osw_cells_seasonal <- sf::st_intersection(cells_seasonal_sf, WEA_wind)
  
  if (nrow(osw_cells_seasonal) > 0) {
    # Summary by season
    osw_seasonal_summary <- osw_cells_seasonal %>%
      sf::st_drop_geometry() %>%
      group_by(season) %>%
      summarise(
        n_cells_with_data = sum(n_all > 0),
        total_records = sum(n_all, na.rm = TRUE),
        mean_records = round(mean(n_all[n_all > 0], na.rm = TRUE), 1),
        median_records = round(median(n_all[n_all > 0], na.rm = TRUE), 1),
        .groups = "drop"
      )
    
    message("\nOSW Areas - Seasonal Summary:")
    print(osw_seasonal_summary, width = Inf)
    
    # Save to CSV
    readr::write_csv(
      osw_seasonal_summary,
      file.path(out_dir, "confidence_summaries", "osw_seasonal_summary.csv")
    )
    message("  ✓ Saved osw_seasonal_summary.csv")
    
  } else {
    message("  Warning: No cells intersect OSW polygons")
  }
}

# ---------------------------
# SEASONAL COMPARISON BAR CHART-------
# ---------------------------

message("\nCreating seasonal comparison charts...")

# Area represented by each grid cell
cell_area_km2 <- (cell_m / 1000)^2

# Calculate records, area with sightings, and records/km² by season
# Seasonal coverage summary
seasonal_cells <- cell_seasonal %>%
  dplyr::group_by(season) %>%
  dplyr::summarise(
    cells_with_data = sum(n_all > 0),
    area_with_data_km2 = cells_with_data * cell_area_km2,
    pct_study_area = area_with_data_km2 / study_area_km2 * 100,
    total_records = sum(n_all, na.rm = TRUE),
    records_per_km2 = total_records / area_with_data_km2,
    .groups = "drop"
  )

# ---- Plot 1: Total records by season -----------------------------------------

p_seasonal_bars <- ggplot(seasonal_cells, aes(x = season, y = total_records, fill = season)) +
  geom_col() +
  geom_text(
    aes(
      label = paste0(
        scales::comma(total_records), " records\n",
        scales::comma(round(area_with_data_km2)), " km² with sightings"
      )
    ),
    vjust = -0.5,
    size = 3.5
  ) +
  scale_fill_manual(values = SEASON_PALETTE) +
  labs(
    title = "Seasonal Distribution of Cetacean Sightings",
    subtitle = paste0(year_min, "–", year_max),
    x = "Season",
    y = "Total sighting records",
    caption = paste0(
      "Area is the approximate area with at least one sighting record, based on ",
      grid_km, " km × ", grid_km, " km grid cells."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  )

ggsave(
  file.path(out_dir, "00_seasonal_comparison_records.png"),
  p_seasonal_bars,
  width = 8,
  height = 6,
  dpi = 300
)

message("  ✓ Saved seasonal comparison chart: records")


# ---- Plot 2: Records per km² with sightings by season -------------------------
p_seasonal_rate <- ggplot(seasonal_cells, aes(x = season, y = records_per_km2, fill = season)) +
  geom_col() +
  
  # n records inside the bar
  geom_text(
    aes(
      y = records_per_km2 / 2,
      label = paste0("n = ", scales::comma(total_records))
    ),
    color = "white",
    fontface = "bold",
    size = 3.6
  ) +
  
  # rate and area above the bar
  geom_text(
    aes(
      label = paste0(
        round(records_per_km2, 3), " records/km²\n",
        scales::comma(round(area_with_data_km2)), " km²"
      )
    ),
    vjust = -0.5,
    size = 3.5
  ) +
  
  scale_fill_manual(values = SEASON_PALETTE) +
  labs(
    title = "Concentration of Cetacean Sighting Records",
    subtitle = paste0(year_min, "–", year_max),
    x = "Season",
    y = "Sighting records per km² with sightings",
    caption = paste0(
      "Sighting records/km² are calculated using grid cells (",
      grid_km, " km × ", grid_km,
      " km) with at least one sighting."
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  ) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.22))
  )

ggsave(
  file.path(out_dir, "00_seasonal_comparison_records_per_km2.png"),
  p_seasonal_rate,
  width = 8,
  height = 6,
  dpi = 300
)

message("  ✓ Saved seasonal comparison chart: records per km²")

#Plot 3: Seasonal coverage----------
p_seasonal_coverage <- ggplot(
  seasonal_cells,
  aes(x = season, y = pct_study_area, fill = season)
) +
  geom_col() +
  
  # n records inside the bar
  geom_text(
    aes(
      y = pct_study_area / 2,
      label = paste0("n = ", scales::comma(total_records))
    ),
    color = "white",
    fontface = "bold",
    size = 3.6
  ) +
  
  # area and percent above the bar
  geom_text(
    aes(
      label = paste0(
        round(pct_study_area, 1), "%\n",
        scales::comma(round(area_with_data_km2)), " km²"
      )
    ),
    vjust = -0.5,
    size = 3.5
  ) +
  
  scale_fill_manual(values = SEASON_PALETTE) +
  labs(
    title = "Seasonal Coverage of Cetacean Sighting Records",
    subtitle = paste0(year_min, "–", year_max, "|", " Study area = ",
      scales::comma(round(study_area_km2)), " km²"),
    x = "Season",
    y = "Study area with at least one sighting (%)",
    caption = paste0(
      "Proportion of the total study area with ≥ 1 sighting in a ",
      grid_km, "km × ", grid_km,
      " km cell."
          )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = .5)
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.18))
  )

ggsave(
  file.path(out_dir, "00_seasonal_coverage_pct_study_area.png"),
  p_seasonal_coverage,
  width = 8,
  height = 6,
  dpi = 300
)

message("  ✓ Saved seasonal coverage chart: percent of study area")

# ---------------------------
# SUMMARY STATISTICS
# ---------------------------
message("\n=== SEASONAL ANALYSIS COMPLETE ===")
message("\nOverall seasonal distribution:")
print(seasonal_summary)

message("\nKey findings:")
message(sprintf("  - Most sampled season: %s (%.1f%%)", 
                seasonal_summary$season[which.max(seasonal_summary$pct)],
                max(seasonal_summary$pct, na.rm = TRUE)))
message(sprintf("  - Least sampled season: %s (%.1f%%)", 
                seasonal_summary$season[which.min(seasonal_summary$pct)],
                min(seasonal_summary$pct, na.rm = TRUE)))

# Calculate ratio
max_season_records <- max(seasonal_summary$n_records, na.rm = TRUE)
min_season_records <- min(seasonal_summary$n_records, na.rm = TRUE)
ratio <- round(max_season_records / min_season_records, 1)

message(sprintf("  - Sampling ratio (max/min): %.1fx more records in peak season", ratio))

message("\nFiles created:")
message("  - 00_seasonal_coverage.png (faceted map)")
message("  - 00_seasonal_comparison.png (bar chart)")
if (!is.null(WEA_wind)) {
  message("  - osw_seasonal_summary.csv")
}
