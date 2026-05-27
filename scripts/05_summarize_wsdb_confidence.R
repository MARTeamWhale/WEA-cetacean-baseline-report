# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-03-12
# Script: 05_summarize_wsdb_confidence.R
# Description: Calculates data confidence summaries for target cetacean species
#              groups within WEA grid cells and the broader study region. Sources
#              objects produced by 04_summarize_wsdb_effort.R (cell_all, dt, targets,
#              target_patterns, WEA_wind). Outputs CSV tables and visualizations
#              to output/figs/Effort_maps/confidence_summaries/.
#
# IMPORTANT: Run 04_summarize_wsdb_effort.R first. This script depends on the following
#            objects remaining in the R environment:
#              cell_all        - gridded confidence data frame
#              dt              - sightings data (data.table or tibble)
#              targets         - named list of regex patterns
#              target_patterns - compiled regex objects
#              target_titles   - human-readable species group names
#              WEA_wind        - sf object for WEA polygons
#              out_dir         - output directory path
#              cell_m          - grid cell size in metres
#              use_datatable   - logical flag for data.table aggregation
#
# Changes:

#   - Cell polygons now correctly centred on xc/yc using +/- cell_m/2.
#   - Added wea_overlap_threshold parameter to configuration block (Section 1).
# ==============================================================================


# ---- 1. Dependency checks --------------------------------------------------------
# All objects below must exist from a completed 04_summarize_wsdb_effort.R run.
# If any are missing, stop with a clear message rather than failing silently.

message("\n=== Confidence Summary Analysis ===")
message("Checking for required objects from 04_summarize_wsdb_effort.R ...\n")

required_objects <- c("cell_all", "dt", "targets", "target_patterns",
                      "target_titles", "WEA_wind", "out_dir", "cell_m",
                      "use_datatable")

missing_objects <- required_objects[!sapply(required_objects, exists)]

if (length(missing_objects) > 0) {
  stop(
    "The following objects are missing. Please run 04_summarize_wsdb_effort.R first:\n  ",
    paste(missing_objects, collapse = ", ")
  )
}

message("  All required objects found.")


# ---- 1b. Configuration -----------------------------------------------------------
# Minimum proportion of a grid cell's area that must overlap the WEA polygon
# for the cell to be included in WEA summaries.
# 25% chosen based on diagnostic: genuine near-WEA cell (27_198) = 25.1% overlap,
# edge artefact cell (24_199) = 0.18% overlap. Threshold cleanly separates them.
wea_overlap_threshold <- 0.25


# ---- 2. Packages -----------------------------------------------------------------

required_pkgs <- c("dplyr", "sf", "tidyr", "readr", "ggplot2", "purrr", "stringr")
invisible(lapply(required_pkgs, library, character.only = TRUE))


# ---- 3. Output directory ---------------------------------------------------------

summary_dir <- file.path(out_dir, "confidence_summaries")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 4. Build grid cell polygons as sf -------------------------------------------
# Cells are constructed as squares centred on xc/yc (+/- cell_m/2).
# WEA cell selection uses a minimum overlap threshold (wea_overlap_threshold)
# rather than any-intersection, to exclude cells that only clip the WEA boundary
# by a tiny corner. Diagnostic confirmed 25% cleanly separates genuine near-WEA
# cells from edge artefacts.

message("\nConverting grid cells to spatial features (centred polygons) ...")

cells_sf <- make_cell_polygons(
  cells = cell_all,
  cell_m = cell_m,
  x_col = "xc",
  y_col = "yc",
  crs = 32620
)

message(sprintf("  Created %d cell polygons", nrow(cells_sf)))
message(sprintf("  Created %d cell polygons", nrow(cells_sf)))


# ---- 4b. WEA cell selection by minimum overlap -----------------------------------
# Calculates the proportion of each cell covered by WEA polygons.
# Only cells meeting wea_overlap_threshold are used for WEA summaries.
# This prevents edge cells with trivial WEA overlap from inflating counts.

get_wea_cells <- function(cells, wea_poly, threshold = wea_overlap_threshold) {
  wea_union  <- sf::st_union(wea_poly)
  cell_area  <- as.numeric(cell_m^2)
  
  # Intersection gives the portion of each cell inside the WEA
  inter <- sf::st_intersection(cells, wea_union) %>%
    mutate(overlap_area = as.numeric(sf::st_area(.)),
           overlap_pct  = overlap_area / cell_area) %>%
    sf::st_drop_geometry() %>%
    select(cell_id, overlap_pct)
  
  # Keep only cells meeting the threshold
  qualifying <- inter %>% filter(overlap_pct >= threshold) %>% pull(cell_id)
  
  message(sprintf(
    "  WEA cell selection: %d cells intersect WEA, %d meet >= %.0f%% overlap threshold",
    nrow(inter), length(qualifying), threshold * 100
  ))
  
  cells %>% filter(cell_id %in% qualifying)
}


# ---- 5. Helper: overall confidence metrics ---------------------------------------
# Summarises data coverage and temporal consistency across a set of grid cells.
# Works on any sf subset (full study region, WEA intersection, etc.).

calculate_confidence_metrics <- function(cells_data, region_name = "Region") {
  
  # Count cells by confidence category
  conf_counts <- cells_data %>%
    sf::st_drop_geometry() %>%
    count(confidence, .drop = FALSE) %>%
    mutate(
      percent = round(n / sum(n) * 100, 1),
      region  = region_name
    )
  
  # Aggregate summary statistics
  summary_stats <- cells_data %>%
    sf::st_drop_geometry() %>%
    summarise(
      region                     = region_name,
      total_cells                = n(),
      cells_with_data            = sum(n_all > 0),
      total_records              = sum(n_all, na.rm = TRUE),
      mean_records_per_cell      = round(mean(n_all[n_all > 0],   na.rm = TRUE), 1),
      median_records_per_cell    = round(median(n_all[n_all > 0], na.rm = TRUE), 1),
      mean_years_per_cell        = round(mean(n_years[n_years > 0],   na.rm = TRUE), 1),
      median_years_per_cell      = round(median(n_years[n_years > 0], na.rm = TRUE), 1),
      pct_high_confidence        = round(sum(confidence == "High",                         na.rm = TRUE) / n() * 100, 1),
      pct_medium_high_confidence = round(sum(confidence %in% c("High", "Medium-High"),     na.rm = TRUE) / n() * 100, 1),
      pct_low_confidence         = round(sum(confidence %in% c("Low", "Low-Medium"),       na.rm = TRUE) / n() * 100, 1)
    )
  
  list(counts = conf_counts, summary = summary_stats)
}


# ---- 6. Helper: target-species confidence metrics --------------------------------
# For a named target group (e.g. NARW), aggregates sightings per grid cell within
# a spatial region, then classifies each cell by evidence strength and overall
# data confidence. Returns counts and summary statistics.

calculate_target_metrics <- function(cells_data, target_name, target_regex,
                                     dt_data, region_name = "Region") {
  
  message(sprintf("  Processing %s ...", target_name))
  
  # Cell IDs present in this spatial region
  region_cell_ids <- cells_data %>%
    sf::st_drop_geometry() %>%
    pull(cell_id)
  
  # ---- 6a. Aggregate sightings per cell ----------------------------------------
  if (use_datatable && "data.table" %in% class(dt_data)) {
    
    target_summary <- dt_data[cell_id %in% region_cell_ids][, .(
      n_all    = sum(w, na.rm = TRUE),
      n_target = sum(w * stringr::str_detect(name_blob, target_regex), na.rm = TRUE)
    ), by = .(cell_id)]
    
    target_years <- dt_data[
      cell_id %in% region_cell_ids &
        stringr::str_detect(name_blob, target_regex),
      .(n_target_years = data.table::uniqueN(year, na.rm = TRUE)),
      by = .(cell_id)
    ]
    
    target_summary <- as_tibble(target_summary)
    target_years   <- as_tibble(target_years)
    
  } else {
    
    target_summary <- dt_data %>%
      filter(cell_id %in% region_cell_ids) %>%
      mutate(is_target = stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(
        n_all    = sum(w, na.rm = TRUE),
        n_target = sum(w[is_target], na.rm = TRUE),
        .groups  = "drop"
      )
    
    target_years <- dt_data %>%
      filter(cell_id %in% region_cell_ids,
             stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(n_target_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
  }
  
  # ---- 6b/6c. Join years, classify evidence, combine with cell confidence ------
  # All steps in one pipeline to avoid column collision from the redundant
  # second confidence join that previously existed in 6c.
  # Starting from cells_data ensures all region cells are in the denominator.
  
  target_summary <- target_summary %>%
    left_join(target_years, by = "cell_id") %>%
    mutate(n_target_years = tidyr::replace_na(n_target_years, 0))
  
  target_summary <- cells_data %>%
    sf::st_drop_geometry() %>%
    select(cell_id, confidence) %>%
    left_join(target_summary, by = "cell_id") %>%
    mutate(
      n_all          = tidyr::replace_na(n_all, 0),
      n_target       = tidyr::replace_na(n_target, 0),
      n_target_years = tidyr::replace_na(n_target_years, 0),
      # Cells with no cetacean data at all get confidence = NA from the left join
      # -- treat these as "No Data" for the target_conf classification
      confidence     = tidyr::replace_na(as.character(confidence), "No Data"),
      evidence = case_when(
        n_target == 0                        ~ "None",
        n_target >= 5 & n_target_years >= 3  ~ "High",
        n_target >= 2 & n_target_years >= 2  ~ "Medium",
        TRUE                                 ~ "Low"
      ),
      evidence = factor(evidence, levels = c("High", "Medium", "Low", "None")),
      target_conf = case_when(
        evidence == "None"                                                  ~ "No evidence",
        confidence %in% c("No Data", "Low")                                ~ "No evidence",
        confidence == "Low-Medium" & evidence == "Low"                     ~ "Low",
        confidence %in% c("Medium", "Medium-High") & evidence == "Medium" ~ "Moderate",
        confidence %in% c("Medium-High", "High") & evidence == "High"     ~ "High",
        TRUE                                                                ~ "Moderate"
      ),
      target_conf = factor(target_conf,
                           levels = c("High", "Moderate", "Low", "No evidence")))
  
  message(sprintf("  %s: n() in target_summary = %d (expect %d)",
                  target_name, nrow(target_summary),
                  nrow(sf::st_drop_geometry(cells_data))))
  
  # ---- 6d. Count cells by target confidence ------------------------------------
  conf_counts <- target_summary %>%
    count(target_conf, .drop = FALSE) %>%
    mutate(
      percent = round(n / sum(n) * 100, 1),
      target  = target_name,
      region  = region_name
    )
  message(sprintf("  %s: nrow before summarise = %d, n_target==0: %d, evidence None: %d, target_conf No evidence: %d",
                  target_name,
                  nrow(target_summary),
                  sum(target_summary$n_target == 0),
                  sum(target_summary$evidence == "None"),
                  sum(target_summary$target_conf == "No evidence")))
  
  # ---- 6e. Summary statistics --------------------------------------------------
  summary_stats <- target_summary %>%
    summarise(
      target                 = target_name,
      region                 = region_name,
      total_cells            = n(),
      cells_with_target      = sum(n_target > 0),
      total_target_records   = sum(n_target, na.rm = TRUE),
      mean_target_per_cell   = round(mean(n_target[n_target > 0],   na.rm = TRUE), 1),
      median_target_per_cell = round(median(n_target[n_target > 0], na.rm = TRUE), 1),
      pct_high_confidence    = round(sum(target_conf == "High",                  na.rm = TRUE) / n() * 100, 1),
      pct_moderate_plus      = round(sum(target_conf %in% c("High", "Moderate"), na.rm = TRUE) / n() * 100, 1),
      pct_low_confidence     = round(sum(target_conf == "Low",                   na.rm = TRUE) / n() * 100, 1),
      pct_no_evidence        = round(sum(target_conf == "No evidence",           na.rm = TRUE) / n() * 100, 1)
    )
  
  list(counts = conf_counts, summary = summary_stats)
}


# ---- 7. Study region metrics -----------------------------------------------------

message("\nCalculating study region confidence metrics ...")

study_metrics <- calculate_confidence_metrics(cells_sf, "Study Region")

message(sprintf("  Total cells:          %d",   study_metrics$summary$total_cells))
message(sprintf("  Cells with data:      %d (%.1f%%)",
                study_metrics$summary$cells_with_data,
                study_metrics$summary$cells_with_data /
                  study_metrics$summary$total_cells * 100))
message(sprintf("  High confidence:      %.1f%%", study_metrics$summary$pct_high_confidence))
message(sprintf("  Medium-High+:         %.1f%%", study_metrics$summary$pct_medium_high_confidence))

# Target-specific metrics — study region
study_target_metrics <- list()
for (nm in names(targets)) {
  study_target_metrics[[nm]] <- calculate_target_metrics(
    cells_sf, nm, target_patterns[[nm]], dt, "Study Region"
  )
}


# ---- 8. WEA area metrics ---------------------------------------------------------
# Uses get_wea_cells() (Section 4b) to select only cells with >= 25% overlap
# with the WEA polygons, excluding edge artefacts.

if (!is.null(WEA_wind)) {
  message("\nCalculating WEA confidence metrics ...")
  
  WEA_cells <- get_wea_cells(cells_sf, WEA_wind, threshold = wea_overlap_threshold)
  
  if (nrow(WEA_cells) > 0) {
    
    WEA_metrics <- calculate_confidence_metrics(WEA_cells, "WEA")
    
    message(sprintf("  WEA cells (>=%.0f%% overlap): %d", wea_overlap_threshold * 100,
                    WEA_metrics$summary$total_cells))
    message(sprintf("  Cells with data:             %d (%.1f%%)",
                    WEA_metrics$summary$cells_with_data,
                    WEA_metrics$summary$cells_with_data /
                      WEA_metrics$summary$total_cells * 100))
    message(sprintf("  High confidence:             %.1f%%", WEA_metrics$summary$pct_high_confidence))
    message(sprintf("  Medium-High+:                %.1f%%", WEA_metrics$summary$pct_medium_high_confidence))
    
    # Target-specific metrics — WEA areas
    WEA_target_metrics <- list()
    for (nm in names(targets)) {
      WEA_target_metrics[[nm]] <- calculate_target_metrics(
        WEA_cells, nm, target_patterns[[nm]], dt, "WEA Areas"
      )
    }
    
  } else {
    message("  WARNING: No cells intersect WEA polygons after correction.")
    message("  Check that cells_sf and WEA_wind share CRS 32620.")
    WEA_metrics       <- NULL
    WEA_target_metrics <- NULL
  }
  
} else {
  message("\nSkipping WEA summary (WEA_wind is NULL).")
  WEA_metrics        <- NULL
  WEA_target_metrics <- NULL
}


# ---- 9. Compile results ----------------------------------------------------------

message("\nCompiling results ...")

overall_conf_summary <- bind_rows(
  study_metrics$summary,
  if (!is.null(WEA_metrics)) WEA_metrics$summary else NULL
)

overall_conf_counts <- bind_rows(
  study_metrics$counts,
  if (!is.null(WEA_metrics)) WEA_metrics$counts else NULL
)

target_conf_summary <- bind_rows(
  lapply(names(study_target_metrics), function(nm) study_target_metrics[[nm]]$summary),
  if (!is.null(WEA_target_metrics))
    lapply(names(WEA_target_metrics), function(nm) WEA_target_metrics[[nm]]$summary)
  else NULL
)

target_conf_counts <- bind_rows(
  lapply(names(study_target_metrics), function(nm) study_target_metrics[[nm]]$counts),
  if (!is.null(WEA_target_metrics))
    lapply(names(WEA_target_metrics), function(nm) WEA_target_metrics[[nm]]$counts)
  else NULL
)


# ---- 10. Save CSV outputs --------------------------------------------------------

readr::write_csv(overall_conf_summary, file.path(summary_dir, "overall_confidence_summary.csv"))
readr::write_csv(overall_conf_counts,  file.path(summary_dir, "overall_confidence_counts.csv"))
readr::write_csv(target_conf_summary,  file.path(summary_dir, "target_confidence_summary.csv"))
readr::write_csv(target_conf_counts,   file.path(summary_dir, "target_confidence_counts.csv"))

# Also write flat copies for RMarkdown / downstream scripts
readr::write_csv(overall_conf_summary, here::here("output/data/overall_conf_summary.csv"))
readr::write_csv(target_conf_summary,  here::here("output/data/target_conf_summary.csv"))

message("  Saved: overall_confidence_summary.csv")
message("  Saved: overall_confidence_counts.csv")
message("  Saved: target_confidence_summary.csv")
message("  Saved: target_confidence_counts.csv")
message("  Saved: output/data/overall_conf_summary.csv (for RMarkdown)")
message("  Saved: output/data/target_conf_summary.csv  (for RMarkdown)")


# ---- 11. Visualisation: region comparison bar chart -----------------------------

message("\nCreating summary visualisations ...")

if (!is.null(WEA_metrics)) {
  
  p_comparison <- ggplot(overall_conf_counts,
                         aes(x = confidence, y = percent, fill = region)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Study Region" = "#2c7fb8", "WEA" = "#ff7f00")) +
    labs(
      title = "Data Confidence Comparison: Study Region vs WEA Areas",
      x     = "Confidence Category",
      y     = "Percentage of Grid Cells (%)",
      fill  = "Region"
    ) +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(
    file.path(summary_dir, "confidence_comparison.png"),
    p_comparison, width = 10, height = 6, dpi = 300
  )
  message("  Saved: confidence_comparison.png")
}


# ---- 12. Visualisation: target confidence heatmap (WEA) -------------------------
# Four confidence categories for target species within WEA cells:
#   High        = strong positive evidence (>=5 records, >=3 years)
#   Moderate+   = cumulative High + Moderate
#   Low         = positive records but poor data quality
#   No evidence = species never recorded in this cell (absence signal)

if (!is.null(WEA_target_metrics)) {
  
  p_target_heatmap <- target_conf_summary %>%
    filter(region == "WEA Areas") %>%
    select(target, pct_high_confidence, pct_moderate_plus,
           pct_low_confidence, pct_no_evidence) %>%
    tidyr::pivot_longer(
      cols      = starts_with("pct_"),
      names_to  = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = case_when(
        metric == "pct_high_confidence" ~ "High",
        metric == "pct_moderate_plus"   ~ "Moderate+",
        metric == "pct_low_confidence"  ~ "Low",
        metric == "pct_no_evidence"     ~ "No evidence"
      ),
      metric     = factor(metric, levels = c("High", "Moderate+", "Low", "No evidence")),
      color_cat  = ifelse(metric == "No evidence", "no_evidence", "has_data"),
      blue_value = case_when(
        value == 0            ~ NA_real_,
        metric == "High"      ~ value * 1.0,
        metric == "Moderate+" ~ value * 0.75,
        metric == "Low"       ~ value * 0.4,
        TRUE                  ~ NA_real_
      ),
      grey_value = ifelse(metric == "No evidence" & value > 0, value, NA_real_)
    ) %>%
    ggplot(aes(x = metric, y = target)) +
    geom_tile(data = ~filter(.x, color_cat == "has_data"),
              aes(fill = blue_value), color = "white") +
    geom_tile(data = ~filter(.x, color_cat == "no_evidence" & !is.na(grey_value)),
              aes(alpha = grey_value), fill = "grey60", color = "white") +
    geom_tile(data = ~filter(.x, value == 0),
              fill = "white", color = "grey80") +
    geom_text(aes(label = sprintf("%.1f%%", value)), size = 3.5) +
    scale_fill_gradient(
      low      = "#deebf7",
      high     = "#08519c",
      name     = "% cells",
      limits   = c(0, 100),
      na.value = "white"
    ) +
    scale_alpha_continuous(
      range    = c(0.2, 0.8),
      limits   = c(0, 100),
      na.value = 0,
      guide    = "none"
    ) +
    labs(
      title    = "Sightings data confidence: WEA areas",
      subtitle = paste0(
        "Percentage of WEA grid cells in each confidence category\n",
        "Blue = species recorded (darker = higher confidence); ",
        "Grey = no records for this species in WEA cells"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.text.y     = element_text(size = 10),
      plot.subtitle   = element_text(size = 9, color = "grey30"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(summary_dir, "target_confidence_heatmap_WEA.png"),
    p_target_heatmap, width = 11, height = 6, dpi = 300
  )
  message("  Saved: target_confidence_heatmap_WEA.png")
  
} else {
  message("  Skipping target heatmap (no WEA data available).")
}


# ---- 13. Print summary tables to console ----------------------------------------

message("\n=== SUMMARY TABLES ===\n")
message("OVERALL CONFIDENCE SUMMARY:")
print(overall_conf_summary, width = Inf)

message("\nTARGET-SPECIFIC CONFIDENCE SUMMARY:")
print(target_conf_summary, width = Inf)


# ---- 14. Completion message ------------------------------------------------------

message("\n", paste(rep("=", 70), collapse = ""))
message("CONFIDENCE SUMMARY COMPLETE")
message(paste(rep("=", 70), collapse = ""))
message("Output directory: ", summary_dir)
message("\nFiles created:")
message("  confidence_summaries/overall_confidence_summary.csv")
message("  confidence_summaries/overall_confidence_counts.csv")
message("  confidence_summaries/target_confidence_summary.csv")
message("  confidence_summaries/target_confidence_counts.csv")
message("  confidence_summaries/confidence_comparison.png")
message("  confidence_summaries/target_confidence_heatmap_WEA.png")
message("  output/data/overall_conf_summary.csv  (for RMarkdown)")
message("  output/data/target_conf_summary.csv   (for RMarkdown)")
message(paste(rep("=", 70), collapse = ""))
