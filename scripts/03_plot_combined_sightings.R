# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-03-11
# Script: 03_plot_combined_sightings.R
# Description: Plot cetacean sightings from two sources — WSDB (vessel) and
#              DFO aerial surveys — using a single pre-harmonised combined
#              CSV (combined_dedup_1km_day.csv). Records are
#              filtered by the 'source' column. Produces per-species maps and
#              a facet map for each source. Palettes and shapes are built once
#              and shared across per-species and facet maps so colour/shape are
#              consistent for a given species within each source.
# Changes from previous version:
#   - Section 6: replaced Spectral ramp with explicit WHALE_PALETTE keyed by
#     common name (no Latin initials). Covers baleen, beaked, and
#     dolphins/odontocetes. Keep in sync with 07_summarize_pam_detections.R.
#   - BUGFIX: Input coords corrected from x_m/y_m (UTM) to lon/lat (WGS84)
#   - BUGFIX: Added source value diagnostic block (section 5b)
#   - BUGFIX: Normalised source column to lowercase at load time
#   - BUGFIX: Guard clauses added before both facet_wrap calls
#   - BUGFIX: Corrected section 12 cat() — was printing wsdb_sf levels
#   - BUGFIX: ns_last() defined once in section 6, reused throughout
# ==============================================================================


# ---- 1. Load spatial helpers and packages ---------------------------------------

suppressWarnings(
  source(here::here("scripts/00_load_helpers.R"))
)
pacman::p_load(
  sf, tidyverse, here, scales, ggspatial, RColorBrewer
)

if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(
    dplyr::filter,
    dplyr::select,
    dplyr::mutate,
    dplyr::arrange,
    dplyr::summarise,
    dplyr::rename,
    lubridate::year,
    lubridate::month,
    .quiet = TRUE
  )
}


# ---- 2. Input file ---------------------------------------------------------------

combined_file <- "combined_dedup_1km_day.csv"

spatial_layers <- load_core_spatial_layers()
ess_study_final <- spatial_layers$study_area
WEA <- spatial_layers$wea
land <- spatial_layers$land
cont <- spatial_layers$contours


# ---- 3. Helpers ------------------------------------------------------------------

safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

# Sort species so NS groups always appear last in legends/facets
ns_last <- function(x) c(sort(x[!grepl("\\(NS\\)", x)]), sort(x[grepl("\\(NS\\)", x)]))


# ---- 4. Map extent ---------------------------------------------------------------

study_extent <- map_study_bbox(ess_study_final)
xlims_utm    <- study_extent$xlim
ylims_utm    <- study_extent$ylim


# ---- 5. Load combined sightings --------------------------------------------------

combined_raw <- read_csv(
  here("output/data/", combined_file),
  col_types = cols(.default = "c"),
  locale    = locale(encoding = "latin1")
) %>%
  mutate(
    lon    = as.numeric(lon),
    lat    = as.numeric(lat),
    source = tolower(trimws(source))
  ) %>%
  filter(!is.na(lon), !is.na(lat))

combined_sf <- st_as_sf(combined_raw, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM20) %>%
  st_make_valid() %>%
  mutate(year = as.integer(substr(date_utc, 1, 4))) %>%
  filter(year >= 2010) %>%
  filter(!grepl("Seal", common_name, ignore.case = TRUE)) %>%
  mutate(common_name = case_when(
    grepl("(NS)", common_name, fixed = TRUE) &
      grepl("Whale|Cetacean|Baleen|Mesoplodont|Beaked", common_name, ignore.case = TRUE) ~ "Cetacean/Whale (NS)",
    grepl("(NS)", common_name, fixed = TRUE) &
      grepl("Dolphin|Porpoise", common_name, ignore.case = TRUE) ~ "Dolphin/Porpoise (NS)",
    TRUE ~ common_name
  ))


# ---- 5b. DIAGNOSTIC: confirm source values ---------------------------------------

cat("\n--- SOURCE DIAGNOSTIC (section 5b) ---\n")
cat("Unique values in 'source' column (after tolower):\n")
print(table(combined_sf$source))
cat("Expected: 'wsdb' and 'aerial'\n")
cat("--------------------------------------\n\n")

if (!any(combined_sf$source == "wsdb") && !any(combined_sf$source == "aerial")) {
  stop(
    "Neither 'wsdb' nor 'aerial' found in source column after normalisation. ",
    "Check the table printed above and update the source filter strings."
  )
}


# ---- 5a. Front-end data summary --------------------------------------------------

n_wsdb_raw   <- sum(combined_sf$source == "wsdb")
n_aerial_raw <- sum(combined_sf$source == "aerial")
aerial_years = combined_sf%>%filter(source == "aerial")%>%distinct(year)
wsdb_sp_all   <- combined_sf %>% st_drop_geometry() %>%
  filter(source == "wsdb")   %>% distinct(common_name) %>% pull(common_name)
aerial_sp_all <- combined_sf %>% st_drop_geometry() %>%
  filter(source == "aerial") %>% distinct(common_name) %>% pull(common_name)
all_sp <- base::union(wsdb_sp_all, aerial_sp_all)

cat("\n=================================================================\n")
cat("  DATA SUMMARY — combined_dedup_1km_day.csv\n")
cat("=================================================================\n")
cat(sprintf("  Year range:          %d – %d\n", min(combined_sf$year), max(combined_sf$year)))
cat(sprintf("  Total records:       %d\n",       nrow(combined_sf)))
cat(sprintf("    WSDB (vessel):     %d\n",       n_wsdb_raw))
cat(sprintf("    Aerial:            %d\n",       n_aerial_raw))
cat(sprintf("    Aerial years:         %d - %d – %d\n",  nrow(aerial_years), min(aerial_years$year), max(aerial_years$year))   )
cat(sprintf("  Species/groups:\n"))
cat(sprintf("    WSDB:      %d  — %s\n", length(wsdb_sp_all),
            paste(sort(wsdb_sp_all), collapse = ", ")))
cat(sprintf("    Aerial:    %d  — %s\n", length(aerial_sp_all),
            paste(sort(aerial_sp_all), collapse = ", ")))
cat(sprintf("    Combined unique:   %d\n", length(all_sp)))
cat(sprintf("    Shared:            %d\n", length(base::intersect(wsdb_sp_all, aerial_sp_all))))
cat(sprintf("    WSDB only:         %d\n", length(base::setdiff(wsdb_sp_all, aerial_sp_all))))
cat(sprintf("    Aerial only:       %d\n", length(base::setdiff(aerial_sp_all, wsdb_sp_all))))
cat("=================================================================\n\n")


# ---- 6. Shared species palette and shapes ---------------------------------------
# WHALE_PALETTE is sourced from scripts/helpers/helper_palettes.R.
# NS groups use a generated grey-to-brown ramp plus cross/x shapes.

base_shapes <- c(21, 22, 23, 24, 25, 15, 8, 20, 18)

all_species   <- ns_last(unique(combined_sf$common_name))
named_species <- all_species[!grepl("\\(NS\\)", all_species)]
ns_species    <- all_species[grepl("\\(NS\\)", all_species)]
n_ns          <- length(ns_species)

# Warn if any species in the data are not covered — add them to WHALE_PALETTE
missing_spp <- named_species[!named_species %in% names(WHALE_PALETTE)]
if (length(missing_spp) > 0) {
  warning("Species missing from WHALE_PALETTE — add them: ",
          paste(missing_spp, collapse = ", "))
}

# Named species pull directly from WHALE_PALETTE by common name
named_colours <- WHALE_PALETTE[named_species]
ns_colours    <- setNames(
  colorRampPalette(c("grey50", "burlywood4"))(n_ns),
  ns_species
)

# Filled shapes for identified species; cross/x for NS groups
named_shapes <- setNames(rep(base_shapes, length.out = length(named_species)), named_species)
ns_shapes    <- setNames(rep(c(3, 4),     length.out = n_ns),                  ns_species)

shared_pal    <- c(named_colours, ns_colours)
shared_shapes <- c(named_shapes,  ns_shapes)


# ==============================================================================
# ---- 7. WSDB — filter, palette, shapes -------------------------------------------
# ==============================================================================

wsdb_threshold <- 25

wsdb_sf <- combined_sf %>%
  filter(source == "wsdb") %>%
  add_count(common_name) %>%
  filter(n >= wsdb_threshold) %>%
  select(-n)

wsdb_species <- ns_last(sort(unique(wsdb_sf$common_name)))
n_wsdb       <- length(wsdb_species)

cat("\nWSDB species retained (>=", wsdb_threshold, "records):\n")
print(wsdb_sf %>% st_drop_geometry() %>% count(common_name, sort = TRUE))

wsdb_pal    <- shared_pal[wsdb_species]
wsdb_shapes <- shared_shapes[wsdb_species]


# ---- 8. WSDB per-species maps ----------------------------------------------------

dir.create(here("output/figs/sights"), showWarnings = FALSE, recursive = TRUE)

for (sp in wsdb_species) {
  
  sp_data <- wsdb_sf %>% filter(as.character(common_name) == sp)
  if (nrow(sp_data) == 0) next
  
  sci_name <- sp_data$scientific_name[!is.na(sp_data$scientific_name)][1]
  sci_name <- if (is.na(sci_name) || length(sci_name) == 0) "Species not identified" else sci_name
  n_sp     <- nrow(sp_data)
  
  subtitle_text <- if (grepl("NS", sp)) {
    paste0(sci_name, " | WSDB opportunistic data (n = ", n_sp, ")")
  } else {
    bquote(italic(.(sci_name)) ~ " | WSDB opportunistic data (n = " ~ .(n_sp) ~ ")")
  }
  
  Map_sp <- ggplot() +
    map_layer_bathy_contours(cont) +
    map_layer_land(land) +
    map_layer_wea(WEA) +
    map_layer_study_area(ess_study_final) +
    geom_sf(
      data  = sp_data,
      aes(color = common_name, fill = common_name, shape = common_name),
      size  = 2, stroke = 0.5, alpha = 0.65
    ) +
    scale_color_manual(values = wsdb_pal[sp],    name = "") +
    scale_fill_manual( values = wsdb_pal[sp],    name = "") +
    scale_shape_manual(values = wsdb_shapes[sp], name = "") +
    labs(title = sp, subtitle = subtitle_text) +
    map_base_theme(base_size = 12, show_axis_text = TRUE, legend_position = "none") +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(hjust = 0.75)) +
    map_coord(xlim = xlims_utm, ylim = ylims_utm) +
    annotation_scale(location = "br", width_hint = 0.25)
  
  ggsave(
    here("output/figs/sights", paste0("WSDB_Sightings_", safe_filename(sp), ".png")),
    Map_sp, height = 7.5, width = 10, units = "in", dpi = 300
  )
}


# ---- 9. WSDB facet map -----------------------------------------------------------

if (nrow(wsdb_sf) == 0 || length(wsdb_species) == 0) {
  message("WSDB facet map skipped — no species met the threshold of ", wsdb_threshold, " records.")
} else {
  
  wsdb_sf <- wsdb_sf %>%
    mutate(common_name = factor(common_name, levels = wsdb_species))
  
  cat("WSDB factor levels:\n"); print(levels(wsdb_sf$common_name))
  
  Map_WSDB_Facet <- ggplot() +
    map_layer_bathy_contours(cont) +
    map_layer_land(land) +
    map_layer_wea(WEA) +
    map_layer_study_area(ess_study_final, linewidth = 0.4) +
    geom_sf(
      data  = wsdb_sf,
      aes(color = common_name, fill = common_name, shape = common_name),
      size  = 1.5, stroke = 0.3, alpha = 0.65
    ) +
    scale_color_manual(values = wsdb_pal,    name = "") +
    scale_fill_manual( values = wsdb_pal,    name = "") +
    scale_shape_manual(values = wsdb_shapes, name = "") +
    facet_wrap(~common_name, labeller = labeller(
      common_name = function(x) stringr::str_to_title(tolower(x))
    )) +
    labs(title = sprintf("WSDB sightings by cetacean species (n ≥ %d records)", wsdb_threshold)) +
    map_base_theme(base_size = 11, legend_position = "none") +
    xlab("") + ylab("") +
    theme(
      strip.text      = element_text(size = 8)
    ) +
    map_coord(xlim = xlims_utm, ylim = ylims_utm)
  
  Map_WSDB_Facet
  
  ggsave(
    here("output/figs/sights/WSDB_Whale_sightings_facet.png"),
    Map_WSDB_Facet,
    height = 10, width = 14, units = "in", dpi = 300
  )
}


# ==============================================================================
# ---- 10. Aerial — filter, palette, shapes ----------------------------------------
# ==============================================================================

aerial_threshold <- 10

aerial_sf <- combined_sf %>%
  filter(source == "aerial") %>%
  add_count(common_name) %>%
  filter(n >= aerial_threshold) %>%
  select(-n)

aerial_species <- ns_last(sort(unique(aerial_sf$common_name)))
n_aerial       <- length(aerial_species)

cat("\nAerial species retained (>=", aerial_threshold, "records):\n")
print(aerial_sf %>% st_drop_geometry() %>% count(common_name, sort = TRUE))

aerial_pal    <- shared_pal[aerial_species]
aerial_shapes <- shared_shapes[aerial_species]


# ---- 11. Aerial per-species maps -------------------------------------------------

for (sp in aerial_species) {
  
  sp_data <- aerial_sf %>% filter(as.character(common_name) == sp)
  if (nrow(sp_data) == 0) next
  
  sci_name <- sp_data$scientific_name[!is.na(sp_data$scientific_name)][1]
  sci_name <- if (is.na(sci_name) || length(sci_name) == 0) "Species not identified" else sci_name
  n_sp     <- nrow(sp_data)
  
  subtitle_text <- if (grepl("NS", sp)) {
    paste0(sci_name, " | DFO systematic surveys (n = ", n_sp, ")")
  } else {
    bquote(italic(.(sci_name)) ~ " | DFO systematic surveys (n = " ~ .(n_sp) ~ ")")
  }
  
  Map_aerial_sp <- ggplot() +
    map_layer_bathy_contours(cont) +
    map_layer_land(land) +
    map_layer_wea(WEA) +
    map_layer_study_area(ess_study_final) +
    geom_sf(
      data  = sp_data,
      aes(color = common_name, fill = common_name, shape = common_name),
      size  = 2, stroke = 0.4, alpha = 0.7
    ) +
    scale_color_manual(values = aerial_pal[sp],    name = "") +
    scale_fill_manual( values = aerial_pal[sp],    name = "") +
    scale_shape_manual(values = aerial_shapes[sp], name = "") +
    labs(title = sp, subtitle = subtitle_text) +
    map_base_theme(base_size = 12, show_axis_text = TRUE, legend_position = "none") +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(hjust = 0.75)) +
    map_coord(xlim = xlims_utm, ylim = ylims_utm) +
    annotation_scale(location = "br", width_hint = 0.25)
  
  ggsave(
    here("output/figs/sights", paste0("Aerial_Sightings_", safe_filename(sp), ".png")),
    Map_aerial_sp, height = 7.5, width = 10, units = "in", dpi = 300
  )
}


# ---- 12. Aerial facet map --------------------------------------------------------

if (nrow(aerial_sf) == 0 || length(aerial_species) == 0) {
  message("Aerial facet map skipped — no species met the threshold of ", aerial_threshold, " records.")
} else {
  
  aerial_sf <- aerial_sf %>%
    mutate(common_name = factor(common_name, levels = aerial_species))
  
  cat("Aerial factor levels:\n"); print(levels(aerial_sf$common_name))
  
  Map_Aerial_Facet <- ggplot() +
    map_layer_bathy_contours(cont) +
    map_layer_land(land) +
    map_layer_wea(WEA) +
    map_layer_study_area(ess_study_final, linewidth = 0.4) +
    geom_sf(
      data  = aerial_sf,
      aes(color = common_name, fill = common_name, shape = common_name),
      size  = 1.5, stroke = 0.3, alpha = 0.65
    ) +
    scale_color_manual(values = aerial_pal,    name = "") +
    scale_fill_manual( values = aerial_pal,    name = "") +
    scale_shape_manual(values = aerial_shapes, name = "") +
    facet_wrap(~common_name, labeller = labeller(
      common_name = function(x) stringr::str_to_title(tolower(x))
    )) +
    labs(title = sprintf("Aerial survey sightings by cetacean species (n ≥ %d records)", aerial_threshold)) +
    map_base_theme(base_size = 11, legend_position = "none") +
    xlab("") + ylab("") +
    theme(
      strip.text      = element_text(size = 8)
    ) +
    map_coord(xlim = xlims_utm, ylim = ylims_utm)
  
  Map_Aerial_Facet
  
  ggsave(
    here("output/figs/sights/Aerial_Whale_sightings_facet.png"),
    Map_Aerial_Facet,
    height = 10, width = 14, units = "in", dpi = 300
  )
}


# ---- 13. Save threshold RDS ------------------------------------------------------
saveRDS(
  list(wsdb = wsdb_threshold, aerial = aerial_threshold),
  here("output/data/sightings_thresholds.rds")
)
