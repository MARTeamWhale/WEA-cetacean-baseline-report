# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-03-10
# Script: 02_prepare_combined_sightings.R
# Description: Combines WSDB vessel and DFO aerial sightings into a single
#              deduplicated, study-area-clipped dataset for downstream analysis.
# Deduplication rules:
#   - WSDB only: max 1 record per species per day per 1 km grid cell
#   - Aerial: untouched — pre-deduplicated at source
#   - Cross-source duplicates (same cell/day) are retained from both sources
# Key outputs:
#   - combined_dedup_1km_day.csv       — analysis-ready combined dataset
#   - deduplication_summary.csv        — WSDB duplicate groups removed
#   - deduplication_audit_summary.csv  — compact duplicate audit by source type
#   - sightings_summary.csv            — record counts for Rmd reporting
#   - sightings_by_source.csv          — per-source summary for Rmd reporting
#   - sightings_by_species.csv         — per-species summary for Rmd reporting
# Changes from previous version:
#   - Added family (Odontocete / Mysticete) lookup hard-coded by
#     species_cd; appended to combined_dedup_1km_day.csv output
# ==============================================================================

pacman::p_load(dplyr, stringr, lubridate, readr, sf, here, tidyr)
suppressWarnings(source(here::here("scripts/00_load_helpers.R")))

# ==============================================================================
# CONFIGURATION ----------------------------------------------------------------
# ==============================================================================

AERIAL_PATH        <- "input/processed_data/aerial/deduplicated_sightings.csv"
WSDB_PATH          <- "input/raw_data/wsdb/WSDB_DFO_Sep2025.csv"
SPECIES_CODES_PATH <- "input/raw_data/speciesCodes.csv"
SPECIES_LOOKUP_PATH <- "input/raw_data/species_lookup_aerial_to_wsdb.csv"
DEDUP_OUTPUT       <- "output/data/combined_dedup_1km_day.csv"
SUMMARY_OUTPUT     <- "output/data/deduplication_summary.csv"
AUDIT_SUMMARY_OUTPUT <- "output/data/deduplication_audit_summary.csv"

PROJ_CRS    <- 32620   # UTM Zone 20N
CELL_SIZE_M <- 1000    # 1 km dedup grid


# ==============================================================================
# HELPER FUNCTIONS -------------------------------------------------------------
# ==============================================================================

# Parse WSDB datetime: "2000-01-12 18:00.00.000000000" -> POSIXct
parse_wsdb_datetime <- function(x) {
  x_clean <- str_replace(as.character(x), "(\\d{2}:\\d{2})\\.(\\d{2}).*$", "\\1:\\2")
  suppressWarnings(ymd_hms(x_clean, tz = "UTC"))
}

# Parse HHMM integer time to "HH:MM:SS" string
parse_hhmm_time <- function(hhmm) {
  x  <- str_pad(str_trim(as.character(hhmm)), width = 4, side = "left", pad = "0")
  x[x == "NA" | is.na(x)] <- NA_character_
  paste0(str_sub(x, 1, 2), ":", str_sub(x, 3, 4), ":00")
}

# 1 km grid cell ID from projected coordinates
grid_cell_id <- function(x_m, y_m) {
  paste0(floor(x_m / CELL_SIZE_M), "_", floor(y_m / CELL_SIZE_M))
}

add_grid_cells <- function(df, crs = PROJ_CRS) {
  df_proj <- df %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
    st_transform(crs)
  coords <- st_coordinates(df_proj)
  df_proj %>%
    st_drop_geometry() %>%
    mutate(
      date_utc      = as.Date(datetime_utc, tz = "UTC"),
      x_m           = coords[, 1],
      y_m           = coords[, 2],
      grid_cell_1km = grid_cell_id(x_m, y_m)
    )
}

# Reformat DB-style names: "WHALE-BLUE" -> "Blue Whale", "DOLPHINS-COMMON" -> "Common Dolphin"
clean_common_name <- function(x) {
  stringr::str_replace(
    stringr::str_to_title(tolower(as.character(x))),
    "^(Whales?|Dolphins?|Porpoise|Seal)s?-\\s*(.+)$", "\\2 \\1"
  ) %>%
    stringr::str_replace("\\bWhales\\b",  "Whale") %>%
    stringr::str_replace("\\bDolphins\\b", "Dolphin") %>%
    stringr::str_replace_all("\\(Ns\\)", "(NS)") %>%
    # move mid-name (NS) to end: "Mesoplodont (NS) Whale" -> "Mesoplodont Whale (NS)"
    stringr::str_replace("^(.+?)\\s*\\(NS\\)\\s*(.+)$", "\\1 \\2 (NS)") %>%
    stringr::str_trim()
}

# ==============================================================================
# 1. LOAD REFERENCE DATA -------------------------------------------------------
# ==============================================================================

message("Loading reference data...")

# Species reference — clean names here so every downstream join inherits them
species_ref <- read_csv(SPECIES_CODES_PATH, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  rename_with(~str_replace_all(., "commonname",  "common_name")) %>%
  rename_with(~str_replace_all(., "scientif",    "scientific_name")) %>%
  mutate(
    species_cd  = as.character(species_cd),
    common_name = clean_common_name(common_name)
  )

# Aerial -> WSDB species code mapping
species_lookup <- read_csv(SPECIES_LOOKUP_PATH, show_col_types = FALSE) %>%
  rename_with(~str_replace(., "sp_code_aerial", "sp_code_raw")) %>%
  select(sp_code_raw, wsdb_species_cd) %>%
  mutate(across(everything(), as.character)) %>%
  filter(!is.na(wsdb_species_cd), wsdb_species_cd != "", wsdb_species_cd != "NA")

if (nrow(species_lookup) == 0) stop("No valid species mappings found in: ", SPECIES_LOOKUP_PATH)

message(sprintf("  ✓ %d species in reference | %d aerial->WSDB mappings",
                nrow(species_ref), nrow(species_lookup)))

# ==============================================================================
# 2. LOAD AND STANDARDIZE WSDB -------------------------------------------------
# ==============================================================================

message("Reading WSDB data...")

wsdb <- read_csv(WSDB_PATH, show_col_types = FALSE,
                 col_types = cols(.default = col_character())) %>%
  mutate(
    datetime_utc = parse_wsdb_datetime(WS_DATETIME_UTC),
    # fallback: rebuild from date + time if datetime_utc failed to parse
    datetime_utc = if_else(
      is.na(datetime_utc),
      suppressWarnings(ymd_hms(paste(WS_DATE, parse_hhmm_time(WS_TIME)), tz = "UTC")),
      datetime_utc
    )
  ) %>%
  transmute(
    source       = "wsdb",
    source_id    = as.character(WS_ID),
    species_cd   = as.character(SPECIES_CD),
    lon          = as.numeric(LONGITUDE),
    lat          = as.numeric(LATITUDE),
    datetime_utc = datetime_utc
  ) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(datetime_utc), !is.na(species_cd))

message(sprintf("  ✓ %d WSDB records loaded", nrow(wsdb)))

# ==============================================================================
# 3. LOAD AND STANDARDIZE AERIAL -----------------------------------------------
# ==============================================================================

message("Reading aerial data...")

aerial <- read_csv(AERIAL_PATH, show_col_types = FALSE,
                   locale = locale(encoding = "Windows-1252")) %>%
  transmute(
    source    = "aerial",
    source_id = coalesce(as.character(ID), as.character(`...1`)),
    sp_code_raw  = as.character(Sp_code),
    lon          = as.numeric(Long_sight),
    lat          = as.numeric(Lat_sight),
    datetime_utc = suppressWarnings(ymd_hms(Date_time_utc, tz = "UTC"))
  ) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(datetime_utc), !is.na(sp_code_raw)) %>%
  left_join(species_lookup, by = "sp_code_raw") %>%
  filter(!is.na(wsdb_species_cd)) %>%
  transmute(source, source_id, species_cd = wsdb_species_cd, lon, lat, datetime_utc)

message(sprintf("  ✓ %d aerial records loaded and mapped", nrow(aerial)))

# ==============================================================================
# 4. AUDIT DUPLICATES BEFORE DEDUPLICATION -------------------------------------
# ==============================================================================

message("Auditing duplicate groups before deduplication...")

combined_pre_dedup <- bind_rows(wsdb, aerial) %>%
  add_grid_cells() %>%
  mutate(dedupe_key = paste(species_cd, date_utc, grid_cell_1km, sep = "|"))

dedup_audit <- combined_pre_dedup %>%
  group_by(species_cd, date_utc, grid_cell_1km, dedupe_key) %>%
  summarise(
    n_records = n(),
    n_wsdb    = sum(source == "wsdb"),
    n_aerial  = sum(source == "aerial"),
    .groups   = "drop"
  ) %>%
  filter(n_records > 1) %>%
  mutate(
    duplicate_type = case_when(
      n_wsdb > 1  & n_aerial == 0 ~ "wsdb_internal",
      n_aerial > 1 & n_wsdb == 0  ~ "aerial_internal",
      n_wsdb >= 1 & n_aerial >= 1 ~ "cross_source",
      TRUE                        ~ "other"
    ),
    excess_records = case_when(
      duplicate_type == "wsdb_internal"   ~ n_wsdb - 1L,
      duplicate_type == "aerial_internal" ~ n_aerial - 1L,
      duplicate_type == "cross_source"    ~ n_records - 1L,
      TRUE                                ~ n_records - 1L
    )
  )

audit_summary <- dedup_audit %>%
  group_by(duplicate_type) %>%
  summarise(
    n_groups       = n(),
    excess_records = sum(excess_records),
    .groups        = "drop"
  ) %>%
  tidyr::complete(
    duplicate_type = c("wsdb_internal", "aerial_internal", "cross_source", "other"),
    fill = list(n_groups = 0L, excess_records = 0L)
  ) %>%
  mutate(note = case_when(
    duplicate_type == "wsdb_internal"   ~ "Removed by this script: WSDB only, max 1 record/species/day/1km cell",
    duplicate_type == "aerial_internal" ~ "Not removed: aerial is expected to be pre-deduplicated",
    duplicate_type == "cross_source"    ~ "Not removed: WSDB and aerial matches are retained as independent sources",
    TRUE                                ~ "Review if non-zero"
  ))

dir.create(dirname(AUDIT_SUMMARY_OUTPUT), showWarnings = FALSE, recursive = TRUE)
write_csv(audit_summary, AUDIT_SUMMARY_OUTPUT)

message(sprintf("  ✓ audit summary saved to %s", AUDIT_SUMMARY_OUTPUT))
print(audit_summary)

# ==============================================================================
# 5. DEDUPLICATE WSDB (max 1 record per species/day/1km cell) ------------------
# ==============================================================================

message("Deduplicating WSDB records...")

wsdb_gridded <- wsdb %>% add_grid_cells()

n_pre <- nrow(wsdb_gridded)

wsdb_dedup <- wsdb_gridded %>%
  arrange(species_cd, date_utc, grid_cell_1km, datetime_utc) %>%
  group_by(species_cd, date_utc, grid_cell_1km) %>%
  mutate(records_in_group = n(), is_duplicate = row_number() > 1) %>%
  slice(1) %>%
  ungroup()

n_removed <- n_pre - nrow(wsdb_dedup)
message(sprintf("  ✓ %d removed | %d retained (%.1f%% reduction)",
                n_removed, nrow(wsdb_dedup), 100 * n_removed / n_pre))

# Save dedup summary for records
dedup_summary <- wsdb_gridded %>%
  group_by(species_cd, date_utc, grid_cell_1km) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  filter(n_records > 1) %>%
  left_join(species_ref %>% select(species_cd, common_name), by = "species_cd") %>%
  arrange(desc(n_records))

write_csv(dedup_summary, SUMMARY_OUTPUT)
message(sprintf("  ✓ %d duplicate groups saved to %s", nrow(dedup_summary), SUMMARY_OUTPUT))

# ==============================================================================
# 6. COMBINE DEDUPED WSDB + AERIAL ---------------------------------------------
# ==============================================================================

message("Combining deduped WSDB and aerial records...")

combined <- bind_rows(
  wsdb_dedup,
  aerial %>% mutate(
    date_utc         = as.Date(datetime_utc, tz = "UTC"),
    x_m              = NA_real_,
    y_m              = NA_real_,
    grid_cell_1km    = NA_character_,
    records_in_group = NA_integer_,
    is_duplicate     = FALSE
  )
) %>%
  mutate(datetime_utc = as.POSIXct(datetime_utc, tz = "UTC")) %>%
  left_join(species_ref %>% select(species_cd, common_name, scientific_name, family),
            by = "species_cd") 

# Recode ambiguous (NS) names to a single consistent label
combined <- combined %>%
  mutate(common_name = case_when(
    str_detect(coalesce(common_name, ""), regex("whale.*\\(NS\\)|cetacean.*\\(NS\\)", ignore_case = TRUE)) ~ "Cetacean (NS)",
    TRUE ~ common_name
  ))

message(sprintf("  ✓ %d total records (WSDB: %d | Aerial: %d)",
                nrow(combined), sum(combined$source == "wsdb"),
                sum(combined$source == "aerial")))

# ==============================================================================
# 7. CLIP TO STUDY AREA --------------------------------------------------------
# ==============================================================================

message("Clipping to study area...")

study_area <- load_spatial_layer("study_area", crs = PROJ_CRS)

combined_clipped <- combined %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE) %>%
  st_transform(PROJ_CRS) %>%
  st_filter(study_area) %>%
  st_drop_geometry()

message(sprintf("  ✓ %d retained after clipping (%d removed)",
                nrow(combined_clipped), nrow(combined) - nrow(combined_clipped)))

# ==============================================================================
# 8. SAVE OUTPUTS --------------------------------------------------------------
# ==============================================================================

message("Saving outputs...")

write_csv(combined_clipped, DEDUP_OUTPUT)
message(sprintf("  ✓ %s", DEDUP_OUTPUT))

# ---- Rmd summary tables (all from combined_clipped: deduped + study-area clipped) ----
# Truncate to analysis period (>= 2010) — applied before KDEs, maps, and all summary exports
DATE_START <- as.Date("2010-01-01")   # Analysis period start (inclusive)

combined_report <- combined_clipped %>%
  filter(date_utc >= DATE_START, 
         !str_detect(coalesce(common_name, ""), regex("seal", ignore_case = TRUE)))

message(sprintf("  ✓ %d retained after date filter (>= 2010)", nrow(combined_report)))

# Overall summary — this is the dataset that feeds KDEs and maps
combined_report %>%
  summarise(
    n_records    = n(),
    n_species    = n_distinct(common_name),
    n_species_ns = n_distinct(common_name[str_detect(coalesce(common_name, ""), "\\(NS\\)")]),
    date_min     = min(date_utc, na.rm = TRUE),
    date_max     = max(date_utc, na.rm = TRUE)
  ) %>%
  write_csv("output/data/sightings_summary.csv")

# Per-source summary — post-deduplication, post-clip (i.e. what goes into KDEs/maps)
combined_report %>%
  group_by(source) %>%
  summarise(
    n_records    = n(),
    n_species    = n_distinct(common_name),
    n_species_ns = n_distinct(species_cd[str_detect(coalesce(common_name, ""), "\\(NS\\)")]),
    date_range   = paste(min(date_utc, na.rm = TRUE), "to",
                         max(date_utc, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  write_csv("output/data/sightings_by_source.csv")

combined_report %>%  
  filter(!grepl("Seal", common_name, ignore.case = TRUE)) %>%
  count(common_name, sort = TRUE) %>%
  slice_head(n = 30) %>%
  write_csv("output/data/sightings_by_species.csv")

message("  ✓ sightings_summary.csv")
message("  ✓ sightings_by_source.csv")

# ==============================================================================
# 9. FINAL SUMMARY -------------------------------------------------------------
# ==============================================================================

message("\n==============================================================================")
message("FINAL SUMMARY")
message("==============================================================================")
combined_clipped %>%
  group_by(source, family) %>%
  summarise(n_records = n(), n_species = n_distinct(species_cd),
            date_range = paste(min(date_utc, na.rm = TRUE), "to",
                               max(date_utc, na.rm = TRUE)),
            .groups = "drop") %>%
  print()
message(sprintf("NOTE: %d WSDB duplicates removed; aerial untouched; both clipped to study area.",
                n_removed))
message("==============================================================================")


# ---- Export parameters for report ----
dedup_params <- tibble::tribble(
  ~Parameter,    ~Value,                    ~Script,                 ~Notes,
  "CELL_SIZE_M", as.character(CELL_SIZE_M), "02_prepare_combined_sightings.R", "Grid for dedup (1 record/species/day/cell)",
  "PROJ_CRS",    as.character(PROJ_CRS),    "02_prepare_combined_sightings.R", "CRS used throughout (UTM Zone 20N)"
)
readr::write_csv(dedup_params, here::here("output/data/params_dedup.csv"))
