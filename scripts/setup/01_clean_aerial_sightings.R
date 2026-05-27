# ==============================================================================
# Marine Mammal Aerial Survey Data Cleaning with Deduplication
# ==============================================================================
# Purpose: Load, check variables for NAs, deduplicate, and summarize aerial sighting data from
#          DFO systematic surveys using Observation_ID and double platform filtering to deduplicate
# ==============================================================================

# Load libraries ---------------------------------------------------------------
pacman::p_load(
  dplyr,
  lubridate,
  sf,
  ggplot2,
  patchwork,
  here,
  readr
)

# Load data --------------------------------------------------------------------
input_file <- "AllPlatforms_DFO_systematic_surveys_SS_BF_CS_Geoloc-20251205.csv"

# Read with UTF-8 (as indicated by file metadata)
whale_data <- read_csv(
  here("input/raw_data/aerial", input_file),
  col_types = cols(.default = "c"),
  locale = locale(encoding = "UTF-8")
)

# Data preparation -------------------------------------------------------------

## Parse datetime fields
df <- whale_data |>
  mutate(
    Date_time_utc = ymd_hms(Date_time_utc, tz = "UTC"),
    Date_time_loc = ymd_hms(Date_time_loc, quiet = TRUE)
  )

## Convert key numeric fields
df <- df |>
  mutate(
    Inclino = as.numeric(Inclino),
    Nb_tot  = as.numeric(Nb_tot)
  )

## Create Observation_id if missing (fallback to row number)
if (!"Observation_id" %in% names(df)) {
  df <- df |>
    mutate(Observation_id = row_number())
}

NA_date_aerial_MM = df%>%filter(is.na(Date_time_utc)) 

# write_csv(NA_date_aerial_MM, "NA_date_aerial_MM.csv")

# Check marine mammal species names ----------------

species_known <- whale_data |>
  group_by(Sp_code, Sp_name) |>
  summarise(
    n_sightings = n_distinct(Observation_id),
    n_animals   = sum(as.numeric(Nb_tot), na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_sightings))

cat("\nKnown species summary:\n")
print(species_known)

## Define marine mammal species codes (excludes fish, sharks, turtles, boats, etc.)
marine_mammal_codes <- c(
  # Baleen whales (Mysticetes)
  "Ba",    # Minke Whale
  "Bb",    # Sei Whale
  "Bm",    # Blue Whale
  "Bmys",  # Bowhead Whale
  "Bp",    # Fin Whale
  "Eg",    # North Atlantic Right Whale
  "Mn",    # Humpback Whale
  
  # Toothed whales (Odontocetes)
  "Dl",    # Beluga Whale
  "Pm",    # Sperm Whale
  "Kb",    # Pygmy Sperm Whale
  "Ha",    # Northern Bottlenose Whale
  "Mbi",   # Sowerby's Beaked Whale
  "Zc",    # Cuvier's Beaked Whale
  "Mm",    # Narwhal
  "Oo",    # Killer Whale
  "Gm",    # Pilot Whale
  
  # Dolphins
  "Dd",    # Common Dolphin
  "Gg",    # Risso's Dolphin
  "La",    # Atlantic White-sided Dolphin
  "Lalbi", # White-beaked Dolphin
  "Sc",    # Striped Dolphin
  "Sf",    # Atlantic Spotted Dolphin
  "Tt",    # Bottlenose Dolphin
  
  # Porpoises
  "Pp",    # Harbour Porpoise
  
  # Seals
  "Cc",    # Hooded Seal
  "Eb",    # Bearded Seal
  "Hg",    # Grey Seal
  "Pg",    # Harp Seal
  "Ph",    # Ringed Seal
  "Pv",    # Harbour Seal
  
  # Other marine mammals
  "Or",    # Walrus
  "Um",    # Polar Bear
  
  # Unknown marine mammals (keep these for analysis)
  "UnBe",  # Unknown beaked whale
  "UnCe",  # Unknown cetacean
  "UnDo",  # Unknown dolphin
  "UnDP",  # Unknown dolphin or porpoise
  "UnMe",  # Unknown mesoplodon
  "UnPo",  # Unknown porpoise
  "UnSe",  # Unknown seal
  "UnWh"   # Unknown whale
)

## Filter to marine mammals only
df <- df |>
  filter(Sp_code %in% marine_mammal_codes)

cat("Filtered to", nrow(df), "marine mammal records\n")

## Standardize ALL species names using Sp_code as the authoritative key
## This completely avoids encoding issues by ignoring Sp_name field
## Also consolidates related unknown species categories
df <- df |>
  mutate(
    Sp_name_clean = case_when(
      # Baleen whales
      Sp_code == "Ba" ~ "Minke Whale / Petit rorqual",
      Sp_code == "Bb" ~ "Sei Whale / Rorqual de Rudolphi",
      Sp_code == "Bm" ~ "Blue Whale / Rorqual bleu",
      Sp_code == "Bmys" ~ "Bowhead Whale / Baleine boreale",
      Sp_code == "Bp" ~ "Fin Whale / Rorqual commun",
      Sp_code == "Eg" ~ "North Atlantic Right Whale / Baleine noire",
      Sp_code == "Mn" ~ "Humpback Whale / Rorqual a bosse",
      
      # Toothed whales
      Sp_code == "Dl" ~ "Beluga Whale / Beluga",
      Sp_code == "Pm" ~ "Sperm Whale / Cachalot",
      Sp_code == "Kb" ~ "Pygmy Sperm Whale / Cachalot pygmee",
      Sp_code == "Ha" ~ "Northern Bottlenose Whale / Baleine a bec commune",
      Sp_code == "Mbi" ~ "Sowerby's Beaked Whale / Baleine a bec de Sowerby",
      Sp_code == "Zc" ~ "Cuvier's Beaked Whale / Baleine a bec de Cuvier",
      Sp_code == "Mm" ~ "Narwhal / Narval",
      Sp_code == "Oo" ~ "Killer Whale / Epaulard",
      Sp_code == "Gm" ~ "Pilot Whale / Globicephale",
      
      # Dolphins
      Sp_code == "Dd" ~ "Common Dolphin / Dauphin commun",
      Sp_code == "Gg" ~ "Risso's Dolphin / Dauphin de Risso",
      Sp_code == "La" ~ "Atlantic White-sided Dolphin / Dauphin a flanc blanc",
      Sp_code == "Lalbi" ~ "White-beaked Dolphin / Dauphin a nez blanc",
      Sp_code == "Sc" ~ "Striped Dolphin / Dauphin bleu et blanc",
      Sp_code == "Sf" ~ "Atlantic Spotted Dolphin / Dauphin tachete",
      Sp_code == "Tt" ~ "Bottlenose Dolphin / Grand dauphin",
      
      # Porpoises
      Sp_code == "Pp" ~ "Harbour Porpoise / Marsouin commun",
      
      # Seals
      Sp_code == "Cc" ~ "Hooded Seal / Phoque a capuchon",
      Sp_code == "Eb" ~ "Bearded Seal / Phoque barbu",
      Sp_code == "Hg" ~ "Grey Seal / Phoque gris",
      Sp_code == "Pg" ~ "Harp Seal / Phoque du Groenland",
      Sp_code == "Ph" ~ "Ringed Seal / Phoque annele",
      Sp_code == "Pv" ~ "Harbour Seal / Phoque commun",
      
      # Other marine mammals
      Sp_code == "Or" ~ "Walrus / Morse",
      Sp_code == "Um" ~ "Polar Bear / Ours polaire",
      
      # Unknown marine mammals - MERGED CATEGORIES
      # Merge UnBe + UnMe -> Unknown Beaked Whale
      Sp_code %in% c("UnBe", "UnMe") ~ "Unknown Beaked Whale / Baleine a bec inconnue",
      
      # Merge UnCe + UnWh -> Unknown Whale/Cetacean
      Sp_code %in% c("UnCe", "UnWh") ~ "Unknown Whale/Cetacean / Baleine/Cetace inconnu",
      
      # Merge UnDo + UnPo + UnDP -> Unknown Dolphin/Porpoise
      Sp_code %in% c("UnDo", "UnPo", "UnDP") ~ "Unknown Dolphin/Porpoise / Dauphin/Marsouin inconnu",
      
      # Unknown Seal (standalone)
      Sp_code == "UnSe" ~ "Unknown Seal / Phoque inconnu",
      
      # Fallback (should not happen if marine_mammal_codes is complete)
      TRUE ~ paste0(Sp_code, " / Unknown")
    )
  )

# Basic data exploration -------------------------------------------------------

## Number of unique species
n_species <- df |>
  summarise(n_species = n_distinct(Sp_code))

## Species lookup table
species_types <- df |>
  distinct(Sp_code, Sp_name)

## Survey date range (UTC)
date_range <- range(df$Date_time_utc, na.rm = TRUE)

## Per-species summary (raw data, marine mammals only)
species_summary_raw <- df |>
  group_by(Sp_code, Sp_name_clean) |>
  summarise(
    n_sightings = n_distinct(Observation_id),
    n_animals   = sum(Nb_tot, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_sightings))

# Print exploration results
cat("Number of species:", n_species$n_species, "\n")
cat("Date range:", as.character(date_range), "\n\n")
print(species_summary_raw)

# ==============================================================================
# DEDUPLICATION: Two-step approach based on email instructions
# ==============================================================================
#
# STEP 1: Remove duplicate Observation_IDs
#   - Observation_ID is assigned to the same individual seen multiple times
#   - This applies to double platform surveys where same animal seen by multiple observers
#   - Keep only the FIRST occurrence of each Observation_ID
#
# STEP 2: Handle double platform configurations (same-side duplicates)
#   - When Effort_type_l or Effort_type_r == "double", multiple observers on same side
#   - Keep only ONE observer per side when in double configuration
#   - DO NOT discard entire sides (would lose half the data)
#   - This addresses same-side observers recording the same individual twice
# ==============================================================================

# STEP 1: Remove duplicate Observation_IDs ------------------------------------
# IMPORTANT: Observation_ID is only assigned for double platform surveys
# where the same individual was seen multiple times. Most records will have NA.
# We only deduplicate records that have a VALID (non-NA) Observation_ID.

# First, check how many records have valid Observation_IDs
n_with_obs_id <- sum(!is.na(df$Observation_id)) #most dont have obsID
n_na_obs_id <- sum(is.na(df$Observation_id))

cat("\nObservation_ID distribution:\n")
cat("  Records WITH Observation_ID:", n_with_obs_id, "\n")
cat("  Records with NA Observation_ID:", n_na_obs_id, "\n")

# Separate data into two groups: with and without Observation_ID
df_with_id <- df |>
  filter(!is.na(Observation_id)) |>
  arrange(Survey, Date_time_utc, Observation_id) |>
  group_by(Observation_id) |>
  slice(1) |>  # Keep first occurrence of each Observation_ID
  ungroup()

df_no_id <- df |>
  filter(is.na(Observation_id))

# Combine back together
df_step1 <- bind_rows(df_with_id, df_no_id) |>
  arrange(Survey, Date_time_utc)

# Calculate removals from Step 1
n_removed_step1 <- n_with_obs_id - nrow(df_with_id)
cat("\nStep 1 - Observation_ID deduplication:\n")
cat("  Duplicate records removed:", n_removed_step1, "\n")
cat("  Records remaining:", nrow(df_step1), "\n")
cat("  (Note: Only records with non-NA Observation_ID were evaluated)\n")

# STEP 2: Remove secondary observers from double platform configurations ------
# For double platform configs: keep only primary observer records, remove secondary (rear seat position)

## First, let's examine the double platform configuration
cat("\n=== STEP 2 DIAGNOSTIC for double platforms ===\n")

# Check how many records have double platform on each side
n_double_left <- sum(!is.na(df_step1$Effort_type_l) & df_step1$Effort_type_l == "double", na.rm = TRUE)
n_double_right <- sum(!is.na(df_step1$Effort_type_r) & df_step1$Effort_type_r == "double", na.rm = TRUE)

cat("Records with double platform LEFT:", n_double_left, "\n")
cat("Records with double platform RIGHT:", n_double_right, "\n")

## Determine which side each observer is on
df_step2 <- df_step1 |>
  mutate(
    obs_side = case_when(
      grepl("l", Observer_pos, ignore.case = TRUE) | 
        grepl("left|babord|port", Observer_pos2, ignore.case = TRUE) ~ "left",
      grepl("r", Observer_pos, ignore.case = TRUE) | 
        grepl("right|tribord|starboard", Observer_pos2, ignore.case = TRUE) ~ "right",
      Observer_pos == "c" | Observer_pos == "AllPositions" ~ "center",
      TRUE ~ "unknown"
    ),
    # Is this side in double platform configuration?
    side_is_double = case_when(
      obs_side == "left" & (!is.na(Effort_type_l) & Effort_type_l == "double") ~ TRUE,
      obs_side == "right" & (!is.na(Effort_type_r) & Effort_type_r == "double") ~ TRUE,
      TRUE ~ FALSE
    ),
    # Is this observer NOT in front position?
    is_secondary = !is.na(Observer_pos2) & Observer_pos2 != "f"
  )

# Diagnostic: show what's being removed
non_front_summary <- df_step2 |>
  filter(side_is_double & is_secondary) |>
  group_by(Survey, obs_side) |>
  summarise(
    n_non_front_removed = n(),
    observers = paste(unique(Observer), collapse = ", "),
    observer_positions = paste(unique(Observer_pos), collapse = ", "),
    observer_pos2_values = paste(unique(Observer_pos2), collapse = ", "),
    .groups = "drop"
  )

cat("\n=== Non-front observer records to be removed from double platforms ===\n")
print(non_front_summary)

# Remove non-front observers from double platform configurations
df_dedup <- df_step2 |>
  filter(!(side_is_double & is_secondary)) |>
  select(-obs_side, -side_is_double, -is_secondary)

# Calculate removals from Step 2
n_removed_step2 <- nrow(df_step1) - nrow(df_dedup)
cat("\nStep 2 - Remove secondary observers from double platforms:\n")
cat("  Records removed:", n_removed_step2, "\n")
cat("  Records remaining:", nrow(df_dedup), "\n")

# Store final deduplicated dataset (no Step 3 needed)
df_dedup_final <- df_dedup

cat("\n=== TOTAL DEDUPLICATION SUMMARY ===\n")
cat("Starting records:   ", nrow(df), "\n")
cat("After Step 1:       ", nrow(df_step1), " (removed:", n_removed_step1, ")\n")
cat("After Step 2:       ", nrow(df_dedup_final), " (removed:", n_removed_step2, ")\n")
cat("Total removed:      ", nrow(df) - nrow(df_dedup_final), 
    " (", round(100 * (nrow(df) - nrow(df_dedup_final)) / nrow(df), 1), "%)\n")

# Summary statistics -----------------------------------------------------------

## Raw data summary
species_all <- df |>
  group_by(Sp_code, Sp_name_clean) |>
  summarise(
    n_sightings = n_distinct(Observation_id),
    n_animals   = sum(Nb_tot, na.rm = TRUE),
    .groups = "drop"
  )

## Deduplicated data summary (final)
species_dedup <- df_dedup_final |>
  group_by(Sp_code, Sp_name_clean) |>
  summarise(
    n_sightings = n_distinct(Observation_id),
    n_animals   = sum(Nb_tot, na.rm = TRUE),
    .groups = "drop"
  )

## Comparison table
comparison <- species_all |>
  rename(sightings_raw = n_sightings, animals_raw = n_animals) |>
  left_join(
    species_dedup |> 
      rename(sightings_dedup = n_sightings, animals_dedup = n_animals),
    by = c("Sp_code", "Sp_name_clean")
  ) |>
  mutate(
    pct_sightings_retained = round(100 * sightings_dedup / sightings_raw, 1),
    pct_animals_retained = round(100 * animals_dedup / animals_raw, 1)
  ) |>
  arrange(desc(sightings_raw))

cat("\n")
print(comparison)


# Create version without unknown species ---------------------------------------
# Filter out all species codes starting with "Un" (unknown species)
df_dedup_known <- df_dedup_final |>
  filter(!grepl("^Un", Sp_code))

cat("\n=== KNOWN SPECIES ONLY ===\n")
cat("Records with unknown species removed:", nrow(df_dedup_final) - nrow(df_dedup_known), "\n")
cat("Records remaining (known species only):", nrow(df_dedup_known), "\n")

## Known species summary
species_known <- df_dedup_known |>
  group_by(Sp_code, Sp_name_clean) |>
  summarise(
    n_sightings = n_distinct(Observation_id),
    n_animals   = sum(Nb_tot, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_sightings))

cat("\nKnown species summary:\n")
print(species_known)

# Spatial visualization --------------------------------------------------------

    ## Helper function to create sf objects
    make_sf <- function(dat) {
      dat |>
        filter(!is.na(Long_sight), !is.na(Lat_sight)) |>
        st_as_sf(coords = c("Long_sight", "Lat_sight"), crs = 4326)
    }

    # Build sf objects for mapping
    g_all   <- make_sf(df_dedup_final)   # all marine mammals, post-dedup
    g_dedup <- make_sf(df_dedup_known)   # known species only, post-dedup

    # Side-by-side maps with unified legend----

    # Get all unique species from both datasets for consistent color mapping
    all_species <- unique(c(g_all$Sp_code, g_dedup$Sp_code))
    
    # Create consistent color scale
    n_species <- length(all_species)
    species_colors <- setNames(
      scales::hue_pal()(n_species),
      all_species
    )
    
    p_all <- ggplot(g_all) +
      geom_sf(size = 1, alpha = 0.5, aes(col = Sp_code)) +
      scale_color_manual(
        name = "Species",
        values = species_colors,
        limits = all_species,
        drop = FALSE
      ) +
      ggtitle(paste0("Deduplicated MM Observations (n = ", nrow(g_all), ")")) +
      theme_minimal() +
      theme(legend.position = "right")
    
    p_dedup <- ggplot(g_dedup) +
      geom_sf(size = 1, alpha = 0.5, aes(col = Sp_code)) +
      scale_color_manual(
        name = "Species",
        values = species_colors,
        limits = all_species,
        drop = FALSE
      ) +
      ggtitle(paste0("Deduplicated species level (n = ", nrow(g_dedup), ")")) +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Combine with unified legend
    p_combined <- (p_all / p_dedup) 
    
    # Display
    print(p_combined)
    
	    # Save as PNG
	    dir.create("output/figs/aerial_deduplication", showWarnings = FALSE, recursive = TRUE)
	    ggsave(
	      filename = "output/figs/aerial_deduplication/deduplication_comparison_map.png",
	      plot = p_combined,
      width = 10,
      height = 12,
      dpi = 300,
      bg = "white"
    )
    
	    cat("Map saved as 'output/figs/aerial_deduplication/deduplication_comparison_map.png'\n")
    
    

# Export deduplicated data (optional) ------------------------------------------
# Uncomment to save:
dir.create(here("input/processed_data/aerial"), showWarnings = FALSE, recursive = TRUE)
write_csv(df_dedup_final, here("input/processed_data/aerial/deduplicated_sightings.csv"))
