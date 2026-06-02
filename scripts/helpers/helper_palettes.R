# scripts/helpers/helper_palettes.R
#
# Shared colour constants for figures. Source via scripts/00_load_helpers.R.

SEASON_PALETTE <- c(
  "Winter" = "#2E86AB",
  "Spring" = "#06A77D",
  "Summer" = "#fd8d3c",
  "Fall"   = "#A23B72"
)

SEASON_LEVELS <- names(SEASON_PALETTE)

MAP_LAYER_STYLE <- list(
  land_fill         = "grey60",
  land_color        = NA,
  bathy_color       = "grey75",
  wea_color         = "#FF6B35",
  study_area_color  = "black",
  kde_contour_color = "lightblue"
)

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

WHALE_PALETTE <- c(
  "Blue Whale"                   = "#3288BD",
  "Fin Whale"                    = "#276B95",
  "Sei Whale"                    = "#2A5857",
  "Fin/Sei Whale"                = "#7AB0C0",
  "Humpback Whale"               = "#6DAFB1",
  "Minke Whale"                  = "#9bc4f8",
  "North Atlantic Right Whale"   = "#7C6FB3",
  "Northern Bottlenose Whale"    = "#B07D62",
  "Sowerby's Beaked Whale"       = "#8B6B4E",
  "Cuvier's Beaked Whale"        = "#C9A97A",
  "True's/Gervais' Beaked Whale" = "#7A6651",
  "Common Dolphin"               = "#C73E5A",
  "Atlantic Bottlenose Dolphin"  = "#A31D3F",
  "Atlantic White-Sided Dolphin" = "#D94A76",
  "Beluga Whale"                 = "#B23A8F",
  "False Killer Whale"           = "#6F2C8F",
  "Killer Whale"                 = "#3B1F5C",
  "White-Beaked Dolphin"         = "#8E2F6B",
  "Striped Dolphin"              = "#E05A8A",
  "Risso's Dolphin"              = "#7A3E8E",
  "Long-Finned Pilot Whale"      = "#4B2E83",
  "Sperm Whale"                  = "#5B1A4B",
  "Harbour Porpoise"             = "#C04B2B"
)
