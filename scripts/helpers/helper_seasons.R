# scripts/helpers/helper_seasons.R
#
# Single shared definition of season boundaries (Northern Hemisphere), used
# by every script in the pipeline that classifies records by season:
#   04_summarize_wsdb_effort.R      - dt$season (WSDB sightings)
#   06_summarize_wsdb_seasonality.R - reuses dt$season from 04
#   07_summarize_pam_detections.R   - whale_data$Season (PAM detections)
#
# Change SEASON_MONTHS here and every seasonal figure/table in the report
# follows automatically - there should be no other place in the pipeline
# that hardcodes which calendar months belong to which season.
#
# Requires SEASON_LEVELS (scripts/helpers/helper_palettes.R) to already be
# loaded, so this file must be sourced after helper_palettes.R (see
# scripts/00_load_helpers.R).

# Winter: Dec, Jan, Feb | Spring: Mar, Apr, May
# Summer: Jun, Jul, Aug | Fall:   Sep, Oct, Nov
SEASON_MONTHS <- list(
  Winter = c(12, 1, 2),
  Spring = c(3, 4, 5),
  Summer = c(6, 7, 8),
  Fall   = c(9, 10, 11)
)

# Guard against SEASON_MONTHS and SEASON_LEVELS drifting apart, and against
# a typo leaving a month unassigned or double-assigned.
stopifnot(
  "SEASON_MONTHS names must match SEASON_LEVELS" =
    setequal(names(SEASON_MONTHS), SEASON_LEVELS),
  "SEASON_MONTHS must cover each of the 12 months exactly once" =
    setequal(unlist(SEASON_MONTHS, use.names = FALSE), 1:12) &&
    length(unlist(SEASON_MONTHS, use.names = FALSE)) == 12
)

#' Classify records into seasons by calendar month
#'
#' @param data       A data frame / tibble to classify.
#' @param date_col   Unquoted date column to derive month from (only used
#'                    when `month_col` is not already present in `data`).
#' @param season_col Name of the output column (string). Defaults to
#'                    "season"; pass e.g. "Season" to match a script's
#'                    existing naming convention.
#' @param month_col  Name of an existing numeric month column to reuse
#'                    instead of recomputing it from `date_col`.
#' @return `data` with `season_col` added as a factor leveled per SEASON_LEVELS.
add_season <- function(data, date_col = date_utc, season_col = "season", month_col = NULL) {

  month_vals <- if (!is.null(month_col) && month_col %in% names(data)) {
    data[[month_col]]
  } else {
    lubridate::month(dplyr::pull(data, {{ date_col }}))
  }

  season_vals <- rep(NA_character_, length(month_vals))
  for (nm in names(SEASON_MONTHS)) {
    season_vals[month_vals %in% SEASON_MONTHS[[nm]]] <- nm
  }

  data[[season_col]] <- factor(season_vals, levels = SEASON_LEVELS)
  data
}
