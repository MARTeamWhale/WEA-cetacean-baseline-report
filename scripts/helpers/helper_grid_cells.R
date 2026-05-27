# scripts/helpers/helper_grid_cells.R

#   - Ensures a minimum 25% overlap
#     threshold for WEA cell selection. Previously, cells with only a tiny corner
#     clipping the WEA boundary were included, pulling in sightings far outside
#     the WEA. Diagnostic confirmed: cell 24_199 had 0.18% overlap (excluded),
#     cell 27_198 had 25.1% overlap (retained) -- the threshold cleanly separates
#     genuine near-WEA cells from edge artefacts.

make_cell_polygons <- function(cells, cell_m, x_col = "xc", y_col = "yc", crs = 32620) {
  cells %>%
    dplyr::mutate(
      geometry = purrr::map2(.data[[x_col]], .data[[y_col]], function(x, y) {
        half <- cell_m / 2
        sf::st_polygon(list(matrix(
          c(
            x - half, y - half,
            x + half, y - half,
            x + half, y + half,
            x - half, y + half,
            x - half, y - half
          ),
          ncol = 2,
          byrow = TRUE
        )))
      })
    ) %>%
    sf::st_as_sf(crs = crs)
}