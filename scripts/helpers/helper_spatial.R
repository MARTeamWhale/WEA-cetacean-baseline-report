# scripts/helpers/helper_spatial.R
#
# NOTE FOR NEW USERS: All paths in SPATIAL_PATHS are relative to the project
# root. They will resolve correctly as long as you open the project via
# WEA-cetacean-baseline-report.Rproj in RStudio (which sets the working
# directory automatically). Do not run scripts by sourcing them from outside
# the project or the shapefiles will not be found.
#
# To point to a different shapefile, update the relevant entry in SPATIAL_PATHS
# below -- all scripts that use that layer will pick up the change automatically.

SPATIAL_CRS_UTM20 <- 32620
UTM20 <- sf::st_crs(SPATIAL_CRS_UTM20)

SPATIAL_PATHS <- list(
  study_area       = "shapefiles/studyArea/ESS_study_area_simple.shp",
  wea              = "shapefiles/WEA/Designated_WEAs_25_07_29.shp",
  wea_5km_buffer   = "shapefiles/WEA_5km_buffer_combined/WEA_5km_buffer_combined.shp",
  land             = "shapefiles/Canada/NE_10mLand.shp",
  bathy_contours   = "shapefiles/bathy/bathymetry_l_v2.shp",
  bathy_study_area = "shapefiles/bathy/ESS_study_area_lt400m.shp",
  eez              = "shapefiles/EEZ/EEZ_can.shp",
  nafo             = "shapefiles/NAFO/NAFO_Divisions_2021_poly_clipped.shp",
  OSW_studyarea = "shapefiles/OSW_Study_Area/Study_Area_NS.shp"
  
)

spatial_path <- function(name) {
  if (!name %in% names(SPATIAL_PATHS)) {
    stop("Unknown spatial path key: ", name)
  }
  SPATIAL_PATHS[[name]]
}

load_spatial_layer <- function(name,
                               crs = SPATIAL_CRS_UTM20,
                               required = TRUE,
                               quiet = TRUE) {
  path <- spatial_path(name)
  if (!file.exists(path)) {
    if (required) stop("Spatial file not found for ", name, ": ", path)
    return(NULL)
  }
  sf::st_read(path, quiet = quiet) |>
    sf::st_make_valid() |>
    sf::st_transform(crs)
}

load_bathy_contours <- function(levels = c(200, 500,1000, 2500),
                                crs = SPATIAL_CRS_UTM20,
                                required = FALSE) {
  contours <- load_spatial_layer("bathy_contours", crs = crs, required = required)
  if (is.null(contours)) return(NULL)

  if ("DEPTH" %in% names(contours) && !"level" %in% names(contours)) {
    contours <- dplyr::rename(contours, level = DEPTH)
  }
  if (!is.null(levels) && "level" %in% names(contours)) {
    contours <- dplyr::filter(contours, .data$level %in% levels)
  }
  contours
}

load_core_spatial_layers <- function(crs = SPATIAL_CRS_UTM20,
                                     bathy_levels = c(200, 500,1000, 2500)) {
  list(
    study_area = load_spatial_layer("study_area", crs = crs),
    wea        = load_spatial_layer("wea", crs = crs),
    land       = load_spatial_layer("land", crs = crs),
    contours   = load_bathy_contours(levels = bathy_levels, crs = crs, required = FALSE)
  )
}

map_study_bbox <- function(study_area,
                           crs = UTM20,
                           buffer_m = 0) {
  bbox <- study_area |>
    sf::st_transform(crs) |>
    sf::st_bbox()

  list(
    xlim = c(bbox["xmin"] - buffer_m, bbox["xmax"] + buffer_m),
    ylim = c(bbox["ymin"] - buffer_m, bbox["ymax"] + buffer_m)
  )
}

map_base_theme <- function(base_size = 11,
                           show_axis_text = FALSE,
                           legend_position = "right") {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = if (show_axis_text) {
        ggplot2::element_text()
      } else {
        ggplot2::element_blank()
      },
      legend.position = legend_position
    )
}

map_coord <- function(xlim,
                      ylim,
                      crs = UTM20,
                      expand = FALSE) {
  ggplot2::coord_sf(xlim = xlim, ylim = ylim, crs = crs, expand = expand)
}

map_layer_land <- function(land,
                           fill = MAP_LAYER_STYLE$land_fill,
                           color = MAP_LAYER_STYLE$land_color,
                           linewidth = 0.2) {
  ggplot2::geom_sf(data = land, fill = fill, color = color, linewidth = linewidth)
}

map_layer_wea <- function(wea,
                          color = MAP_LAYER_STYLE$wea_color,
                          linewidth = .7, 
                          alpha = .05) {
  ggplot2::geom_sf(data = wea, fill = NA, color = color, linewidth = linewidth,
                   alpha = alpha)
}

map_layer_study_area <- function(study_area,
                                 color = MAP_LAYER_STYLE$study_area_color,
                                 linewidth = 0.7,
                                 linetype = "dashed") {
  ggplot2::geom_sf(data = study_area, fill = NA, color = color,
                   linewidth = linewidth, linetype = linetype)
}

map_layer_bathy_contours <- function(contours,
                                     levels = c(200, 500,1000, 3000),
                                     color = MAP_LAYER_STYLE$bathy_color,
                                     linewidth = 0.45) {
  if (is.null(contours)) {
    return(ggplot2::geom_blank())
  }
  if ("level" %in% names(contours)) {
    contours <- dplyr::filter(contours, .data$level %in% levels)
  }
  ggplot2::geom_sf(data = contours, color = color, linewidth = linewidth)
}

map_layer_kde_contour <- function(contour_sf,
                                  color = MAP_LAYER_STYLE$kde_contour_color,
                                  linewidth = 0.3) {
  ggplot2::geom_sf(data = contour_sf, fill = NA, color = color, linewidth = linewidth)
}

map_annotation_label <- function(x,
                                 y,
                                 label,
                                 vjust,
                                 hjust = -0.05,
                                 size = 2) {
  ggplot2::annotate(
    "label",
    x = x, y = y, label = label,
    hjust = hjust, vjust = vjust,
    size = size, color = "darkblue",
    fill = "white", label.size = 0.15,
    alpha = 0.85
  )
}
