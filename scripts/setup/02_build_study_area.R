# ==============================================================================
# Laura Joan Feyrer
# Date Updated: 2026-05-26
# Script: setup/02_build_study_area.R
# Description: Builds the Eastern Scotian Shelf study area polygon by
#              intersecting NAFO 4V/4W divisions with the OSW regional
#              assessment area, keeping the largest polygon, filling holes,
#              simplifying, and saving to shapefiles/studyArea/.
#              Run once per project, or whenever boundary inputs change.
# Changes from previous version:
#   - Fixed OSW_studyarea load: now uses spatial_path("OSW_studyarea") via
#     load_spatial_layer() instead of bare undefined variable.
#   - Fixed map section: replaced undefined cont_UTM/landUTM with
#     load_core_spatial_layers() calls.
#   - Added ggplot2 to package loading.
# ==============================================================================

# Load packages ----------------------------------------------------------------
pacman::p_load(sf, dplyr, here, terra, smoothr, ggplot2)
suppressWarnings(source(here::here("scripts/00_load_helpers.R")))

# ---- 1. Read in existing shapefiles ------------------------------------------
# All paths are managed in scripts/helpers/helper_spatial.R
# To point to different files, edit the SPATIAL_PATHS list there.

nafo          <- st_read(spatial_path("nafo"),    quiet = TRUE)
eez           <- st_read(spatial_path("eez"),     quiet = TRUE)
land          <- st_read(spatial_path("land"),    quiet = TRUE)
WEA           <- st_read(spatial_path("wea"),     quiet = TRUE)
OSW_studyarea <- st_read(spatial_path("OSW_studyarea"), quiet = TRUE)

plot(st_geometry(OSW_studyarea))


# Pick a common CRS (WGS84 here, use UTM20 or whatever you are mapping in)
target_crs <- 32620

nafo <- st_transform(nafo, target_crs)
eez  <- st_transform(eez, target_crs)
OSW_studyarea  <- st_transform(OSW_studyarea, target_crs)
WEA  <- st_transform(WEA, target_crs)



# --2. Filter NAFO to 3V and 3W, union to one----------------------

# Adjust DIV / DIVISION to match your attribute name
nafo_ess <- nafo |>
  filter(Division %in% c("4V", "4W"))

nafo_ess_union <- nafo_ess |>
  st_make_valid() |>
  st_union() |>
  st_as_sf() |>
  mutate(id = 1)  # simple id so it has an attribute

plot(st_geometry(nafo_ess_union))

# --- 3. Clip to Regional assessment area to stay inside Canada------------------------------------------


ess_nafo_RA <- st_intersection(nafo_ess_union, OSW_studyarea)

 plot(st_geometry(ess_nafo_RA))


# # ---- 4. Remove deep water as not economically feasible for development---------------------------------------------
# 
#  # Create a polygon mask for depth ≤ 400 m
# 
#  # Note: depth values are negative, so threshold is -400
#  shallow_mask <- bathy <= -400   # returns TRUE for deep water; we invert
#  
#  # Create shallow area raster (depth < 400 m)
#  shallow_area <- bathy >= -400   # TRUE where water is shallower than 400m
#  
#  # Convert raster mask to polygons
#  shallow_poly <- as.polygons(shallow_area, dissolve = TRUE) |> 
#    st_as_sf() |> 
#    st_make_valid()
#  
#  # Keep only TRUE pixels
#  shallow_poly <- shallow_poly |> filter(gebco_2020 == 1)
#  
#  # Clip ESS study area to <400 m --------------------------------
#  
#  ess_shallow <- st_intersection(ess_nafo_RA, shallow_poly) |> 
#    st_make_valid() |> 
#    mutate(depth_filter = "<400m")
#  
#  plot(st_geometry(ess_shallow))
#  
 
 # Clip to main polygon (largest area)-----
 
 # 2. Explode multipolygons to single polygons
 ess_parts_utm <- st_cast(ess_nafo_RA, "MULTIPOLYGON")
 ess_parts_utm <- st_cast(ess_parts_utm, "POLYGON")
 
 # 3. Compute area in km2 OUTSIDE any pipe
 #    Use sf::st_area explicitly so we don't hit some other st_area method
 areas <- sf::st_area(ess_parts_utm)
 ess_parts_utm$area_km2 <- as.numeric(areas) / 1e6
 
 # 4. Keep only the largest polygon
 ess_parts_utm <- ess_parts_utm[order(ess_parts_utm$area_km2, decreasing = TRUE), ]
 
 # First row is largest
 ess_main_utm <- ess_parts_utm[1, ]
 
 # Union just in case
 ess_main_utm <- st_union(ess_main_utm)
 ess_main_utm <- st_as_sf(ess_main_utm)
 
 plot(st_geometry(ess_main_utm))
 
 
 # 4. Remove interior holes
 #    threshold is in square metres (CRS is UTM). Set it huge so all holes are filled.
 ess_main_utm <- smoothr::fill_holes(
   ess_main_utm,
   threshold = 1e12   # 1,000,000 km^2, safely bigger than your shelf polygon
 )
 plot(st_geometry(ess_main_utm))
 
 # 5. Simplify
 ess_simple_utm <- st_simplify(
   ess_main_utm,
   dTolerance = 2000,
   preserveTopology = TRUE
 )
 plot(st_geometry(ess_simple_utm))
 
 # # 6. Optional smoothing
 # ess_smooth_utm <- smoothr::smooth(
 #   ess_simple_utm,
 #   method = "ksmooth",
 #   smoothness = 2
 # )
 # plot(st_geometry(ess_smooth_utm))
 
 # 7. Save
 ess_study_final <- ess_simple_utm
 ess_study_final$study_area <- "ess_simple_utm"
 
 st_write(
   ess_study_final,
   spatial_path("study_area"),
   delete_dsn = TRUE
 )
 
 
 # Make a map to visualize extent ---------------------------------------------
 # Load land and bathy contours via the helper (same layers used in all scripts)
 base_layers <- load_core_spatial_layers()
 land_map    <- base_layers$land
 cont_map    <- base_layers$contours   # NULL if shapefile absent -- handled below

 bbox <- st_bbox(ess_simple_utm)

 xmin <- bbox["xmin"]
 xmax <- bbox["xmax"]
 ymin <- bbox["ymin"]
 ymax <- bbox["ymax"]

 xlims <- c(xmin - 100000, xmax)
 ylims <- c(ymin, ymax)

 # For the legend
 OSW_studyarea$type   <- "Offshore wind regional assessment area"
 ess_study_final$type <- "Eastern Scotian Shelf study area"
 WEA$type             <- "Wind Energy AOIs"

 gg_map <- ggplot() +
    theme_bw() +

    # bathy contours (skipped gracefully if file is absent)
    {if (!is.null(cont_map))
       geom_sf(
         data = dplyr::filter(cont_map, level %in% c(-200, -400, -1000, -2500, -3200)),
         col = MAP_LAYER_STYLE$bathy_color, linewidth = 0.2
       )
    } +

    # land
    geom_sf(data = land_map, color = MAP_LAYER_STYLE$land_color,
            fill = MAP_LAYER_STYLE$land_fill) +
    
    # OSW study area
    geom_sf(
       data = OSW_studyarea,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    # ESS < 400 m
    geom_sf(
       data = ess_study_final,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    # Designated wind energy areas
    geom_sf(
       data = WEA,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    scale_fill_manual(
       name   = "",
       values = c(
          "Offshore wind regional assessment area" = "grey90",
          "Eastern Scotian Shelf study area"              = "lightblue",
          "Wind Energy AOIs"        = MAP_LAYER_STYLE$wea_color
       )
    ) +
    scale_color_manual(
       name   = "",
       values = c(
          "Offshore wind regional assessment area" = "black",
          "Eastern Scotian Shelf study area"              = "blue",
          "Wind Energy AOIs"        = MAP_LAYER_STYLE$wea_color
       )
    ) +
    
    coord_sf(
       lims_method = "orthogonal",
       xlim = xlims,
       ylim = ylims,
       crs  = UTM20,
       expand = TRUE
    ) +
    
    labs(x = "", y = "") +
    
    theme(
       legend.position      = c(0.98, 0.02),  # bottom right inside plot
       legend.justification = c("right", "bottom"),
       legend.background    = element_rect(fill = "white", color = "grey80"),
       legend.key.size      = unit(0.6, "cm"),
       legend.title = element_blank()
    )
 
 
 #visualize & save plot
 
 gg_map
 ggsave("output/figs/proposed_StudyAreamap.png", gg_map, dpi = 300)
 
 
