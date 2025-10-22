#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(terra)
  library(geosphere)
  library(nngeo)
  library(dplyr)
  # …etc…
})

calc_building_proximities <- function(building_fprints, tree_crowns) {
  # Validate inputs
  if (!inherits(tree_crowns, "sf") || !inherits(building_fprints, "sf")) {
    stop("Both 'tree_crowns' and 'building_fprints' must be sf objects.")
  }
  if (is.null(sf::st_crs(tree_crowns)) || is.null(sf::st_crs(building_fprints)) ||
      sf::st_crs(tree_crowns) != sf::st_crs(building_fprints)) {
    building_fprints <- st_transform(building_fprints, st_crs(tree_crowns))
    #stop("CRS of 'tree_crowns' and 'building_fprints' must be defined and identical.")
  }
  
  # Set up progress bar
  total_steps <- 10
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0
  
  message("1/10: Calculating centroids...")
  tree_crowns$centroid <- sf::st_centroid(tree_crowns$geometry)
  building_fprints$centroid <- sf::st_centroid(building_fprints$geometry)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("2/10: Finding nearest building to focal building...")
  nn <- nngeo::st_nn(building_fprints$centroid, building_fprints$centroid, k = 2, maxdist = 100)
  nearest_build <- sapply(nn, function(idxs) if (length(idxs) >= 2) idxs[2] else NA_integer_)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("3/10: Calculating distance to closest building...")
  building_fprints$dist_to_nearest_building <- sf::st_distance(
    building_fprints,
    building_fprints[nearest_build, ],
    by_element = TRUE
  )
  setTxtProgressBar(pb, step <- step + 1)
  
  message("4/10: Calculating bearings to closest building...")
  build_cents <- sf::st_transform(building_fprints[nearest_build, ]$centroid, crs = 4326)
  centroids_ll <- sf::st_transform(building_fprints$centroid, crs = 4326)
  orig_coords    <- sf::st_coordinates(centroids_ll)
  nearest_coords <- sf::st_coordinates(build_cents)
  valid <- stats::complete.cases(orig_coords, nearest_coords)
  bearings <- rep(NA_real_, nrow(orig_coords))
  bearings[valid] <- geosphere::bearingRhumb(orig_coords[valid, ], nearest_coords[valid, ])
  building_fprints$bearing_to_closest_building <- bearings
  setTxtProgressBar(pb, step <- step + 1)
  
  message("5/10: Finding nearest tree neighbors...")
  nearest_tree <- sf::st_nearest_feature(building_fprints, tree_crowns)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("6/10: Calculating distance to nearest tree...")
  building_fprints$dist_to_nearest_tree <- sf::st_distance(
    building_fprints,
    tree_crowns[nearest_tree, ],
    by_element = TRUE
  )
  setTxtProgressBar(pb, step <- step + 1)
  
  message("7/10: Calculating bearings to nearest tree...")
  nearest_tree_cents <- sf::st_transform(tree_crowns[nearest_tree, ]$centroid, crs = 4326)
  centroids_ll <- sf::st_transform(building_fprints$centroid, crs = 4326)
  orig_coords    <- sf::st_coordinates(centroids_ll)
  nearest_coords <- sf::st_coordinates(nearest_tree_cents)
  valid <- stats::complete.cases(orig_coords, nearest_coords)
  bearings <- rep(NA_real_, nrow(orig_coords))
  bearings[valid] <- geosphere::bearingRhumb(orig_coords[valid, ], nearest_coords[valid, ])
  building_fprints$bearing_to_nearest_tree <- bearings
  setTxtProgressBar(pb, step <- step + 1)
  
  message("8/10: Counting trees within buffers...")
  building_fprints$num_trees_2m  <- building_fprints %>% st_buffer(2) %>% 
    st_intersects(tree_crowns) %>% 
    lengths()
  building_fprints$num_trees_5m  <- building_fprints %>% st_buffer(5) %>% 
    st_intersects(tree_crowns) %>% 
    lengths()
  building_fprints$num_trees_10m <- building_fprints %>% st_buffer(10) %>% 
    st_intersects(tree_crowns) %>% 
    lengths()
  setTxtProgressBar(pb, step <- step + 1)
  message("9/10: Calculating tree canopy areas in buffers...")
  tcv <- terra::vect(tree_crowns)
  fpv <- vect(building_fprints)
  cu   <- aggregate(tcv, cores = 4)          # union
  fb2  <- buffer(fpv, width = 2)  # 2 m buffer
  iv2 <- terra::intersect(fb2, cu)
  area2 <- terra::expanse(iv2)
  area_tree_df_2m  <- data.frame(UID = iv2$UID, area_tree_2m = area2)
  # 5m buffer
  fb5 <- terra::buffer(fpv, width = 5)
  iv5 <- terra::intersect(fb5, cu)
  area5 <- terra::expanse(iv5)
  area_tree_df_5m  <- data.frame(UID = iv5$UID, area_tree_5m = area5)
  # 10m buffer
  fb10 <- terra::buffer(fpv, width = 10)
  iv10 <- terra::intersect(fb10, cu)
  area10 <- terra::expanse(iv10)
  area_tree_df_10m <- data.frame(UID = iv10$UID, area_tree_10m = area10)
  # join and compute extents
  fprints_df <- as.data.frame(building_fprints)
  area_tree_df <- fprints_df %>% 
    dplyr::left_join(area_tree_df_2m, by = "UID") %>%
    dplyr::left_join(area_tree_df_5m,  by = "UID") %>% 
    dplyr::left_join(area_tree_df_10m, by = "UID") %>% 
    select(UID, area_tree_2m, area_tree_5m, area_tree_10m)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("10/10: Calculating building areas in buffers...")
  building_fprints$footprint_area <- st_area(building_fprints$geometry)
  
  # Helper: calculate building area within any buffer distance
  calc_building_area <- function(building_fprints, dist) {
    fpv <- terra::vect(building_fprints)
    fu  <- terra::aggregate(fpv, cores = 4)
    fb  <- terra::buffer(fpv, width = dist)
    iv  <- terra::intersect(fb, fu)
    areas <- terra::expanse(iv)
    col_name <- paste0("total_build_", dist, "m")
    df <- data.frame(UID = iv$UID)
    df[[col_name]] <- areas
    return(df)
  }
  
  area_build_df_2m <- calc_building_area(building_fprints, 2)
  
  area_build_df_10m <- calc_building_area(building_fprints, 10)
  
  area_build_df_25m <- calc_building_area(building_fprints, 25)
  
  area_build_df_50m <- calc_building_area(building_fprints, 50)
  
  area_build_df_100m <- calc_building_area(building_fprints, 100)
  
    area_build_df <- left_join(as.data.frame(building_fprints), area_build_df_2m, by = "UID") %>%
    dplyr::left_join(area_build_df_10m, by = "UID") %>% 
    dplyr::left_join(area_build_df_25m, by = "UID") %>% 
    dplyr::left_join(area_build_df_50m, by = "UID") %>%
    dplyr::left_join(area_build_df_100m, by = "UID") %>% 
      mutate(area_ext_build_2 = round(total_build_2m - as.numeric(footprint_area), digits = 0), 
             area_ext_build_10 = round(total_build_10m - as.numeric(footprint_area), digits = 0),
             area_ext_build_25 = round(total_build_25m - as.numeric(footprint_area), digits = 0), 
             area_ext_build_50 = round(total_build_50m - as.numeric(footprint_area), digits = 0), 
             area_ext_build_100 = round(total_build_100m - as.numeric(footprint_area), digits = 0)) %>% 
      select(UID, area_ext_build_2, area_ext_build_10, area_ext_build_25, area_ext_build_50, area_ext_build_100)
    
  building_fprints <- merge(building_fprints, area_tree_df,  by = "UID")
  building_fprints <- merge( building_fprints, area_build_df, by = "UID")
  building_fprints <- building_fprints %>% rename(dist2n_build = dist_to_nearest_building, dist2n_tree = dist_to_nearest_tree)
  setTxtProgressBar(pb, step <- step + 1)
  
  close(pb)
  message("Done calculating proximities.")
  return(building_fprints)
}

# define CLI options
option_list <- list(
  make_option(c("-t","--trees"),   type="character", help="Path to tree crowns (sf‐readable)"),
  make_option(c("-b","--buildings"),type="character", help="Path to building footprints"),
  make_option(c("-o","--output"),   type="character", default="proximities.geojson",
              help="Output path [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# validate
if (is.null(opt$trees) || is.null(opt$buildings)) {
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}

# read
tree_sf  <- sf::st_read(opt$trees,   quiet=TRUE)
build_sf <- sf::st_read(opt$buildings,quiet=TRUE)

# run your function
result <- calc_building_proximities(building_fprints = build_sf, tree_crowns = tree_sf)



# write
sf::st_write(result, opt$output, append = FALSE)

