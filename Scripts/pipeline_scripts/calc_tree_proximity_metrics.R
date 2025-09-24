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


calc_tree_proximities <- function(tree_crowns, building_fprints) {
  # Validate inputs
  if (!inherits(tree_crowns, "sf") || !inherits(building_fprints, "sf")) {
    stop("Both 'tree_crowns' and 'building_fprints' must be sf objects.")
  }
  if (is.null(sf::st_crs(tree_crowns)) || is.null(sf::st_crs(building_fprints)) ||
      sf::st_crs(tree_crowns) != sf::st_crs(building_fprints)) {
    building_fprints <- st_transform(building_fprints, st_crs(tree_crowns))
    stop("CRS of 'tree_crowns' and 'building_fprints' must be defined and identical.")
  }
  
  # Set up progress bar
  total_steps <- 11
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  step <- 0
  
  message("1/11: Calculating tree centroids...")
  tree_crowns$centroid <- sf::st_centroid(tree_crowns$geometry)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("2/11: Finding nearest building for each tree...")
  nearest <- sf::st_nearest_feature(tree_crowns, building_fprints)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("3/11: Calculating distance to closest building...")
  tree_crowns$dist_to_building <- sf::st_distance(
    tree_crowns,
    building_fprints[nearest, ],
    by_element = TRUE
  )
  setTxtProgressBar(pb, step <- step + 1)
  
  message("4/11: Calculating building centroids...")
  building_fprints$centroid <- sf::st_centroid(building_fprints$geometry)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("5/11: Calculating bearings to closest building...")
  build_cents <- sf::st_transform(building_fprints[nearest, ]$centroid, crs = 4326)
  centroids_ll <- sf::st_transform(tree_crowns$centroid, crs = 4326)
  orig_coords    <- sf::st_coordinates(centroids_ll)
  nearest_coords <- sf::st_coordinates(build_cents)
  valid <- stats::complete.cases(orig_coords, nearest_coords)
  bearings <- rep(NA_real_, nrow(orig_coords))
  bearings[valid] <- geosphere::bearingRhumb(orig_coords[valid, ], nearest_coords[valid, ])
  tree_crowns$bearing_to_closest_building <- bearings
  setTxtProgressBar(pb, step <- step + 1)
  
  message("6/11: Finding nearest tree neighbors...")
  nn <- nngeo::st_nn(tree_crowns$centroid, tree_crowns$centroid, k = 2, maxdist = 100)
  nearest_tree <- sapply(nn, function(idxs) if (length(idxs) >= 2) idxs[2] else NA_integer_)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("7/11: Calculating distance to nearest tree...")
  tree_crowns$dist2n_tree <- sf::st_distance(
    tree_crowns,
    tree_crowns[nearest_tree, ],
    by_element = TRUE
  )
  setTxtProgressBar(pb, step <- step + 1)
  
  message("8/11: Calculating bearings to nearest tree...")
  nearest_tree_cents <- sf::st_transform(tree_crowns[nearest_tree, ]$centroid, crs = 4326)
  centroids_ll <- sf::st_transform(tree_crowns$centroid, crs = 4326)
  orig_coords    <- sf::st_coordinates(centroids_ll)
  nearest_coords <- sf::st_coordinates(nearest_tree_cents)
  valid <- stats::complete.cases(orig_coords, nearest_coords)
  bearings <- rep(NA_real_, nrow(orig_coords))
  bearings[valid] <- geosphere::bearingRhumb(orig_coords[valid, ], nearest_coords[valid, ])
  tree_crowns$bearing_to_nearest_tree <- bearings
  setTxtProgressBar(pb, step <- step + 1)
  
  message("9/11: Counting trees within buffers...")
  tree_crowns$num_trees_2m  <- lengths(sf::st_intersects(sf::st_buffer(tree_crowns, 2),  tree_crowns)) - 1
  tree_crowns$num_trees_5m  <- lengths(sf::st_intersects(sf::st_buffer(tree_crowns, 5),  tree_crowns)) - 1
  tree_crowns$num_trees_10m <- lengths(sf::st_intersects(sf::st_buffer(tree_crowns, 10), tree_crowns)) - 1
  setTxtProgressBar(pb, step <- step + 1)
  
  message("10/11: Calculating tree canopy areas in buffers...")
  tcf <- terra::vect(tree_crowns)
  tu  <- terra::aggregate(tcf, cores = 4)
  # 2m buffer
  tb2 <- terra::buffer(tcf, width = 2)
  iv2 <- terra::intersect(tb2, tu)
  area2 <- terra::expanse(iv2)
  area_tree_df_2m  <- data.frame(fcl_mdn = iv2$fcl_mdn, total_tree_2m = area2)
  # 5m buffer
  tb5 <- terra::buffer(tcf, width = 5)
  iv5 <- terra::intersect(tb5, tu)
  area5 <- terra::expanse(iv5)
  area_tree_df_5m  <- data.frame(fcl_mdn = iv5$fcl_mdn, total_tree_5m = area5)
  # 10m buffer
  tb10 <- terra::buffer(tcf, width = 10)
  iv10 <- terra::intersect(tb10, tu)
  area10 <- terra::expanse(iv10)
  area_tree_df_10m <- data.frame(fcl_mdn = iv10$fcl_mdn, total_trees_10m = area10)
  # join and compute extents
  tree_crowns$canopy_area <- sf::st_area(tree_crowns$geometry)
  area_tree_df <- left_join(as.data.frame(tree_crowns), area_tree_df_2m, by = "fcl_mdn") %>%
    dplyr::left_join(area_tree_df_5m,  by = "fcl_mdn") %>% 
    dplyr::left_join(area_tree_df_10m, by = "fcl_mdn") %>%
    dplyr::mutate(
      area_ext_tree_2  = round(total_tree_2m  - as.numeric(canopy_area), 0),
      area_ext_tree_5  = round(total_tree_5m  - as.numeric(canopy_area), 0),
      area_ext_tree_10 = round(total_trees_10m - as.numeric(canopy_area), 0)
    ) %>%
    dplyr::select(fcl_mdn, area_ext_tree_2, area_ext_tree_5, area_ext_tree_10)
  setTxtProgressBar(pb, step <- step + 1)
  
  message("11/11: Calculating building areas in buffers...")
  fpv <- terra::vect(building_fprints)
  fu  <- terra::aggregate(fpv, cores = 4)
  tb2 <- terra::buffer(tcf, width = 2); iv  <- terra::intersect(tb2, fu); areas <- terra::expanse(iv)
  area_build_df_2m  <- data.frame(fcl_mdn = iv$fcl_mdn, total_build_2m = areas)
  tb5 <- terra::buffer(tcf, width = 5); iv  <- terra::intersect(tb5, fu); areas <- terra::expanse(iv)
  area_build_df_5m  <- data.frame(fcl_mdn = iv$fcl_mdn, total_build_5m = areas)
  tb10 <- terra::buffer(tcf, width = 10); iv <- terra::intersect(tb10, fu); areas <- terra::expanse(iv)
  area_build_df_10m <- data.frame(fcl_mdn = iv$fcl_mdn, total_build_10m = areas)
  tree_crowns_df <- as.data.frame(tree_crowns)
  area_build_df <- tree_crowns_df %>% 
    dplyr::left_join(area_build_df_2m, by = "fcl_mdn") %>%
    dplyr::left_join(area_build_df_5m,  by = "fcl_mdn") %>%
    dplyr::left_join(area_build_df_10m, by = "fcl_mdn") %>% 
    dplyr::select(fcl_mdn, total_build_2m, total_build_5m, total_build_10m)
  tree_crowns2 <- merge(tree_crowns, area_tree_df,  by = "fcl_mdn")
  tree_crowns3 <- merge(tree_crowns2, area_build_df, by = "fcl_mdn")
  setTxtProgressBar(pb, step <- step + 1)
  
  close(pb)
  message("Done calculating proximities.")
  return(tree_crowns3)
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
result <- calc_tree_proximities(tree_sf, build_sf)

# write
sf::st_write(result, opt$output, append = FALSE)

  
  