#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(terra)
  library(dplyr)
  library(spatialEco)
  # …etc…
})


extract_kernel_density <- function(tree_crowns, fprints) {
  # 1. Compute centroids as standalone sf objects
  tree_cents  <- st_centroid(tree_crowns)
  build_cents <- st_centroid(fprints)
  
  # 2. Kernel‐density of tree centroids
  tree_dens_rast <- sf.kde(tree_cents, res = 30, bw = 200)
  
  # 3. Extract tree density at crown polygons
  crown_tree_density <- terra::extract(
    tree_dens_rast, tree_crowns,
    fun  = mean,
    bind = TRUE
  ) %>%
    st_as_sf() %>%
    select(fcl_mdn, tree_dens = z)
  
  # 5. Extract tree density at footprints
  footprint_tree_density <- terra::extract(
    tree_dens_rast, fprints,
    fun  = mean,
    bind = TRUE
  ) %>%
    st_as_sf() %>%
    select(UID, tree_dens = z)
  
  # 6. Kernel‐density of building centroids
  build_dens_rast <- sf.kde(build_cents, res = 30, bw = 200)
  
  # 7. Extract building density at crowns
  crown_build_density <- terra::extract(
    build_dens_rast, tree_crowns,
    fun  = mean,
    bind = TRUE
  ) %>%
    st_as_sf() %>%
    select(fcl_mdn, build_dens = z)
  
  # 8. Extract building density at footprints
  footprint_build_density <- terra::extract(
    build_dens_rast, fprints,
    fun  = mean,
    bind = TRUE
  ) %>%
    st_as_sf() %>%
    select(UID, build_dens = z)
  
  # 9. Return a named list
  return(list(
    crown_tree_density     = crown_tree_density,
    footprint_tree_density = footprint_tree_density,
    crown_build_density    = crown_build_density,
    footprint_build_density= footprint_build_density
  ))
}

# define CLI options
option_list <- list(
  make_option(c("-t","--trees"),   type="character", help="Path to tree crowns (sf‐readable)"),
  make_option(c("-b","--buildings"),type="character", help="Path to building footprints"),
  make_option(c("-B","--buildings_output"),   type="character", default="proximities.geojson"),
  make_option(c("-T","--tree_output"),   type="character", default="proximities.geojson",
              help="Output path [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# validate
if (is.null(opt$trees) || is.null(opt$buildings)) {
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}

tree_crowns <- st_read(opt$trees, quiet = TRUE)

fprints <- st_read(opt$buildings, quiet = TRUE)

message("extracting kernel densities")
kds <- extract_kernel_density(tree_crowns = tree_crowns, fprints = fprints)


tree_crowns_dens <- as.data.frame(kds$crown_tree_density) %>% left_join(as.data.frame(kds$crown_build_density), by = "fcl_mdn") %>% 
  select(fcl_mdn, tree_dens, build_dens)

names(tree_crowns_dens)

names(tree_crowns)

message("merging with original tree crowns")
tree_crowns_dens <- base::merge(tree_crowns, tree_crowns_dens, by = "fcl_mdn")

fprints_dens <- as.data.frame(kds$footprint_tree_density) %>% left_join(as.data.frame(kds$footprint_build_density), by = "UID") %>% 
  select(tree_dens, build_dens, "UID")

fprints_dens <- merge(fprints, fprints_dens, by = "UID")

st_write(tree_crowns_dens, opt$tree_output, append=FALSE)

st_write(fprints_dens, opt$buildings_output, append=FALSE)


