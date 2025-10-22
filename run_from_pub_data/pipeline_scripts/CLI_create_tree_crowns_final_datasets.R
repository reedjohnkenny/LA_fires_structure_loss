#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(terra)
  library(dplyr)
  # …etc…
})




process_tree_crowns <- function(crowns, fprints, perim, study_area, burn_level) {
  if (burn_level == "unburned") {
    message("retain only crowns within 50m of structure")
    crowns <- st_join(crowns, fprints, join = st_is_within_distance, dist = 50, left = FALSE)
    message("retain only crowns outside of fire perimeter")
    crowns <- crowns %>% filter(!duplicated(fcl_mdn)) %>% st_difference(perim)
  } else {
    perim_minus_50 <- st_buffer(perim, -50)
    message("crop crowns to study area")
    crowns <- st_join(crowns, study_area, left = FALSE)
    message("retain only crowns within 50m of structure")
    crowns <- st_join(crowns, fprints, join = st_is_within_distance, dist = 50, left = FALSE)
    message("retain only crowns inside of fire perimeter")
    crowns <- st_join(crowns, perim_minus_50, left = FALSE)
    crowns <- crowns %>% filter(!duplicated(fcl_mdn))
  }
  crowns$centroid <- st_centroid(crowns$geometry)
  crowns$canopy_area <- st_area(crowns)
  return(crowns)
}


# define CLI options
option_list <- list(
  make_option(c("-t","--trees"),   type="character", help="Path to tree crowns (sf‐readable)"),
  make_option(c("-b","--buildings"),   type="character", help="Path to building footprints (sf‐readable)"),
  make_option(c("-f","--fire"),type="character", help="name of fire, must be eaton or palisades"),
  make_option(c("-s","--study_area"),type="character", help="path to study area (sf-readable)"),
  make_option(c("-l","--burn_level"),type="character", help="must be burned or unburned"),
  make_option(c("-p","--burn_perimeter", type="character", help="path_to_burn_perimeters")),
  make_option(c("-o","--output"),   type="character", default="proximities.geojson", help="Output path [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))


# ——— Validate required args ———
required <- c("trees","buildings","fire","study_area","burn_level")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing)) {
  message("Missing arguments: ", paste(missing, collapse=", "))
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}

# ——— Read inputs ———
crowns     <- st_read(opt$trees,      quiet=TRUE) %>% st_transform(32611)
buildings  <- st_read(opt$buildings,  quiet=TRUE) %>% st_transform(32611)
study_area <- st_read(opt$study_area, quiet=TRUE) %>% st_transform(32611)
perim_path <- opt$burn_perimeter

# ——— Load the correct perimeter based on fire name ———
if (tolower(opt$fire) == "eaton") {
  perim <- st_read(
    perim_path) %>%
    filter(poly_Incid == "Eaton") %>%
    st_transform(32611) %>%
    select(poly_Incid, poly_Featu)
} else {
  perim <- st_read(
    perim_path) %>%
    filter(poly_Incid == "PALISADES") %>%
    st_transform(32611) %>%
    select(poly_Incid, poly_Featu)
}

# ——— Process & write ———
crowns_final <- process_tree_crowns(
  crowns, buildings, perim,
  study_area = study_area,
  burn_level = tolower(opt$burn_level)
)

st_write(crowns_final, opt$output, delete_dsn = TRUE)