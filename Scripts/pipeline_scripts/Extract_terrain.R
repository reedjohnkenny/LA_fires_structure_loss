#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(terra)
  library(dplyr)
})

extract_terrain <- function(dtm, fprints) {
  # Compute slope and aspect from the DTM
  slp <- terrain(dtm, v = "slope",  unit = "degrees")
  asp <- terrain(dtm, v = "aspect", unit = "degrees")

  # Stack elevation, slope, aspect
  terrain_stack <- c(dtm, slp, asp)
  names(terrain_stack) <- c("elevation", "slope", "aspect")

  # Extract mean values per polygon
  terrain_vals <- terra::extract(
    terrain_stack, fprints,
    fun  = mean, na.rm = TRUE,
    bind = TRUE
  ) %>%
    st_as_sf()

  return(terrain_vals)
}

# CLI options
option_list <- list(
  make_option(c("-d", "--dtm"),       type = "character", help = "Path to DTM raster"),
  make_option(c("-b", "--buildings"), type = "character", help = "Path to building footprints"),
  make_option(c("-o", "--output"),    type = "character", default = "tmp_data/eaton_burned_fprints_terrain.shp",
              help = "Output path [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dtm) || is.null(opt$buildings)) {
  print_help(OptionParser(option_list = option_list))
  quit(status = 1)
}

dtm     <- rast(opt$dtm)
dtm_m <- dtm * 0.3048
fprints <- st_read(opt$buildings, quiet = TRUE)


message("Extracting terrain variables (elevation, slope, aspect)...")
fprints_terrain <- extract_terrain(dtm_m, fprints)

st_write(fprints_terrain, opt$output, append = FALSE)
message("Done. Output written to: ", opt$output)
