#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(dplyr)
  library(tidyr)
  # …etc…
})

gen_mod_input <- function(status, comp, build_dens, tree_dens, build_vols, tree_prox, build_prox, chm_diff, ndvi_diff, pre_height) {
    mod_input <- as.data.frame(build_prox) %>% 
      mutate(lon = st_coordinates(st_transform(st_centroid(geometry), crs = 4326))[,1],
             lat = st_coordinates(st_transform(st_centroid(geometry), crs = 4326))[,2]) %>% 
      left_join(build_dens, by = "UID") %>% 
      mutate(utm_x = st_coordinates(st_centroid(geometry))[,1], 
             utm_y = st_coordinates(st_centroid(geometry))[,2], 
             area_of_trees_2m = replace_na(ar_tr_2, 0),
             area_of_trees_5m = replace_na(ar_tr_5, 0), 
             area_of_trees_10m = replace_na(ar_t_10, 0)) %>% 
      select(UID, lon, lat, utm_x, utm_y, DAMAGE = DAMAGE.x, distance_to_nearest_building = dst2n_b, bearing_to_nearest_building = brng_t_c_, distance_to_nearest_tree = dst2n_t, bearing_to_nearest_tree = brng_t_n_, number_of_trees_2m = nm_tr_2, number_of_trees_5m = nm_tr_5, number_of_trees_10m = nm_t_10, area_of_trees_2m, area_of_trees_5m, area_of_trees_10m, area_ext_build_2 = ar_x__2, area_ext_build_10 = ar___10, area_ext_build_25 = ar___25, area_ext_build_50 = ar___50, area_ext_build_100 = a___100, tree_dens, build_dens)

  return(mod_input)
}

# define CLI options
option_list <- list(
  make_option(c("-b","--buildings_prox"),type="character", help="Path to building prox"),
  make_option(c("-D","--build_density"),type="character", help="Path to building dens"),
  make_option(c("-o","--output"),   type="character", default="proximities.geojson",
              help="Output path [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# validate
if ( is.null(opt$buildings_prox)) {
  print_help(OptionParser(option_list=option_list))
  quit(status=1)
}



build_dens <- st_read(opt$build_density, quiet = T)


build_prox <- st_read(opt$buildings_prox, quiet = T)


mod_input <- gen_mod_input(build_dens = build_dens, build_prox = build_prox)

write.csv(mod_input, opt$output)

