library(terra)
library(sf)
library(spatialEco)

# 1) Common reference grid over the union bbox (or your study area mask)
ref <- rast(ext(vect(st_union(st_geometry(Eaton_tree_crowns)))), 
            resolution = 30, crs = st_crs(Eaton_tree_crowns)$wkt)

# 2) KDE as probability density (units ~ 1/m^2), with clean scaling
kde_all  <- sf.kde(st_centroid(Eaton_tree_crowns),  bw = 200, res = 30,
                   ref = ref, standardize = FALSE, scale.factor = 1)
kde_burn <- sf.kde(st_centroid(eaton_burn_crowns),  bw = 200, res = 30,
                   ref = ref, standardize = FALSE, scale.factor = 1)

# 3) Convert pdf to intensity (points per m^2): lambda(x) = n * pdf(x)
n_all  <- nrow(st_as_sf(Eaton_tree_crowns))
n_burn <- nrow(st_as_sf(eaton_burn_crowns))

int_all  <- kde_all  * n_all
int_burn <- kde_burn * n_burn

# Now maxima are comparable as intensity (points/m^2) on the same grid
global(int_all, "max", na.rm = TRUE)
global(int_burn, "max", na.rm = TRUE)




library(sf)
library(terra)
library(spatialEco)

# --- inputs (sf POINTS in an equal-area CRS) ---
# cityA_pts, cityB_pts : sf POINT objects for each city
# cityA_boundary, cityB_boundary : sf POLYGONs for each city boundary (same CRS)

# 1) Choose a common bandwidth (meters)
bw_m <- 200  # or pick via spatstat (see section B) and reuse

# 2) Build identical-resolution rasters over each city's boundary
make_ref <- function(boundary, res = 30) {
  rast(ext(vect(boundary)), resolution = res, crs = st_crs(boundary)$wkt)
}

ref_eaton <- make_ref(Eaton_tree_crowns, res = 30)
ref_pal <- make_ref(pal_tree_crowns, res = 30)

# 3) KDE as probability density (disable arbitrary scaling)
k_eaton_pdf <- sf.kde(
  x = st_centroid(Eaton_tree_crowns),
  bw = bw_m, ref = ref_eaton, standardize = FALSE, scale.factor = 1, mask = FALSE, res = 30
)

k_pal_pdf <- sf.kde(
  x = st_centroid(pal_tree_crowns),
  bw = bw_m, ref = ref_pal, standardize = FALSE, scale.factor = 1, mask = FALSE, res = 30
)

# 4) Convert to intensity (points per m^2): lambda(x) = n * pdf(x)
nA <- nrow(cityA_pts)
nB <- nrow(cityB_pts)

kA_int_m2 <- kA_pdf * nA
kB_int_m2 <- kB_pdf * nB

# 5) Convert to points per hectare
kA_tpha <- kA_int_m2 * 10000
kB_tpha <- kB_int_m2 * 10000

# --- Sanity check: integral of intensity ≈ number of points ---
cell_area <- res(refA)[1] * res(refA)[2]        # m^2
sumA <- global(kA_int_m2, "sum", na.rm = TRUE)[,1] * cell_area
sumB <- global(kB_int_m2, "sum", na.rm = TRUE)[,1] * cell_area
c(checkA = sumA, pointsA = nA, checkB = sumB, pointsB = nB)

