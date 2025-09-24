library(terra)
library(stringr)

# --- Earth–Sun distance (AU) ---
earth_sun_distance_au <- function(date_time){
  doy <- as.integer(strftime(date_time, "%j"))
  g <- 2*pi*(doy - 1)/365
  1.00011 + 0.034221*cos(g) + 0.00128*sin(g) + 0.000719*cos(2*g) + 0.000077*sin(2*g)
}

# --- Parse IMD robustly (WV2/WV3) ---
parse_imd_robust <- function(imd_path){
  txt <- readLines(imd_path, warn = FALSE)
  txt <- paste(txt, collapse = "\n")
  
  # satId (WV02/WV03)
  satId <- str_match(txt, "satId\\s*=\\s*\"([^\"]+)\"")[,2]
  
  pat <- regex(
    "BEGIN_GROUP\\s*=\\s*BAND_([A-Z0-9]+)\\s*(.*?)\\s*END_GROUP\\s*=\\s*BAND_\\1",
    dotall = TRUE
  )
  
  blocks <- str_match_all(txt, pat)[[1]]
  
  stopifnot(nrow(blocks) > 0)
  
  cal <- lapply(seq_len(nrow(blocks)), function(i){
    code  <- blocks[i,2]
    block <- blocks[i,3]
    acf <- as.numeric(str_match(block, "absCalFactor\\s*=\\s*([0-9eE\\.\\+\\-]+)")[,2])
    ebw <- as.numeric(str_match(block, "effectiveBandwidth\\s*=\\s*([0-9eE\\.\\+\\-]+)")[,2])
    if (is.na(acf) || is.na(ebw))
      stop("Missing absCalFactor/effectiveBandwidth for band ", code)
    list(code = code, absCalFactor = acf, effectiveBandwidth = ebw)
  })
  names(cal) <- vapply(cal, `[[`, "", "code")
  
  # Sun elevation: try several keys
  sunElev <-
    suppressWarnings(as.numeric(str_match(txt, "sunElevation\\s*=\\s*([0-9eE\\.\\+\\-]+)")[,2]))
  if (is.na(sunElev))
    sunElev <- suppressWarnings(as.numeric(str_match(txt, "meanSunEl\\s*=\\s*([0-9eE\\.\\+\\-]+)")[,2]))
  if (is.na(sunElev))
    sunElev <- suppressWarnings(as.numeric(str_match(txt, "maxSunEl\\s*=\\s*([0-9eE\\.\\+\\-]+)")[,2]))
  if (is.na(sunElev)) stop("Could not find sun elevation in IMD.")
  
  # Acquisition time: firstLineTime or earliestAcqTime
  tstr <- str_match(txt, "firstLineTime\\s*=\\s*([0-9T:\\.\\-]+)Z")[,2]
  if (is.na(tstr))
    tstr <- str_match(txt, "earliestAcqTime\\s*=\\s*([0-9T:\\.\\-]+)Z")[,2]
  if (is.na(tstr)) stop("Could not parse acquisition time from IMD.")
  acq <- as.POSIXct(tstr, tz = "UTC", format = "%Y-%m-%dT%H:%M:%OS")
  
  list(cal = cal, sunElev = sunElev, acq = acq, satId = satId)
}

# ESUN tables (W·m^-2·µm^-1) for VIS-NIR
ESUN_WV2 <- c(C=1758.222, B=1970.376, G=1856.410, Y=1738.479,
              R=1559.455, RE=1342.069, N=1069.730, N2=861.286)
ESUN_WV3 <- c(C=1758.0,   B=2004.6,   G=1820.7,   Y=1551.3,
              R=1204.8,   RE=1053.1,  N=858.8,    N2=732.3)

# --- DN -> TOA reflectance ---
dn_to_toa_reflectance_maxar <- function(
    in_tif, imd, out_tif,
    band_codes = c("C","B","G","Y","R","RE","N","N2"),
    ESUN = NULL, overwrite = TRUE
){
  meta <- parse_imd_robust(imd)
  
  # pick ESUN automatically if not supplied
  if (is.null(ESUN)) {
    if (!is.null(meta$satId) && toupper(meta$satId) == "WV02") {
      ESUN <- ESUN_WV2
    } else if (!is.null(meta$satId) && toupper(meta$satId) == "WV03") {
      ESUN <- ESUN_WV3
    } else {
      stop("Unknown satId in IMD; please supply ESUN manually.")
    }
  }
  
  d <- earth_sun_distance_au(meta$acq)
  cos_theta <- cos((90 - meta$sunElev) * pi/180)
  
  r <- rast(in_tif)
  if (nlyr(r) != length(band_codes))
    stop("Input has ", nlyr(r), " layers but band_codes has ", length(band_codes))
  
  out <- rast(r)
  for (i in seq_len(nlyr(r))){
    code <- band_codes[i]
    if (is.null(meta$cal[[code]]))
      stop("Band code ", code, " not found in IMD. Parsed bands: ",
           paste(names(meta$cal), collapse = ", "))
    acf <- meta$cal[[code]]$absCalFactor
    ebw <- meta$cal[[code]]$effectiveBandwidth
    
    # DN -> radiance (W m^-2 sr^-1 µm^-1)
    L <- app(r[[i]], function(x) x * (acf/ebw), cores = 1)
    
    # radiance -> TOA reflectance (unitless)
    rho <- app(L, function(x) (pi * x * d^2) / (ESUN[code] * cos_theta), cores = 1)
    names(rho) <- names(r[[i]])
    out[[i]] <- rho
  }
  
  writeRaster(out, out_tif, datatype="FLT4S",
              gdal=c("COMPRESS=DEFLATE","PREDICTOR=3","TILED=YES"),
              overwrite=overwrite)
  invisible(list(out=out_tif, meta=meta))
  return(out)
}

# correct Palisades pre

toa <- dn_to_toa_reflectance_maxar(
  in_tif  = "~/Desktop/Urban_tree_fire/landscape_analysis/Data/Maxar_raw/Maxar_5_99km2___California__USA_25318O6E-1_DAY_California-USA/200008990764_01/200008990764_01_P001_MUL/Palisades_unburn_MUL_ortho.TIF",
  imd     = "~/Desktop/Urban_tree_fire/landscape_analysis/Data/Maxar_raw/Maxar_5_99km2___California__USA_25318O6E-1_DAY_California-USA/200008990764_01/200008990764_01_P001_MUL/24OCT20184322-M2AS-200008990764_01_P001.IMD",
  out_tif = "~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Palisades/Palisades_pre/Multi_spec/24OCT20184322-M2AS-200008990764_01_P001_TOA.TIF",
  band_codes = c("C","B","G","Y","R","RE","N","N2")  # your order
)



# correct eaton post

dn_to_toa_reflectance_maxar(
  in_tif = "~/Desktop/Urban_tree_fire/landscape_analysis/Data/Maxar_raw/Maxar_9_5km2___California__USA_2519990J-1_DAY_California-USA/200007937065_01_200007937065_01_P001_MUL_25JAN14183036-M2AS-200007937065_01_P001.TIF", 
  imd = "~/Desktop/Urban_tree_fire/landscape_analysis/Data/Maxar_raw/Maxar_9_5km2___California__USA_2519990J-1_DAY_California-USA/200007937065_01_200007937065_01_P001_MUL_25JAN14183036-M2AS-200007937065_01_P001.IMD",
  out_tif = "~/OneDrive - Cal Poly/Advanced_GIS_data/Maxar/Eaton/POST/Multi-spec/200007937065_01_200007937065_01_P001_MUL_25JAN14183036-M2AS-200007937065_01_P001_TOA.TIF",
  band_codes = c("C","B","G","Y","R","RE","N","N2") 
)

