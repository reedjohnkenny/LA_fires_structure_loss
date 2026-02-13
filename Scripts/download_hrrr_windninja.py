#!/usr/bin/env python3
"""
Download HRRR f00 surface GRIB2 files for the Palisades fire extent
for use as WindNinja input.

Time range: 01/07/2025 11:00 AM – 01/08/2025 8:00 PM PST
            (01/07/2025 19:00 – 01/09/2025 04:00 UTC)

Extent derived from palisades_burned_dtm_2023_epsg32611.tif
"""

import os
from pathlib import Path

import pandas as pd
import rasterio
from rasterio.warp import transform_bounds
from herbie import Herbie

# -- Configuration --
DTM_PATH = "raw_data/palisades_burned_dtm_2023_epsg32611.tif"
OUT_DIR = Path("wind_data/palisades_hrrr")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Time range in UTC (PST + 8h)
START_UTC = "2025-01-07 19:00"
END_UTC = "2025-01-09 04:00"

# Buffer (degrees) around the DTM extent so WindNinja has context
BUFFER = 0.1

# -- Get domain extent from DTM --
with rasterio.open(DTM_PATH) as src:
    west, south, east, north = transform_bounds(src.crs, "EPSG:4326", *src.bounds)

west -= BUFFER
south -= BUFFER
east += BUFFER
north += BUFFER

print(f"Domain (buffered): W={west:.4f} S={south:.4f} E={east:.4f} N={north:.4f}")

# -- Download HRRR f00 for each hour --
dates = pd.date_range(START_UTC, END_UTC, freq="1h")

for date in dates:
    print(f"\n--- {date} UTC ---")
    H = Herbie(date, model="hrrr", product="sfc", fxx=0)

    # Download the full surface file subsetted to the domain
    try:
        path = H.download(
            search="(?:UGRD|VGRD):10 m|(?:TMP):2 m|(?:TCDC)|(?:HGT):surface",
            save_dir=OUT_DIR,
        )
        print(f"  Saved: {path}")
    except Exception as e:
        print(f"  ERROR: {e}")

print(f"\nDone. {len(dates)} timesteps downloaded to {OUT_DIR}/")
print(
    "\nTo convert to NetCDF for WindNinja, run:\n"
    "  for f in wind_data/palisades_hrrr/*.grib2; do\n"
    '    wgrib2 "$f" -netcdf "${f%.grib2}.nc"\n'
    "  done"
)
