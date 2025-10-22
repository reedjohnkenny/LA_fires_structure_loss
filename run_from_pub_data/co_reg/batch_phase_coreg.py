## Script for aligning lidar and satellite imagery 

# Required arguments --input_pan: Input sattelite image that lidar will be aligned to
#                    --input_dem: Input digital elevation map derived from lidar, recomeded to be a canopy height model
#                    --output: file path for output lidar DEM that is aligned with satellite. 

# eaxample usage: python3 batch_phase_coreg.py --input_pan /path/to/pan.tif --input_dem /path/to/DEM --output /path_to_ouput_folder/DEM_coreg.tif




import argparse
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from skimage.registration import phase_cross_correlation
from scipy.ndimage import shift

def safe_normalize(img):
    img = np.nan_to_num(img, nan=np.nanmean(img))
    return (img - np.mean(img)) / np.std(img)

def run_coreg_phase(input_pan, input_dem, output_path):
    # Read WorldView PAN
    with rasterio.open(input_pan) as src_pan:
        pan = src_pan.read(1)
        pan_meta = src_pan.meta.copy()
        pan_transform = src_pan.transform
        pan_crs = src_pan.crs

    # Read and resample LiDAR DEM
    with rasterio.open(input_dem) as src_dem:
        dem = src_dem.read(1)
        
         # Set all NA (NaN) values to 0
        dem = np.nan_to_num(dem, nan=0)  # Replace NaN with 0
        
        dem_resampled = np.empty_like(pan, dtype=np.float32)
        reproject(
            source=dem,
            destination=dem_resampled,
            src_transform=src_dem.transform,
            src_crs=src_dem.crs,
            dst_transform=pan_transform,
            dst_crs=pan_crs,
            resampling=Resampling.bilinear
        )
        
        max_dem_resamp = dem_resampled.max()
    
        print(f"resamp dem max {max_dem_resamp}")

    # Normalize and estimate shift
    pan_norm = safe_normalize(pan)
    dem_norm = safe_normalize(dem_resampled)
    
    max_dem_norm = dem_norm.max()
    
    print(f"max dem value {max_dem_norm}")

    shift_yx, _, _ = phase_cross_correlation(pan_norm, dem_norm, upsample_factor=10)
    print(f"Estimated shift (y, x): {shift_yx}")

    # Apply shift
    dem_shifted = shift(dem_resampled, shift=shift_yx)

    # Save result
    pan_meta.update(dtype=rasterio.float32)
    with rasterio.open(output_path, 'w', **pan_meta) as dst:
        dst.write(dem_shifted.astype(np.float32), 1)

    print(f"✅ Output saved to: {output_path}")

# === CLI ===
def main():
    parser = argparse.ArgumentParser(description="Co-register one LiDAR DEM to one WorldView PAN using phase correlation.")
    parser.add_argument('--input_pan', required=True, help='Path to WorldView PAN .tif file')
    parser.add_argument('--input_dem', required=True, help='Path to LiDAR DEM .tif file')
    parser.add_argument('--output', required=True, help='Path to save co-registered DEM')

    args = parser.parse_args()
    run_coreg_phase(args.input_pan, args.input_dem, args.output)

if __name__ == "__main__":
    main()
