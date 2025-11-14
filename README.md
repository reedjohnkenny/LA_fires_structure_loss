# LA_fires_structure_loss

This repository contains the workflow and scripts used to process spatial data and analyze how the distribution of urban trees and buildings relates to structure loss in recent Los Angeles wildfires.

## Overview

The analysis is implemented as a Snakemake pipeline that prepares spatial datasets, computes wind metrics, and generates model inputs for subsequent statistical modeling in R (v4.4.2). All required R packages can be installed using `install_packages.R`.

Because some source data are not publicly available, the workflow begins at the earliest point where shareable inputs can be used.


## Required Input Files

*Items marked “public” can be obtained from publicly available sources.*

- Tree canopy polygons (burned areas): `{fire}_burned_crowns_polygons_v2.shp`
- Building footprints: `LARIAC6_Buildings_2020_{fire}.shp`
- Study boundary polygons: `{fire}_burned_bound.shp`
- CA fire perimeters, 2025 (public): `Perimeters.shp`
- DINS structure points (public): `DINS_2025_{fire}_Public_View.geojson`
- 

You will need to place these in the expected directory structure or update file paths accordingly.

## Running the Workflow

### 1. Prepare building and tree datasets

If you obtained the archived input data from `{archive}`, update file paths in:

- `pipeline_scripts/Create_buildings_final_dataset.R`
- `pipeline_Snakefile`

Then execute the main pipeline:

snakemake -s pipeline_Snakefile -j 1 -p

This step generates the model-input datasets used in later analyses.

### 2. Compute wind metrics

Wind statistics are derived using two steps:

1. Run `sites_for_wind.R` after correcting input paths.  
2. Run `wind_data.ipynb`, ensuring the path to `all_burned_wind_sites.csv` is correct.

### 3. Generate model inputs and model selection

- Run `gen_altinputs.R` to create alternative input variants.  
- Run `model_selection.R` to fit candidate models and select the preferred specification.

### 4. Produce figures and tables

Run `Figures.Rmd` to reproduce the figures and summary tables used in the analysis.

### 5. Additional scripts

- `prediction_vis.R`: Produces tables summarizing predicted changes in structure-loss probability across mean versus minimum values of predictors.
- `validation/Confusion_matrix`: Performs accuracy assessment of the tree canopy polygons.

## Data Dictionary for Model Inputs

### Categorical variables

- `DAMAGE`  
- `destroyed` (binary)

### Numerical variables

`utm_x`, `utm_y`, `distance_to_nearest_building`, `distance_to_nearest_tree`,  
`area_of_trees_2m`, `tree_dens`, `tree_dens_300`, `tree_dens_500`,  
`build_dens`, `build_dens_300`, `build_dens_500`,  
`angular_diff_build`, `angular_diff_tree`

---

## Variable Definitions

- **DAMAGE**: One of  
  `Inaccessible`,  
  `No Damage`,  
  `Affected (1–9%)`,  
  `Minor (10–25%)`,  
  `Major (26–50%)`,  
  `Destroyed (>50%)`

- **destroyed**:  
  - `0` = “not destroyed” (`No Damage`, `Affected (1–9%)`)  
  - `1` = “destroyed” (`Minor`, `Major`, `Destroyed (>50%)`)  
  
- **utm_x / utm_y**: UTM coordinates in EPSG:32611.

- **distance_to_nearest_building**: Euclidean distance (m) from each structure to the nearest neighboring structure.

- **distance_to_nearest_tree**: Distance (m) from the structure boundary to the nearest tree‐canopy polygon.

- **area_of_trees_2m**: Total canopy area (m²) within 2 m of a structure polygon.

- **tree_dens**, **tree_dens_300**, **tree_dens_500**: Kernel-density estimates of tree canopy using 200 m, 300 m, and 500 m bandwidths.

- **build_dens**, **build_dens_300**, **build_dens_500**: Kernel-density estimates of building footprints using the same respective bandwidths.

- **angular_diff_build**: Absolute angular difference (degrees) between mean wind direction and the bearing from a structure to its nearest building.

- **angular_diff_tree**: Same definition as above, but for the nearest tree.

---

## Pre-processing of licensed data

Pre-processing of licensed data (LiDAR) and sattelite imagery is documented by `lidR_create_chms.R` as well as scripts in the `co_reg` and `Preprocess_sat` directories. Point cloud processing of LiDAR to create canopy height models is done in `lidar_create_chms.R`. Alignment of LiDAR and satellite imagery is done in `co-reg/batch_phase_coreg.py` and `coreg_snakefile` then tiles are stitched using `Stich_lidar.R`. Satellite imagery is stiched together and standardized for NDVI calculations by running the `Stich_maxar.R`, `TOA_reflectance.R` and `correct_NDVI.Rmd` scripts in the `Preprocess_sat` directory. Raw tree canopy polygons are generated in the `canopy_seg` rule of the `pipeline_snakefile` script using the `Canopy_seg.R` script. 

