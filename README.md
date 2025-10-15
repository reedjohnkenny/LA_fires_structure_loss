# LA_fires_structure_loss
Analysis of the effect of spatial distribution of urban trees and structures on structure loss in the LA fires

## Files needed to reproduce analyses 
*(some data files cannot be shared as they are not public. This workflow begins with files that can be shared)
- Raw tree canop polygons: {fire}_burned_crowns_polygons_v2.shp
- Buildings polygons: LARIAC6_Buildings_2020_{fire}.shp
- Study boundaries: {fire}_burned_bound.shp
- CA fire perimeters 2025 (publicly available): Perimeters.shp
- DINS points (publicly available): DINS_2025_{fire}_Public_View.geojson

## Order in which to run scripts

Assuming you have downloaded the above data archived at {archive} you will need to add the appropriate file paths to pipeline_scripts/Create_buildings_final_dataset.R, as well as pipeline_Snakefile. You can then run pipeline_Snakefile using a command like this "snakemake -s pipeline_Snakefile -j 1 -p". That will generate the model inputs. 

Next you need to calculate the wind speeds and directions. First run sites_for_wind.R after correcting the file paths, then run wind_data.ipynb (making sure that the path to all_burned_wind_sites.csv is correct). 

Now you can generate the models. Run create_models.R, to generate the figures, run Figures.Rmd. prediction_vis.R is a helper script that generates tables showing the changes in probability of structure loss as individual predictor variables change. 

