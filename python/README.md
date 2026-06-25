## Mobility Processing Code

This folder contains code for downloading the Advan Neighborhood Patterns Plus dataset, and processing this data into sewershed-level activity data.


### Steps

1. `ww_mobility_download.sh` is a SLURM batch script which calls `download.py` 
2. `download.py` will download the Advan Neighborhood Patterns Plus dataset, which consists of individual CSVs for each month/year. 
3. `2026_04_20_mobility_processing.ipynb` collects the individual CSVs and processes them into monthly and daily level datasets for input into our models.
4. (Unused) `2026_05_07_detrend_daily.ipynb` is a notebook which experiments with two detrending methods for the mobility data. Neither method was used in our final analysis, and this file is just included for completeness.

