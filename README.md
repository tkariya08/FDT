# First Derivate Threshold (FDT) method

## Main code

* `Loess_FDT.R` : Rscript for generating 500 loess fitted residual bootstrap samples from discrete CO<sub>2</sub> flask measurement data at a given site and calculating FDT and zero-crossing derived CUP.

## Input Data:

CSV files containing CO<sub>2</sub> dry air mole fraction data from different NOAA flask measurement sites used in this study have the following file name format: `Flask_data_CO2_site.csv`, where ‘site’ represents the station name. For e.g., `Flask_data_CO2_mlo.csv`, has CO<sub>2</sub> dry air mole fraction data for MLO (Mauna Loa).

## Output:

The `Loess_FDT.R` program will generate FDT and zero-crossing derived CUP estimates for an ensemble of bootstrap samples.
