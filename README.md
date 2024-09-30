# Arctic-TESLA

[![DOI](https://zenodo.org/badge/13863196.svg)](https://zenodo.org/doi/10.5281/zenodo.13863196)

ArcticTESLA is a collection of Python3 and Matlab functions for generating stochastic wave and water level scenarios at coastal locations in Northern Alaska on the Chukchi and Beaufort Seas.
The package creates new time series of forcing conditions by generating new possible synoptic weather chronologies. 
The workflow identifies historical synoptic weather patterns, and the meteorologic and oceanic conditions that occurred during those weather systems.
Markov chains for each pattern, as well as their likelihood of occurrence conditional dependent on large scale climate indices, are then used in monte carlo simulations.
All codes are used in Anderson and Cohn (in review), with some adopted from Anderson et al. (2019) and the references within. 
Many of the codes within this repo are functions from https://github.com/teslakit/teslakit with the goal of this repo eventually being merged into teslaKit.

## Main contents

### Modules:

Data Download and Extraction
- [era5](./dataDownloads/era5metOceanDownloads.py) Automated downloads from Copernicus (follow steps outlined below in 'Downloading Data' to ensure account access)
- [noaa](./dataDownloads/export_local_tides_alaska.m) Automated downloads of the three tide gauges available in Northern Alaska from NOAA Tides and Currents
- [nsidc](./dataDownloads/downloadIce2.py) Automated downloading of the NSIDC daily SIC dataset.

Main Scripts
- [awts](./awts.py) All codes relevant to synthesizing ERSSTv5 into annual sea surface temperature patterns and projecting future climate variability.
- [sic](./sic.py) All codes relevant to synthesizing daily SIC fields and wave basin sizes and projecting future variability.
- [mjo](./mjo.py) All codes relevant to synthesizing outgoing longwave radiation and projecting variability.
- [dwts](./dwts.py) All codes relevant to synthesizing SLPs into classified daily weather patterns.
- [temps](./temps.py) All codes relevant to synthesizing ERA5 surface temperatures and projecting future Arctic air temperatures.

Simulations
- [historical SLP](./simulations/historicalSimulations.py) Application of the ALR model for DWTs to the historical time frame to validate model fits.
- [future SLP](./simulations/futureSLPsimulations.py) Application of the ALR model for DWTs out to 2075.
- [downscaling](./simulations/historicalHydrographsInterpolated.py) Statistically downscaling environmental variables by stretching hydrographs to fit simulating SLP time series.

Functions
- [alr](./functions/alr.py) AutoRegressive Logistic Model customized wrapper
- [mda](./functions/mda.py) Maximum Dissimilarity Algorithm (MDA) code used in multiple scripts.
- [metOcean](./functions/metOcean.py) All codes to separate sea states and define joint probabilities

# Downloading Data
### Sea Level Pressure Fields
Python3 and Matlab functions are provided for loaded CFSR SLP fields that are downloaded to a local directory. The user must download all monthly files since 1979 found at the following links:

https://rda.ucar.edu/datasets/ds093.1/#description
https://rda.ucar.edu/datasets/ds094.1/#description

You will need to chose 'All available' under the data tab, at full resolution in lat/lon, and click the checkbox to convert the download from grib to netcdf inorder to work with the built-in functions in this library.

### Sea Ice Concentration Fields
The user will need to download all northern hemipshere .bin SIC files from the National Snow and Ice Data Center at https://nsidc.org/data/nsidc-0081/versions/2. 
The nc2bin_siconc.py function within dataDownloads will then convert all files to the format expected by sic.py.

### Local Waves, Winds, and Temperatures
To use ERA5 hindcasts, this package requires access to the online Thredds server hosted by Copernicus.

1. Create an account with Copernicus by signing up here.
2. Once you have an account, sign in to your Copercius account here and note the UID and API key at the bottom of the page.
3. Paste the code snippet below into your terminal, replacing 'UID' and 'API' with those from step 2:

(echo 'url: https://cds.climate.copernicus.eu/api/v2';
  echo 'key: UID:API';
  echo 'verify: 0';
   ) >> ~/.cdsapirc

The above command creates the file ~/.cdsapirc with your API key, which is necessary to use the CDS API. As a sanity check, use more ~/.cdsapirc to ensure everything appears correct.
Now you can run 'era5metOceanDownloads.py'.

### Water Levels
Use export_local_tides_alaska.m to download the latest NOAA tide gauge data from around Alaska and interpolate it to Point Hope, AK.

### Climate Variables: SST, OLR, T2M
Large-scale climate is accounted for by spatial patterns of the nearby ocean's sea surface temperature (SST) pattern at an annual scale. The awt.py file expects data in the yearly format available at https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NCDC/.ERSST/.version5/index.html?Set-Language=en. 
And the mjo.py script expects data in the format available at: https://iridl.ldeo.columbia.edu/SOURCES/.BoM/.MJO/.RMM/index.html?Set-Language=en.
Download monthly global surface temperatures for every node at 2-meters (T2M) from ERA5 through Copernicus. Save this as 'ERA5MonthlyTemps.nc'.

# Methods
The TESLA framework requires the below steps to be run sequentially as saved outputs from earlier steps become the inputs for later steps.

1. After downloading CFSR SLPs, provide the data directory to CFSR_extractSLPs_rectify_cropLand.m to save 'slps.mat', run dwts.py which will save 'dwts.pickle'.
2. After downloading EFA5 montly surface temperatures, provide the correct directory and run temps.py to create 'predictedArcticTemps.pickle'.
3. After downloading NSIDC SICs and converting with nc2bin_siconc.py, run sic.py which will create 'ice.pickle'.
4. Run 'export_local_tides_alaska.m' to create 'noaaAlaskaTides.mat', then run 'waterlevel.py' to create 'historicalNTR.pickle'.
5. After downloading ERA5 waves and NOAA tides to create 'waves.pickle', run hydrographs.py
5. After downloading ERSSTv5 SSTs, run awt.py to create 'pastAWTs.pickle' and 'futureAWTs.pickle'
6. After downloading BOM MJO indices, run mjo.py
7. Run copulas.py
8. Run futureIceSimulations.py
9. Run futureSLPsimulations.py to create 'dwtFutureSimulation.pickle'


Anderson, D. and N. Cohn (in review) Future coastal tundra loss due to compounding environmental changes in Alaska.

Anderson, D., A. Rueda, L. Cagigal, J. Antolinez, F. Mendez, and P. Ruggiero. (2019) Time-varying Emulator for Short and Long-Term Analysis of Coastal Flood Hazard Potential. Journal of Geophysical Research: Oceans, 124(12), 9209-9234. https://doi.org/10.1029/2019JC015312