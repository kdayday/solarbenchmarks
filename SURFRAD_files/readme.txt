To jump right into forecasting, the pre-processed SURFRAD data is already available in the GHI_files/
folder. However, if you would like to modify the pre-processing itself, the raw SURFRAD data files
will need to be downloaded into these folders. As examples, one of the files is already provided
in each of the site/year/ folders. The remainder can be downloaded from:
ftp://aftp.cmdl.noaa.gov/data/radiation/surfrad/

The 2018 folders will need all 365 days downloaded (day 1 is provided as the example in each folder); 
the 2017 folders will need the last 20 days of the year downloaded to generate the 20-day persistence
ensemble forecasts (day 346 is provided as the example in each folder).

After downloading and unzipping the CAMS McClear folder, preprocess_SURFRAD.R can be run from the 
solarbenchmars/ folder to re-generate the NetCDF files in the GHI_files/ folder; this preprocessing
does the 5-minute or hourly average of the GHI and clear-sky GHI values to use in the very short-
term and short-term forecasts, respectivley. 