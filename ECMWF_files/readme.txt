This folder contains a framework to download and pre-process ECMWF data for the SURFRAD locations. 
This requires 3 steps:
1. Downloading ECMWF data in .GRIB format from their MARS file system, over entire CONUS 
   for the year 2018
2. Preprocessing from .GRIB format to NetCDF format (done for ease of use in R)
3. Preprocessing from accumulated irradiance to hourly irradiance and spatially interpolate 
   to get hourly forecasts at SURFRAD locations

These scripts assume a Windows rather than Linux system. The scripts will require some 
modification; for instance, the ECMWF_requests.bat file must be updated to your local 
installation of MARS (line 40):
https://confluence.ecmwf.int/display/SUP/2019/02/15/ECMWF+software+under+Windows

Instructions on each of these steps:
1. Refer to "GRIB files\ECMWF_requests.bat" to download daily files of the 50 members of 
   the perturbed forecast and "GRIB control files\ECMWF_ctrl_requests.bat" to download 
   the 1 control forecast member. Due to ECMWF limits to 20 requests per user at a time, 
   these 2 sets of 365 files are requested in batches. Within each .bat files, the user 
   can define the current batch through the upper and lower limits on the index i (l_limit 
   and u_limit), where i is the number of the day of the year (i.e., i is from 1 to 365). 
   These will populate the "GRIB files\" and "GRIB control files\" (alternatively, if the 
   computer loses connection with MARS during the request process, the files may need to 
   be manually downloaded from the MARS website and renamed appropriately). Due to wait 
   times for ECMWF's MARS file system, the entire process may take several weeks.
2. Refer to the "grib_to_nc.bat" and "grib_ctrl_to_nc.bat" batch scripts to translate the 
   GRIB files into NetCDF; like the request scripts, this is also accomplished in parallel 
   batches using l_limit and u_limit indices. This process is on the order of minutes, 
   not weeks. This batch scripts will populate the "Raw/" folder with both the perturbed
   and control .nc files. The NetCDF file format can be easily read into R with the ncdf4
   package. 
3. Finally, run the preprocess_ECMWF.R script in the main solarbenchmarks/ folder to do 
   time averaging and spatial averaging in preparation for the forecast comparison. This
   will populate the "SURFRAD_sites\" folder with one year-long .nc file for each of the 
   7 SURFRAD locations. 