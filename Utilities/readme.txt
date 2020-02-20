This folder contains example batch scripts that show how ECMWF data was requested and processed to NetCDF, based on a Windows rather than Linux system.
The scripts will require some modification; for instance, the ECMWF_requests.bat file must be updated to your local installation of MARS (line 40):
https://confluence.ecmwf.int/display/SUP/2019/02/15/ECMWF+software+under+Windows

These batch scripts request and process GRIB files for each day of the year over 2018. Due to ECMWF limits to 20 requests per user at a time, these 365 files are requested in batches.
The user can define the current batch through the upper and lower limits on the index i (l_limit and u_limit), where i is the number of the day of the year (i.e., i is from 1 to 365). 

Once the user has downloaded the entire year's worth of GRIB files, the files can be translated NetCDF file format using the grib_to_nc.bat batch file (also accomplished in 
parallel batches). The NetCDF file format can be easily read into R with the ncdf4 package. 

Following these two steps, the data can be preprocessed with the preprocess_ECMWF.R script to do time averaging and spatial averaging in preparation for the forecast comparison. 