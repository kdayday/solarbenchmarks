# Author: Kate Doubleday
# Last updated: 10/2/2019
# This script pre-processes SURFRAD files to make hourly averages that match the hourly average ECMWF forecasts. 
# -----------------------------------------------------------------
# Load dependencies

library(here)
library(ncdf4)

# -----------------------------------------------------------------
# Define constants

input_directory <- here("SURFRAD_files", "Raw")
output_directory <- here("SURFRAD_files", "Hourly")
dir.create(output_directory, showWarnings = FALSE)

site_names <- c("Bondville_IL",
                "Boulder_CO",
                "Desert_Rock_NV",
                "Fort_Peck_MT",
                "Goodwin_Creek_MS",
                "Penn_State_PA",
                "Sioux_Falls_SD")

# Parameters for new NetCDF file
dims <- mapply(ncdim_def, name=c('Day', 'Hour'), units=c('', ''), 
               vals=list(1:365, 1:24), 
               longname=c("Day number starting Jan 1, 2018", "Hour of day"), 
               SIMPLIFY = FALSE)
pvar <- ncvar_def("irradiance", "W/m^2", dims, missval=NA, compression = 9)

# -----------------------------------------------------------------
#' A helper function to import each daily SURFRAD file and average hourly
#' @param site Site name
#' @param month Month number
#' @param day Day of month number
#' @return an array of forecasts
get_daily_irr <- function(fname, site) {
  all_data <- read.table(file.path(input_directory, site, fname), skip=2, col.names=c(
    "year","jday","month","day","hour","min",'dt', "zen","dw_solar","qc_dwsolar","uw_solar","qc_uwsolar","direct_n",
    "qc_direct_n","diffuse","qc_diffuse","dw_ir", "qc_dwir","dw_casetemp","qc_dwcasetemp","dw_dometemp",
    "qc_dwdometemp","uw_ir","qc_uwir","uw_casetemp", "qc_uwcasetemp","uw_dometemp","qc_uwdometemp","uvb",
    "qc_uvb","par","qc_par","netsolar","qc_netsolar","netir", "qc_netir","totalnet","qc_totalnet",
    "temp","qc_temp","rh", "qc_rh","windspd","qc_windspd","winddir","qc_winddir","pressure","qc_pressure"))
  
  GHI_min <- all_data[,"dw_solar"]
  # Clean: truncate small negative values to 0
  GHI_min[GHI_min<0] <- 0

  # Average to hourly to match ECMWF forecasts
  GHI <- sapply(1:24, FUN=function(i) {mean(GHI_min[(60*(i-1)+1):(60*i)])})
  
  return(GHI)  
}

# -----------------------------------------------------------------
#' A helper function to gather the daily files into a single irradiance array for a given site and export to a new NetCDF
#' @param site Site name
export_site_netcdf <- function(site) {
  daily_files <- list.files(file.path(input_directory, site), pattern=".dat")
  
  # Get irradiance array for this site
  irr <- t(sapply(daily_files, FUN=get_daily_irr, site=site, simplify="array"))
  
  # Export new NetCDF file
  nc_new <- nc_create(file.path(output_directory, paste(site, ".nc", sep='')), pvar)
  ncvar_put(nc_new, pvar, irr, count=pvar[['varsize']])
  nc_close(nc_new)
}

# -----------------------------------------------------------------
# Cycle through the sites, exporting their site-specific irradiance NetCDF's
for (site in site_names) export_site_netcdf(site)
