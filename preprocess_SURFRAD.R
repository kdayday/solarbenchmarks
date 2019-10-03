# Author: Kate Doubleday
# Last updated: 10/3/2019
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
ivar <- ncvar_def("irradiance", "W/m^2", dims, missval=NA, compression = 9)
svar <- ncvar_def("sun_up", "", dims, compression = 9, prec="integer")

# -----------------------------------------------------------------
#' A helper function to import each daily SURFRAD file and average hourly
#' @param site Site name
#' @param month Month number
#' @param day Day of month number
#' @return an array of forecasts
get_daily_data <- function(fname, site) {
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
  GHI <- sapply(1:24, FUN=function(i) {round(mean(GHI_min[(60*(i-1)+1):(60*i)]))})
  
  zen_min <- all_data[, "zen"]
  sun_up <- sapply(1:24, FUN=function(i) {any(zen_min[(60*(i-1)+1):(60*i)]<90)})
  
  return(matrix(c(GHI, sun_up), ncol=2))  
}

# -----------------------------------------------------------------
#' A helper function to gather the daily files into a single irradiance array for a given site and export to a new NetCDF
#' @param site Site name
export_site_netcdf <- function(site) {
  daily_files <- list.files(file.path(input_directory, site), pattern=".dat")
 
  # Get hourly irradiance and indicator if sun is up, based on solar zenith angle
  hourly_data <- sapply(daily_files, FUN=get_daily_data, site=site, simplify="array")
  irr <- t(hourly_data[,1,])
  sun_up <- t(hourly_data[,2,])
   
  # Export new NetCDF file
  nc_new <- nc_create(file.path(output_directory, paste(site, ".nc", sep='')), list(ivar, svar))
  ncvar_put(nc_new, ivar, irr, count=ivar[['varsize']])
  ncvar_put(nc_new, svar, sun_up, count=svar[['varsize']])
  nc_close(nc_new)
}

# -----------------------------------------------------------------
# Cycle through the sites, exporting their site-specific irradiance NetCDF's
for (site in site_names) export_site_netcdf(site)
