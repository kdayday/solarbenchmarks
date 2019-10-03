# Author: Kate Doubleday
# Last updated: 10/3/2019
# This script pre-processes ECMWF NETCDF files to extract their 51 ensemble forecasts of GHI from the SSRD (surface solar radiation downwards) ECMWF parameter.
# Estimates for the 7 SURFRAD sites are interpolated from the nearest 4 grid points, given a 0.2x0.2 degree lat/lon grid. 
# GHI forecast time resolution is hourly average.
# -----------------------------------------------------------------
# Load dependencies
library(here)
library(ncdf4)
library(pracma)

# -----------------------------------------------------------------
# Define constants

site_names <- c("Bondville",
                "Boulder",
                "Desert_Rock",
                "Fort_Peck",
                "Goodwin_Creek",
                "Penn_State",
                "Sioux_Falls")

site_coord <- list(Bondville=list(N=40.05192, W=88.37309),
                            Boulder=list(N=40.12498, W=105.23680),
                            Desert_Rock=list(N=36.62373, W=116.01947),
                            Fort_Peck=list(N=48.30783, W=105.10170),
                            Goodwin_Creek=list(N=34.2547, W=89.8729),
                            Penn_State=list(N=40.72012, W=77.93085),
                            Sioux_Falls=list(N=43.73403, W=96.62328))

input_directory <- here("NetCDF_files", "Raw")
output_directory <- here("NetCDF_files", "SURFRAD_sites")
dir.create(output_directory, showWarnings = FALSE)

day_indx <- 1:31
month_days <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

# Parameters for new NetCDF file
dims <- mapply(ncdim_def, name=c('Day', 'IssueTime', "LeadTime", 'Member'), units=c('', 'Hour', "Hours", ""), 
               vals=list(day_indx, seq(0, 18, by=6), 1:24, 1:50), 
               longname=c("Day number starting Jan 1, 2018", "Forecast time of issue", "Forecast lead time", "Ensemble member"), 
               SIMPLIFY = FALSE)
pvar <- ncvar_def("irradiance", "W/m^2", dims, missval=NA, compression = 9)

# -----------------------------------------------------------------
#' A helper function to import each daily netCDF, spatially interpolate, and transform to irradiance
#' @param site Site name
#' @param month Month number
#' @param day Day of month number
#' @return an array of forecasts
get_daily_irr <- function(day, month, site) {
  xp <- site_coord[[site]]$N
  yp <- site_coord[[site]]$W
    
  # Extract closest spatial data from NetCDF file
  nc <- nc_open(file.path(input_directory, paste("north_america_", ifelse(month<10, paste(0, month, sep=''), month), ifelse(day<10, paste(0, day, sep=''), day), ".nc", sep='')))
  lon_indx <- which(round(nc$dim$longitude$vals, 1)==-ceiling(yp*5)/5)
  lat_indx <- which(round(nc$dim$latitude$vals, 1)==ceiling(xp*5)/5)
  ssrd_grid <- ncvar_get(nc, varid="ssrd", start=c(lon_indx,lat_indx,1,1,1), count=c(2,2,-1,-1,-1))
  nc_close(nc)
  
  # Spatially interpolate to SURFRAD coordinates. 
  x <- nc$dim$latitude$vals[c(lat_indx+1, lat_indx)] # lat coordinates must be nondecreasing
  y <- nc$dim$longitude$vals[c(lon_indx, lon_indx+1)]
  # interp2 requires Z to length(y) x length(x) rather than length(x) x length(y). 
  ssrd <- apply(ssrd_grid, MARGIN = c(3, 4, 5), FUN=function(z) (interp2(x, y, matrix(z[c(3,4,1,2)], ncol=2), xp, -yp, method="linear")))
  
  # Calculate irradiance from ssrd 
  irr <- apply(ssrd, MARGIN=c(1,3), FUN=diff)/3600
  irr <- apply(irr, MARGIN=c(1,2,3), FUN=round)
  # Clean: truncate small negative values to 0
  irr[irr<0] <- 0
  return(irr)
}

# -----------------------------------------------------------------
#' A helper function to gather the daily files into a single irradiance array for a given site and export to a new NetCDF
#' @param site Site name
export_site_netcdf <- function(site) {
  # Get irradiance array for this site
  irr <- sapply(day_indx, FUN=function(d) {get_daily_irr(d, month=min(which(d <=cumsum(month_days))), site=site)}, simplify="array")
  # Reformat to dimensions: [day]x[issue time]x[forecast horizon]x[members]
  irr <- aperm(irr, c(4,3,1,2))
  
  # Export new NetCDF file
  nc_new <- nc_create(file.path(output_directory, paste(site, ".nc", sep='')), pvar)
  ncvar_put(nc_new, pvar, irr, count=pvar[['varsize']])
  nc_close(nc_new)
}

# -----------------------------------------------------------------
# Cycle through the sites, exporting their site-specific irradiance NetCDF's
for (site in site_names) export_site_netcdf(site)
