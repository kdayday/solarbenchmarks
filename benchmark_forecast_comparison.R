# Author: Kate Doubleday
# Last updated: 10/10/2019
# -----------------------------------------------------------------
# Load dependencies

library(here)
library(ncdf4)
library(pracma)
library(truncnorm)
library(ggplot2)
library(reshape2)
source("evaluation_functions.R")
source("forecast_methods.R")

# -----------------------------------------------------------------
# Define constants

telemetry_directory <- here("SURFRAD_files", "Hourly")
forecast_directory <- here("NetCDF_files", "SURFRAD_sites")
output_directory <- here("Results")
dir.create(output_directory, showWarnings = FALSE)

site_names <- c("Bondville",
                "Boulder",
                "Desert_Rock",
                "Fort_Peck",
                "Goodwin_Creek",
                "Penn_State",
                "Sioux_Falls")

percentiles <- seq(0.01, 0.99, by=0.01)
intervals <- seq(0.1, 0.9, by=0.1) # Central intervals for sharpness evaluation
num_peen <- 20 # Number of members in the persistence ensemble

forecast_names<- c("Climatology", "Ch-PeEn", "ECMWF Ensemble", "PeEn", "Gaussian")

metric_names <- c("CRPS", "Left-tail weighted CRPS", "Center weighted CRPS", "Right-tail weighted CRPS")

# Prepare results data frames
metrics_df <- array(dim=c(length(site_names), length(forecast_names), length(metric_names)), dimnames=list(site_names, forecast_names, metric_names))
reliability_df <- array(dim=c(length(site_names), length(forecast_names), length(percentiles)), dimnames=list(site_names, forecast_names, percentiles))
interval_width_df <- array(dim=c(length(site_names), length(forecast_names), length(intervals)), dimnames=list(site_names, forecast_names, intervals))

# -----------------------------------------------------------------
get_site_data <- function(site, metrics_df, reliability_df, interval_width_df) {
  # Load GHI ECMWF forecasts
  nc <- nc_open(file.path(forecast_directory, list.files(file.path(forecast_directory), pattern=site)))
  nwp <- ncvar_get(nc, varid="irradiance")
  nc_close(nc)
  
  # Load GHI telemetry and sun_up indicator
  nc <- nc_open(file.path(telemetry_directory, list.files(file.path(telemetry_directory), pattern=glob2rx(paste("*", site, "*2018*", sep="")))))
  tel <- ncvar_get(nc, varid="irradiance", count=c(dim(nwp)[1], -1))
  sun_up <- ncvar_get(nc, varid="sun_up", count=c(dim(nwp)[1], -1))
  nc_close(nc)
  
  sun_up <- apply(sun_up, MARGIN = c(1,2), FUN = as.logical)
  
  # Load GHI telemetry and sun_up indicator from end of 2017 for PeEn forecast
  nc <- nc_open(file.path(telemetry_directory, list.files(file.path(telemetry_directory), pattern=glob2rx(paste("*", site, "*2017*", sep="")))))
  oos_tel <- ncvar_get(nc, varid="irradiance")
  oos_sun_up <- ncvar_get(nc, varid="sun_up")
  nc_close(nc)
  
  # -----------------------------------------------------------------
  # Conduct forecasts and get metrics 
  
  for (method in forecast_names) {
    if (method=="Climatology") {
      fc <- forecast_climatology(tel, percentiles, sun_up)    
    } else if (method=="ECMWF Ensemble") {
      fc <- forecast_NWP(nwp, percentiles, sun_up)    
    } else if (method=="Ch-PeEn") {
      fc <- forecast_Ch_PeEn(tel, percentiles, sun_up)    
    } else if (method=="PeEn") {
      fc <- forecast_PeEn(tel, percentiles, sun_up, num_peen, oos_tel)    
    } else if (method=="Gaussian"){
      fc <- forecast_Gaussian(nwp, tel, percentiles, sun_up)
    } else stop(paste("Forecast method", method, "not recognized"))
  
    metrics_df[site, method, "CRPS"] <- qwCRPS(fc, as.vector(t(tel)), as.vector(t(sun_up)), weighting = "none")
    metrics_df[site, method, "Left-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(tel)), as.vector(t(sun_up)), weighting = "left")
    metrics_df[site, method, "Center weighted CRPS"] <- qwCRPS(fc, as.vector(t(tel)), as.vector(t(sun_up)), weighting = "center")
    metrics_df[site, method, "Right-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(tel)), as.vector(t(sun_up)), weighting = "right")
    
    reliability_df[site, method,] <- reliability(fc, as.vector(t(tel)), as.vector(t(sun_up)), percentiles)
    interval_width_df[site, method,] <- interval_width(fc, as.vector(t(sun_up)), intervals = intervals)$widths
  }
  
  return(list(metrics_df=metrics_df, reliability_df=reliability_df, interval_width_df=interval_width_df))
}

for (site in site_names) {
  results <- get_site_data(site, metrics_df, reliability_df, interval_width_df)
  metrics_df <- results$metrics_df
  reliability_df <- results$reliability_df
  interval_width_df <- results$interval_width_df
}

# -----------------------------------------------------------------
# Export graphs

for (site in site_names) {
  # Reliability plot
  df <- melt(reliability_df[site,,], varnames=c("Method", "levels"))
  
  ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
    xlab("Nominal proportion") + ylab("Observed proportion") + 
    scale_shape_manual(values = c(17, 0, 16, 5, 4)) + 
    geom_line(data = data.frame(x=c(0,1), y=c(0,1)), mapping=aes(x=x,y=y), col="black") + 
    theme(legend.justification=c(1,0), legend.position=c(0.95,0.05), text = ggplot2::element_text(size=14))
  ggsave(file.path(output_directory, paste(site, "reliability.pdf", sep="_")))
  
  # Sharpness plot
  df <- melt(interval_width_df[site,,], varnames=c("Method", "levels"))
  
  ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
    xlab("Central Interval (%)") + ylab(expression(paste("Average Width (W/m"^"2", ")"))) + 
    scale_shape_manual(values = c(17, 0, 16, 5, 4)) + 
    theme(legend.justification=c(0,1), legend.position=c(0.05,0.95), text = ggplot2::element_text(size=14)) + 
    scale_x_continuous(breaks=intervals, labels=intervals*100)
  ggsave(file.path(output_directory, paste(site, "interval_width.pdf", sep="_")))
}

# -----------------------------------------------------------------
# Export metrics table

sink(file.path(output_directory, 'CRPS_results.csv'))
for (metric in dimnames(metrics_df)[[3]]) {
  cat(paste(metric, "\n"))
  write.csv(metrics_df[,,metric])
  cat("\n")
}
sink()
