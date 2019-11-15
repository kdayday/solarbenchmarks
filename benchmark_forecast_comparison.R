# Author: Kate Doubleday
# Last updated: 11/14/2019
# -----------------------------------------------------------------
# Load dependencies

library(here)
library(ncdf4)
library(pracma)
library(truncnorm)
library(ggplot2)
library(ggfan)
library(reshape2)
library(tidyr)
library(tibble)
source("evaluation_functions.R")
source("forecast_methods.R")

# -----------------------------------------------------------------
# Define constants

telemetry_directory <- here("SURFRAD_files", "Yearlong")
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
num_peen <- 20 # Number of members in the hourly persistence ensemble
intrahour_training_hours <- 2 # Number of hours of preceeding data to use for training intra-hour PeEn and intra-hour Gaussian methods
resolution <- c("Hourly", "Intrahour")

forecast_names <- list(Hourly=c("Climatology", "Ch-PeEn", "PeEn", "ECMWF Ensemble", "ECMWF Gaussian"),
                      Intrahour=c("Climatology", "Ch-PeEn", "PeEn", "Smart persistence Gaussian"))

metric_names <- c("CRPS", "Left-tail weighted CRPS", "Center weighted CRPS", "Right-tail weighted CRPS")

# -----------------------------------------------------------------
# Add option to export .R graph objects
Rexport <- F

# -----------------------------------------------------------------
get_site_data <- function(res, site, metrics_df, reliability_df, interval_width_df) {
  # Load GHI ECMWF forecasts
  
  nc <- nc_open(file.path(forecast_directory, list.files(file.path(forecast_directory), pattern=site)))
  nwp <- ncvar_get(nc, varid="irradiance")
  nc_close(nc)
  
  # Load GHI telemetry and sun_up indicator
  nc <- nc_open(file.path(telemetry_directory, res, list.files(file.path(telemetry_directory, res), pattern=glob2rx(paste("*", site, "*2018*", sep="")))))
  GHI <- ncvar_get(nc, varid="irradiance", count=c(dim(nwp)[1], -1))
  sun_up <- ncvar_get(nc, varid="sun_up", count=c(dim(nwp)[1], -1))
  clearsky_GHI <- ncvar_get(nc, varid="clearsky_irradiance", count=c(dim(nwp)[1], -1))
  nc_close(nc)
  
  sun_up <- apply(sun_up, MARGIN = c(1,2), FUN = as.logical)
  
  # Load GHI telemetry and sun_up indicator from end of 2017 for PeEn forecast
  nc <- nc_open(file.path(telemetry_directory, res, list.files(file.path(telemetry_directory, res), pattern=glob2rx(paste("*", site, "*2017*", sep="")))))
  GHI_2017 <- ncvar_get(nc, varid="irradiance")
  sun_up_2017 <- ncvar_get(nc, varid="sun_up")
  nc_close(nc)
  
  # -----------------------------------------------------------------
  # Conduct forecasts and get metrics 
  ts_per_hour <- ncol(GHI)/24
  
  for (method in forecast_names[[res]]) {
    if (method=="Climatology") {
      fc <- forecast_climatology(GHI, percentiles, sun_up)    
    } else if (method=="ECMWF Ensemble") {
      fc <- forecast_NWP(nwp, percentiles, sun_up)    
    } else if (method=="Ch-PeEn") {
      fc <- forecast_Ch_PeEn(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour)    
    } else if (method=="PeEn") {
      if (res == "Hourly") {
        fc <- forecast_PeEn_hourly(GHI, percentiles, sun_up, num_peen, GHI_2017)      
      } else fc <- forecast_PeEn_intrahour(as.vector(t(GHI)), percentiles, as.vector(t(sun_up)), as.vector(t(clearsky_GHI)), ts_per_hour=ts_per_hour, nhours=intrahour_training_hours)
    } else if (method=="ECMWF Gaussian"){
      fc <- forecast_Gaussian_hourly(nwp, GHI, percentiles, sun_up, clearsky_GHI)
    } else if (method=="Smart persistence Gaussian") { 
      fc <- forecast_Gaussian_intrahour(as.vector(t(GHI)), percentiles, as.vector(t(sun_up)), as.vector(t(clearsky_GHI)), ts_per_hour=ts_per_hour, nhours=intrahour_training_hours)
    } else stop(paste("Forecast method", method, "not recognized"))
  
    metrics_df[site, method, "CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "none")
    metrics_df[site, method, "Left-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "left")
    metrics_df[site, method, "Center weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "center")
    metrics_df[site, method, "Right-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "right")
    
    reliability_df[site, method,] <- reliability(fc, as.vector(t(GHI)), as.vector(t(sun_up)), percentiles)
    interval_width_df[site, method,] <- interval_width(fc, as.vector(t(sun_up)), intervals = intervals)$widths
    
    # Export sample forecasts
    g <- plot_fanplot(fc, as.vector(t(tel)), c(1, 72))
    ggsave(file.path(output_directory, paste(site, method, "sample.pdf", sep="_")), height=4, width=6)
    save(g, file=file.path(output_directory, paste(site, method, "sample.R", sep="_")))
  }
  
  return(list(metrics_df=metrics_df, reliability_df=reliability_df, interval_width_df=interval_width_df))
}

for (res in resolution) {
  # Prepare results data frames
  metrics_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(metric_names)), dimnames=list(site_names, forecast_names[[res]], metric_names))
  reliability_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(percentiles)), dimnames=list(site_names, forecast_names[[res]], percentiles))
  interval_width_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(intervals)), dimnames=list(site_names, forecast_names[[res]], intervals))
  
  
  for (site in site_names) {
    results <- get_site_data(res, site, metrics_df, reliability_df, interval_width_df)
    metrics_df <- results$metrics_df
    reliability_df <- results$reliability_df
    interval_width_df <- results$interval_width_df
  
    # Export graphs
    # ------------------------------------------------------------------
    shapes <- c(17, 0, 16, 5, 4)
    
    # Reliability plot
    df <- melt(reliability_df[site,,], varnames=c("Method", "levels"))
    
    g <- ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
      xlab("Nominal proportion") + ylab("Observed proportion") + 
      scale_shape_manual(values = shapes[1:nrow(reliability_df)]) + 
      geom_line(data = data.frame(x=c(0,1), y=c(0,1)), mapping=aes(x=x,y=y), col="black") + 
      theme(legend.justification=c(1,0), legend.position=c(0.98,0.02), text = ggplot2::element_text(size=14))
    ggsave(file.path(output_directory, paste(site, res, "reliability.pdf", sep="_")), height=4, width=8)
    if (Rexport) {
      save(g, file=file.path(output_directory, paste(site, res, "reliability_plot.R", sep="_")))
    }
    
    # Sharpness plot
    df <- melt(interval_width_df[site,,], varnames=c("Method", "levels"))
    
    g <- ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
      xlab("Central Interval (%)") + ylab(expression(paste("Average Width (W/m"^"2", ")"))) + 
      scale_shape_manual(values = shapes[1:nrow(interval_width_df)]) + 
      theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), text = ggplot2::element_text(size=14)) + 
      scale_x_continuous(breaks=intervals, labels=intervals*100)
    ggsave(file.path(output_directory, paste(site, res, "interval_width.pdf", sep="_")), height=4, width=8)
    if (Rexport) {
      save(g, file=file.path(output_directory, paste(site, res, "interval_width_plot.R", sep="_")))  
    }
  }

  # -----------------------------------------------------------------
  # Export metrics table
  
  sink(file.path(output_directory, paste('CRPS', res, 'results.csv', sep="_")))
  for (metric in dimnames(metrics_df)[[3]]) {
    cat(paste(metric, "\n"))
    write.csv(metrics_df[,,metric])
    cat("\n")
  }
  sink()
}
