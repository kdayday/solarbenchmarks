# Author: Kate Doubleday
# Last updated: 5/13/2020
# -----------------------------------------------------------------
# Execute this script to generate a Results/ folder containing 200+ plots comparing the benchmarks
# for each of the 7 SURFRAD sites, .csv files containing the forecast quantiles at each
# time-step through the year 2018, and .csv files comparing the CRPS summary statistics
# for each site across the methods at each temporal scale. If the ECMWF data is not available,
# the two ECMWF-based forecasts will be skipped. 
# -----------------------------------------------------------------
# Load dependencies

library(solarbenchmarks)

# -----------------------------------------------------------------

print("Benchmark beginning, and will run for a while. See Results directory and its sub-directories.")

# -----------------------------------------------------------------
# Define constants

telemetry_directory <- here::here("GHI_files")
forecast_directory <- here::here("ECMWF_files", "SURFRAD_sites")
output_directory <- here::here("Results")
dir.create(output_directory, showWarnings = FALSE)
pit_directory <- here::here("Results", "PIT_histograms")
dir.create(pit_directory, showWarnings = FALSE)
ts_directory <- here::here("Results", "Quantile_time_series")
dir.create(ts_directory, showWarnings = FALSE)
qs_directory <- here::here("Results", "Quantile_score_plots")
dir.create(qs_directory, showWarnings = FALSE)

site_names <- c("Bondville",
                "Boulder",
                "Desert_Rock",
                "Fort_Peck",
                "Goodwin_Creek",
                "Penn_State",
                "Sioux_Falls")

percentiles <- seq(0.01, 0.99, by=0.01)
intervals <- seq(0.1, 0.9, by=0.1) # Central intervals for sharpness evaluation
num_peen_days <- 20 # Number of members in the hourly persistence ensemble
intrahour_training_hours <- 2 # Number of hours of preceeding data to use for training intra-hour PeEn and intra-hour Gaussian methods
mcm_training_days <- 20 # Number of days of preceeding data to use for training intra-hour MCM method
histogram_bins <- 20 # Number of bins to use in PIT histogram

resolution <- c("Hourly", "Intrahour")

# Choose the appropriate subset of forecast methods, based on temporal resolution and ECMWF data availability (or not).
if (all(sapply(site_names, FUN=function(site) paste(site, ".nc", sep="")) %in% list.files(file.path(forecast_directory)))) {
  forecast_names <- list(Hourly=c("Climatology", "Ch-PeEn", "PeEn", "ECMWF Ensemble", "ECMWF Gaussian"))
} else forecast_names <- list(Hourly=c("Climatology", "Ch-PeEn", "PeEn"))
forecast_names$Intrahour<- c("Climatology", "Ch-PeEn", "PeEn", "Smart persistence Gaussian", "MCM")
  
metric_names <- c("CRPS", "Left-tail weighted CRPS", "Center weighted CRPS", "Right-tail weighted CRPS")

# -----------------------------------------------------------------
# Data export options
# Option to export a .csv of the forecast quantiles
csv_export <- T

# Option to export .R graph objects
R_graph_export <- F

# -----------------------------------------------------------------
get_site_data <- function(res, site, metrics_df, reliability_df, interval_width_df) {
  
  # Load GHI ECMWF forecasts
  if (file.exists(file.path(forecast_directory, paste(site, ".nc", sep="")))) {
    nc <- ncdf4::nc_open(file.path(forecast_directory, paste(site, ".nc", sep="")))
    nwp <- ncdf4::ncvar_get(nc, varid="irradiance")
    ncdf4::nc_close(nc)
    ndays <- nrow(nwp) # For loading in the same number of days of SURFRAD data, starting Jan. 1st. Useful during testing
  } else{ # If data is unavailable
    ndays <- 365
  }
  
  # Load GHI telemetry and sun_up indicator
  nc <- ncdf4::nc_open(file.path(telemetry_directory, res, list.files(file.path(telemetry_directory, res), pattern=glob2rx(paste("*", site, "*2018*", sep="")))))
  GHI <- ncdf4::ncvar_get(nc, varid="irradiance", count=c(ndays, -1))
  sun_up <- ncdf4::ncvar_get(nc, varid="sun_up", count=c(ndays, -1))
  clearsky_GHI <- ncdf4::ncvar_get(nc, varid="clearsky_irradiance", count=c(ndays, -1))
  ncdf4::nc_close(nc)
  
  sun_up <- apply(sun_up, MARGIN = c(1,2), FUN = as.logical)
  
  # Load GHI telemetry and sun_up indicator from end of 2017 for PeEn forecast
  nc <- ncdf4::nc_open(file.path(telemetry_directory, res, list.files(file.path(telemetry_directory, res), pattern=glob2rx(paste("*", site, "*2017*", sep="")))))
  GHI_2017 <- ncdf4::ncvar_get(nc, varid="irradiance")
  clearsky_GHI_2017 <- ncdf4::ncvar_get(nc, varid="clearsky_irradiance")
  sun_up_2017 <- ncdf4::ncvar_get(nc, varid="sun_up")
  ncdf4::nc_close(nc)
  sun_up_2017 <- apply(sun_up_2017, MARGIN = c(1,2), FUN = as.logical)
  
  # -----------------------------------------------------------------
  # Conduct forecasts and get metrics 
  ts_per_hour <- ncol(GHI)/24
  
  for (method in forecast_names[[res]]) {
    if (method=="Climatology") {
      fc <- forecast_climatology(GHI, percentiles, sun_up)    
    } else if (method=="ECMWF Ensemble") {
      fc <- forecast_NWP(nwp, percentiles, sun_up)    
    } else if (method=="Ch-PeEn") {
      fc <- forecast_CH_PeEn(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour)    
    } else if (method=="PeEn") {
      if (res == "Hourly") {
        fc <- forecast_PeEn_hourly(GHI, percentiles, sun_up, num_peen_days, GHI_2017)      
      } else fc <- forecast_PeEn_intrahour(as.vector(t(GHI)), percentiles, as.vector(t(sun_up)), 
                                           as.vector(t(clearsky_GHI)), ts_per_hour=ts_per_hour, nhours=intrahour_training_hours)
    } else if (method=="ECMWF Gaussian"){
      fc <- forecast_Gaussian_hourly(nwp, GHI, percentiles, sun_up, clearsky_GHI)
    } else if (method=="Smart persistence Gaussian") { 
      fc <- forecast_Gaussian_intrahour(as.vector(t(GHI)), percentiles, as.vector(t(sun_up)), 
                                        as.vector(t(clearsky_GHI)), ts_per_hour=ts_per_hour, nhours=intrahour_training_hours)
    } else if (method=="MCM") { 
      fc <- forecast_mcm(as.vector(t(GHI)), as.vector(t(GHI_2017)), 
                         percentiles, as.vector(t(sun_up)), as.vector(t(sun_up_2017)), 
                         as.vector(t(clearsky_GHI)), as.vector(t(clearsky_GHI_2017)), 
                         ts_per_hour, mcm_training_days)
    } else stop(paste("Forecast method", method, "not recognized"))
  
    metrics_df[site, method, "CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "none")
    metrics_df[site, method, "Left-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "left")
    metrics_df[site, method, "Center weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "center")
    metrics_df[site, method, "Right-tail weighted CRPS"] <- qwCRPS(fc, as.vector(t(GHI)), as.vector(t(sun_up)), weighting = "right")
    
    reliability_df[site, method,] <- reliability(fc, as.vector(t(GHI)), as.vector(t(sun_up)))
    interval_width_df[site, method,] <- interval_width(fc, as.vector(t(sun_up)), intervals = intervals)$widths
    
    # If desired, export this forecast matrix of quantiles to a .csv file for future use.
    if (csv_export) {
      write.table(fc, file=file.path(ts_directory, paste(site, res, method, "quantiles.csv", sep=" ")), sep=",", row.names=F)
    }
    
    # Plot PIT histogram
    plot_PIT_histogram(fc, as.vector(t(GHI)), as.vector(t(sun_up)), histogram_bins, site, res, method, pit_directory, R_graph_export)
    # Plot unweighted and weighted quantile score comparisons
    plot_quantile_score(fc, as.vector(t(GHI)), as.vector(t(sun_up)), percentiles, site, res, method, qs_directory, R_graph_export)
    
    # Export sample forecasts: Days 133-135, May 13-15th
    window <- c(132*24*ts_per_hour+1 + 6*ts_per_hour, 135*24*ts_per_hour + 6*ts_per_hour)
    g <- plot_fanplot(fc, as.vector(t(GHI)), window, ts_per_hour, output_directory, site, res, method, R_graph_export)
  }
  return(list(metrics_df=metrics_df, reliability_df=reliability_df, interval_width_df=interval_width_df))
}

# Repeat the analysis for each temporal resolution: Hourly forecast or intra-hourly forecast
for (res in resolution) {
  # Prepare results data frames
  metrics_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(metric_names)), dimnames=list(site_names, forecast_names[[res]], metric_names))
  reliability_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(percentiles)), dimnames=list(site_names, forecast_names[[res]], percentiles))
  interval_width_df <- array(dim=c(length(site_names), length(forecast_names[[res]]), length(intervals)), dimnames=list(site_names, forecast_names[[res]], intervals))
  
  # Repeat the analysis for each of the 7 SURFRAD sites
  for (site in site_names) {
    # Conduct the forecast validation for this site, saving the results for each method to the three data frames
    results <- get_site_data(res, site, metrics_df, reliability_df, interval_width_df)
    metrics_df <- results$metrics_df
    reliability_df <- results$reliability_df
    interval_width_df <- results$interval_width_df
  
    # Export reliability and interval width graphs comparing the methods for each site
    # ------------------------------------------------------------------
    shapes <- c(17, 0, 16, 5, 4) # Select point shapes for graphs
    plot_reliability(reliability_df, site, res, shapes, output_directory, R_graph_export)
    plot_interval_width(interval_width_df, site, res, shapes, output_directory, R_graph_export)
  }

  # -----------------------------------------------------------------
  # Export CRPS metrics table, comparing the methods across the sites
  sink(file.path(output_directory, paste('CRPS', res, 'results.csv', sep="_")))
  for (metric in dimnames(metrics_df)[[3]]) {
    cat(paste(metric, "\n"))
    write.csv(metrics_df[,,metric])
    cat("\n")
  }
  sink()
}
