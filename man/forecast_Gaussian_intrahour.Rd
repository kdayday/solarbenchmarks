% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forecast_methods.R
\name{forecast_Gaussian_intrahour}
\alias{forecast_Gaussian_intrahour}
\title{Do intra-hourly Gaussian error distribution forecast}
\usage{
forecast_Gaussian_intrahour(GHI, percentiles, sun_up, clearsky_GHI,
  ts_per_hour, nhours)
}
\arguments{
\item{GHI}{A vector of the telemetry}

\item{percentiles}{A vector of the percentiles corresponding to the desired
forecast quantiles}

\item{sun_up}{A vector of logicals, indicating whether the sun is up}

\item{clearsky_GHI}{a vector of clear-sky irradiance estimates}

\item{ts_per_hour}{Time-steps per hour, e.g., 12 for a 5-minute resolution
forecast}

\item{nhours}{Number of preceeding hours to collect training errors from,
e.g., 1 or 2}
}
\value{
a matrix of quantile forecasts at each valid time in the input data
}
\description{
Generates an intra-hourly Gaussian error distribution forecast using a doubly
truncated Gaussian distribution. The distribution's mean is a smart
persistence forecast, based on the most recent clear-sky index before each
hourly issue time. The distribution's standard deviation is calculated as the
standard deviation of the smart persistence errors over the past few hours.
}
\details{
Valid for intra-hourly forecasts. For hourly-resolution forecasts, see
\code{\link{forecast_Gaussian_hourly}}. Distribution is truncated at 0 on the
low end and clear-sky GHI on the upper end. While training data is being
collected at the beginning of the day, forecast starts as deterministic
clear-sky forecast and gathers forecast errors starting after the 1st hour.
This function is for hindcasting only, and uses the same clear-sky GHI
estimates for both the historical and forecast values.
}
\seealso{
Other forecast functions: \code{\link{forecast_CH_PeEn}},
  \code{\link{forecast_Gaussian_hourly}},
  \code{\link{forecast_NWP}},
  \code{\link{forecast_PeEn_hourly}},
  \code{\link{forecast_PeEn_intrahour}},
  \code{\link{forecast_climatology}},
  \code{\link{forecast_mcm}}
}
\concept{forecast functions}
