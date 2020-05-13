#' Do climatology forecast
#'
#' Generates a static climatology forecast using all GHI measurements in the
#' data set, transformed to a full probabilistic forecast using an empirical
#' CDF.
#'
#' Valid for intra-hourly and hourly forecast resolutions
#' @family forecast functions
#'
#' @param GHI A [day x time-step] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A [day x time-step] matrix of logicals, indicating whether the
#'   sun is up
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_climatology <- function(GHI, percentiles, sun_up) {
  xseq <- stats::quantile(as.vector(GHI[sun_up]), probs=percentiles, names=F, na.rm=T, type=1)
  fc <- t(sapply(as.vector(t(sun_up)), FUN=function(s) if (isTRUE(s)) xseq else rep(0, times=length(percentiles)), simplify="array"))
  colnames(fc) <- percentiles
  return(fc)
}

#' Do raw numerical weather prediction ensemble forecast
#'
#' Generates a probabilistic forecast from a NwP ensemble by applying an
#' empirical CDF.
#'
#' Valid at the resolution of the NWP ensemble (e.g., hourly)
#' @family forecast functions
#'
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble
#'   forecasts
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A [day x hour] matrix of logicals, indicating whether the sun
#'   is up
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_NWP <- function(nwp, percentiles, sun_up) {
  if ((dim(nwp)[3]/dim(nwp)[2])%%1!=0 | (dim(nwp)[3]/dim(nwp)[2])==1) stop("Unknown handling for number of issue and lead times")
  # Simplify to most recent available forecast
  nwp_rolling <- matrix(aperm(nwp[,,1:(dim(nwp)[3]/dim(nwp)[2]),], c(3,2,1,4)), ncol=dim(nwp)[4])
  if (nrow(nwp_rolling) !=length(sun_up)) stop("Given incompatible number of forecasts and sun_up times")
  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {if (isTRUE(as.vector(t(sun_up))[i])) (stats::quantile(nwp_rolling[i,], probs=percentiles, names=F, na.rm=T, type=1)) 
    else rep(0, times=length(percentiles))}, simplify="array"))
  colnames(fc) <- percentiles
  return(fc)
}

#' Do hourly persistence ensemble forecast
#'
#' Generates an hourly persistence ensemble forecast, using GHI measurements
#' from the same hour-of-day over the last few days. The GHI ensemble is
#' transformed to a full probabilistic forecast using an empirical CDF.
#'
#' Valid for hourly resolution forecasts. For intra-hourly forecasts, see
#' \code{\link{forecast_PeEn_intrahour}}.
#' @family forecast functions
#'
#' @param GHI A [day x time] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun
#'   is up
#' @param num_days  Number of days of data to take, indicating size of
#'   persistence ensemble
#' @param lead_up_GHI A [day x time] matrix of out-of-sample telemetry for the
#'   days leading up to the start of the sample period. Number of days must be
#'   >= num_days
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_PeEn_hourly <- function(GHI, percentiles, sun_up, num_days, lead_up_GHI) {
  if (dim(lead_up_GHI)[1] < num_days) stop(paste("Must have at least", num_days, "days of data leading up to the sample people available. Given", dim(lead_up_GHI)[1]))
  
  all_GHI <- rbind(lead_up_GHI, GHI)
  
  fc <- sapply(seq(nrow(lead_up_GHI)+1, nrow(all_GHI)), FUN=function(day) {
    sapply(seq_len(ncol(GHI)), FUN=function(hr) {if (isTRUE(sun_up[day-nrow(lead_up_GHI), hr])) {
      stats::quantile(all_GHI[(day-num_days):(day-1),hr], probs=percentiles, names=F, na.rm=T, type=1)} else 
        rep(0, times=length(percentiles))})
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))
  
  colnames(fc) <- percentiles
  return(fc)
}

#' Do intra-hourly persistence ensemble forecast
#'
#' Generates an intra-hourly persistence ensemble forecast, using a few recent
#' hours to generate a clear-sky index ensemble. The clear-sky index ensemble is
#' multiplied by the upcoming clear-sky GHI to generate a GHI ensemble,
#' transformed to a full probabilistic forecast using an empirical CDF.
#'
#' Valid for intra-hourly forecasts. For hourly resolution forecasts, see
#' \code{\link{forecast_PeEn_hourly}}. While training data is being collected at
#' the beginning of the day, forecast starts as deterministic clear-sky forecast
#' and gathers CSI's for the ensemble starting after the 1st hour. This function is for hindcasting only, and
#' uses the same clear-sky GHI estimates for both the historical and forecast
#' values.
#' @family forecast functions
#'
#' @param GHI A vector of telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution
#'   forecast
#' @param nhours Number of preceeding hours to collect training errors from,
#'   e.g., 1 or 2
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_PeEn_intrahour <- function(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour, nhours) {
  
  # Forecast is updated hourly; find indices of start of each hour
  update_times <-seq(from=1, to=length(GHI), by=ts_per_hour) 
  
  fc <- sapply(update_times, FUN=function(i) {
    # On an hourly basis, collect last hour's ensemble of CSI's
    CSI_subset <- (GHI[max(i-nhours*ts_per_hour, 0):max(i-1, 0)]/clearsky_GHI[max(i-nhours*ts_per_hour, 0):max(i-1, 0)])[sun_up[max(i-nhours*ts_per_hour, 0):max(i-1, 0)]]
    # If there is no training data available, either because at start of vector or start of day, assume clear sky. Ignore division by 0's. 
    if (!any(is.finite(CSI_subset))) CSI_subset <- 1 # Force deterministic if no valid CSI's available yet
    sapply(seq_len(ts_per_hour), FUN=function(j, CSI_subset) {
      if (isTRUE(sun_up[i+j-1])){
        return(stats::quantile(CSI_subset[is.finite(CSI_subset)]*clearsky_GHI[i+j-1], probs=percentiles, names=F, na.rm=T, type=1))
      } else return(rep(0, times=length(percentiles)))}, CSI_subset=CSI_subset, simplify="array")
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))

  colnames(fc) <- percentiles
  return(fc)
}

#' Do complete-history persistence ensemble forecast
#'
#' Generates a CH-PeEn forecast using all clear-sky indices at the same hour of
#' day, for both intra-hourly and hourly-resolution forecasts. The clear-sky
#' index ensemble is multiplied by the forecast clear-sky GHI to get a GHI
#' ensemble, which is transformed to a full probabilistic forecast through an
#' empirical CDF.
#'
#' Valid for both intra-hourly and hourly-resolution forecasts. This function is
#' for hindcasting only, and uses the same clear-sky GHI estimates for both the
#' historical and forecast values.
#' @family forecast functions
#'
#' @param GHI A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun
#'   is up
#' @param clearsky_GHI a [day x time] matrix of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution
#'   forecast
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_CH_PeEn <- function(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour) {
  
  # Clear sky indices
  CSI <- GHI/clearsky_GHI
  fc <- sapply(seq_len(nrow(GHI)), FUN=function(day) {
    sapply(seq_len(ncol(GHI)), FUN=function(t) {if (isTRUE(sun_up[day, t])) {
      # For both hourly and subhourly forecasts, CSI's are grouped hourly. First find the subset of time columns for this hour:
      cols <- ((pracma::ceil(t/ts_per_hour)-1)*ts_per_hour)+1:ts_per_hour
      stats::quantile(as.vector( CSI[,cols][sun_up[,cols]])*clearsky_GHI[day,t] , probs=percentiles, names=F, na.rm=T, type=1)} else 
      rep(0, times=length(percentiles))})
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))

  colnames(fc) <- percentiles
  return(fc)
}

#' Do hourly-resolution Gaussian error distribution forecast
#'
#' Generates an hourly Gaussian error distribution forecast using a doubly
#' truncated Gaussian distribution. The distribution's mean is the deterministic
#' forecast from the control member (1st member) in the NWP ensemble. The
#' distribution's standard deviation is calculated as the standard deviation of
#' the deterministic forecast errors at the same hour-of-day over the sample
#' set.
#'
#' Valid for hourly resolution forecasts. For intra-hourly forecasts, see
#' \code{\link{forecast_Gaussian_intrahour}}. Distribution is truncated at 0 on
#' the low end and clear-sky GHI on the upper end.
#'
#' @family forecast functions
#' 
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble forecasts; First member is treated as the control member.
#' @param GHI A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
#' @param clearsky_GHI a [day x time] matrix of clear-sky irradiance estimates
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_Gaussian_hourly <- function(nwp, GHI, percentiles, sun_up, clearsky_GHI) {
  if ((dim(nwp)[3]/dim(nwp)[2])%%1!=0 | (dim(nwp)[3]/dim(nwp)[2])==1) stop("Unknown handling for given number of issue and lead times")
  if (any(dim(GHI)!=dim(sun_up))) stop("Given incompatible number of forecasts and sun_up times")
  
  n <- ncol(GHI)
  nwp_ctrl <- matrix(aperm(nwp[,,1:(dim(nwp)[3]/dim(nwp)[2]),1], c(3,2,1)), nrow=nrow(nwp), byrow=T) # Reshape to [day x hour]
  # Fit a standard deviation for each hour of the day, based on all residuals from that hour
  daily_sd <- sapply(1:n, FUN=function(i) stats::sd(nwp_ctrl[,i] - GHI[,i], na.rm=T))

  # Transform from matrices to vectors for simplicity of for loop
  sun_up_vector <- as.vector(t(sun_up))
  ctrl_vector <- as.vector(t(nwp_ctrl))
  clearsky_vector <- as.vector(t(clearsky_GHI))
  
  fc <- t(sapply(seq_along(sun_up_vector), 
                 FUN=function(i) {if (isTRUE(sun_up_vector[i])) truncnorm::qtruncnorm(p=percentiles, a=0, b=clearsky_vector[i], mean=ctrl_vector[i], sd=daily_sd[(i-1)%%n+1])
    else rep(0, times=length(percentiles))}))
  colnames(fc) <- percentiles
  return(fc)
}

#' Do intra-hourly Gaussian error distribution forecast
#'
#' Generates an intra-hourly Gaussian error distribution forecast using a doubly
#' truncated Gaussian distribution. The distribution's mean is a smart
#' persistence forecast, based on the most recent clear-sky index before each
#' hourly issue time. The distribution's standard deviation is calculated as the
#' standard deviation of the smart persistence errors over the past few hours.
#'
#' Valid for intra-hourly forecasts. For hourly-resolution forecasts, see
#' \code{\link{forecast_Gaussian_hourly}}. Distribution is truncated at 0 on the
#' low end and clear-sky GHI on the upper end. While training data is being
#' collected at the beginning of the day, forecast starts as deterministic
#' clear-sky forecast and gathers forecast errors starting after the 1st hour.
#' This function is for hindcasting only, and uses the same clear-sky GHI
#' estimates for both the historical and forecast values.
#'
#' @family forecast functions
#'
#' @param GHI A vector of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution
#'   forecast
#' @param nhours Number of preceeding hours to collect training errors from,
#'   e.g., 1 or 2
#' @return a matrix of quantile forecasts at each valid time in the input data
#' @export
forecast_Gaussian_intrahour <- function(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour, nhours) {
  
  # First, calculate deterministic smart persistence on an hourly basis:
  # Last measurement of the preceeding hour is persisted for the next hour
  update_times <-seq(from=1, to=length(GHI), by=ts_per_hour) 
  hourly_CSI <- c(1, GHI[(update_times-1)]/clearsky_GHI[(update_times-1)]) # assume first value is clear sky
  # If no CSI is available (i.e., night-time), assume clear sky:
  hourly_CSI[!is.finite(hourly_CSI)] <- 1
  # Persist CSI on an hourly basis over 5-minute clear-sky estimates
  hourly_smart_persistence <- as.vector(sapply(seq_along(update_times), FUN=function(i) {hourly_CSI[i]*clearsky_GHI[update_times[i]:(update_times[i] + ts_per_hour-1)]}))
  # Calculate smart persistence errors
  errors <- hourly_smart_persistence - GHI
  
  fc <- sapply(update_times, FUN=function(i) {
    # On an hourly basis, collect last hour's errors and calculate a new standard deviation for the next hour
    error_subset <- errors[max(i-nhours*ts_per_hour, 0):max(i-1, 0)][sun_up[max(i-nhours*ts_per_hour, 0):max(i-1, 0)]]
    if (!any(is.finite(error_subset))) error_subset <- rep(0, times=nhours*ts_per_hour) # Force deterministic if no errors have accumulated yet
    std_dev <- stats::sd(error_subset, na.rm=T)
    sapply(seq_len(ts_per_hour), FUN=function(j, std_dev) {
      if (isTRUE(sun_up[i+j-1])){
      return(truncnorm::qtruncnorm(p=percentiles, a=0, b=clearsky_GHI[i+j-1], mean=hourly_smart_persistence[i+j-1], sd=std_dev))
    } else return(rep(0, times=length(percentiles)))}, std_dev=std_dev, simplify="array")
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))
  
  colnames(fc) <- percentiles
  return(fc)
}


#' Do intra-hourly Markov-chain mixture distribution forecast
#'
#' Valid for intra-hourly forecasts. This generates a clear-sky index transition
#' matrix based on the previous 20 days of data before each hourly issue time.
#' Using that transition matrix, clear-sky index is forecasted over the next D
#' time steps (i.e., up to 12 for a 5-minute resolution forecast.) Clear-sky
#' indices above 2 are treated as outliers and are set equal to 1; <1% of times
#' have CSI > 2. This removes outliers that can reduce resolution in the area
#' with the bulk of the probability.
#'
#' Modified from the main.R script at
#' https://github.com/SheperoMah/MCM-distribution-forecasting, Reported in: J.
#' Munkhammar, J. Widén, D. W. van der Meer, Probabilistic forecasting of
#' high-resolution clear-sky index time-series using a Markov-chain mixture
#' distribution model, Solar Energy vol. 184, pp. 688-695, 2019.
#'
#' @family forecast functions
#'
#' @param GHI A vector of the telemetry
#' @param lead_up_GHI A vector of out-of-sample telemetry for the days leading
#'   up to the start of the sample period. Number of days must be >= num_days
#' @param percentiles A vector of the percentiles corresponding to the desired
#'   forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param lead_up_clearsky_GHI A vector of out-of-sample clear-sky irradiance
#'   estimates for the days leading up to the start of the sample period,
#'   corresponding to lead_up_GHI
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution
#'   forecast
#' @param num_days  Number of days of training data
#' @param numBins (optional) Number of bins to use in the MCM matrix M. Defaults
#'   to number of percentiles + 1.
#' @param numSamples (optional) Number of samples to take from MCM model to
#'   generate empirical CDF. Defaults to 1000.
#' @param h_per_day (optional) Hours per day = 24 (useful for testing)
#' @export
forecast_mcm <- function(GHI, lead_up_GHI, percentiles, sun_up, clearsky_GHI, 
                         lead_up_clearsky_GHI, ts_per_hour, num_days, numBins=length(percentiles)+1,
                         numSamples=1000, h_per_day=24) {
  if (length(lead_up_GHI)/(h_per_day*ts_per_hour) < num_days) stop(paste("Must have at least", num_days, 
                                                 "days of data leading up to the forecast period. Given", length(lead_up_GHI)))
  
  all_GHI <- c(lead_up_GHI, GHI)
  all_clearsky <- c(lead_up_clearsky_GHI, clearsky_GHI)
  CSI <- all_GHI/all_clearsky
  # If no CSI is available (i.e., night-time), assume clear sky
  # Remove outlier CSI's -- force to CSI=1
  CSI[!is.finite(CSI) | CSI > 2] <- 1
  
  update_times <- seq(length(lead_up_GHI)+1, length(all_GHI), by=ts_per_hour)
  fc <- sapply(update_times, FUN=function(i) {
    training_CSI <- CSI[(i-num_days*h_per_day*ts_per_hour):(i-1)]
    sapply(seq_len(ts_per_hour), FUN=function(j, training_CSI) {
      if (isTRUE(sun_up[i-length(lead_up_GHI)+j-1])){
        p <- mcmFit(training_CSI, numBins, numStepAhead=j)
        forecastToolsList  <- mcmForecast(p, min(training_CSI), max(training_CSI), CSI[i-1])
        # Sample, then estimate uniform percentiles through empirical CDF
        forecastSamples <- mcmRnd(forecastToolsList$binStartingValues, forecastToolsList$transitionProbs, numSamples)
        xseq <- stats::quantile(forecastSamples, probs=percentiles, names=F, na.rm=T, type=1)*all_clearsky[i+j-1]
        return(xseq)
      } else return(rep(0, times=length(percentiles)))}, training_CSI=training_CSI, simplify="array")
  }, simplify="array")

  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))
  colnames(fc) <- percentiles
  return(fc)
}
