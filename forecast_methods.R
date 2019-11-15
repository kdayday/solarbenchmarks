
#' Climatology method using in-sample data
#' 
#' @param GHI A [day x hour] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x hour] matrix of logicals, indicating whether the sun is up
forecast_climatology <- function(GHI, percentiles, sun_up) {
  xseq <- stats::quantile(as.vector(GHI[sun_up]), probs=percentiles, names=F, na.rm=T, type=1)
  fc <- t(sapply(as.vector(t(sun_up)), FUN=function(s) if (isTRUE(s)) xseq else rep(0, times=length(percentiles)), simplify="array"))
  colnames(fc) <- percentiles
  return(fc)
}

#' Raw numerical weather prediction ensemble method
#' 
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble forecasts
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x hour] matrix of logicals, indicating whether the sun is up
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

#' Persistence ensemble method: hourly
#' 
#' @param GHI A [day x time] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
#' @param num_peen  Number of persistence ensemble members to take
#' @param lead_up_GHI A [day x time] matrix of out-of-sample telemetry for the days leading up to the start of the sample period. Number of days must be >= num_peen
forecast_PeEn_hourly <- function(GHI, percentiles, sun_up, num_peen, lead_up_GHI) {
  if (dim(lead_up_GHI)[1] < num_peen) stop(paste("Must have at least", num_peen, "days of data leading up to the sample people available. Given", dim(lead_up_GHI)[1]))
  
  all_GHI <- rbind(lead_up_GHI, GHI)
  
  fc <- sapply(seq(nrow(lead_up_GHI)+1, nrow(all_GHI)), FUN=function(day) {
    sapply(seq_len(ncol(GHI)), FUN=function(hr) {if (isTRUE(sun_up[day-nrow(lead_up_GHI), hr])) {
      stats::quantile(all_GHI[(day-num_peen):(day-1),hr], probs=percentiles, names=F, na.rm=T, type=1)} else 
        rep(0, times=length(percentiles))})
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))
  
  colnames(fc) <- percentiles
  return(fc)
}

#' Persistence ensemble method: intra-hour
#' While training data is being collected at the beginning of the day, forecast starts as deterministic clear-sky forecast and gathers CSI's in 2nd hour. 
#' 
#' @param GHI A vector of telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution forecast
#' @param nhours Number of preceeding hours to collect training errors from
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

#' Complete-history persistence ensemble method from aggregating clear-sky indices on an hourly basis
#' 
#' @param GHI A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
#' @param clearsky_GHI a [day x time] matrix of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution forecast
forecast_Ch_PeEn <- function(GHI, percentiles, sun_up, clearsky_GHI, ts_per_hour) {
  
  # Clear sky indices
  CSI <- GHI/clearsky_GHI
  fc <- sapply(seq_len(nrow(GHI)), FUN=function(day) {
    sapply(seq_len(ncol(GHI)), FUN=function(t) {if (isTRUE(sun_up[day, t])) {
      # For both hourly and subhourly forecasts, CSI's are grouped hourly. First find the subset of time columns for this hour:
      cols <- ((ceil(t/ts_per_hour)-1)*ts_per_hour)+1:ts_per_hour
      stats::quantile(as.vector( CSI[,cols][sun_up[,cols]])*clearsky_GHI[day,t] , probs=percentiles, names=F, na.rm=T, type=1)} else 
      rep(0, times=length(percentiles))})
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))

  colnames(fc) <- percentiles
  return(fc)
}

#' Gaussian error dressing: hourly, using hourly errors of the deterministic NWP forecast to fit truncated normal distributions
#' Distribution is truncated at 0 on the low end and clear-sky GHI on the upper end
#' 
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble forecasts; First member is treated as the control member.
#' @param GHI A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
#' @param clearsky_GHI a [day x time] matrix of clear-sky irradiance estimates
forecast_Gaussian_hourly <- function(nwp, GHI, percentiles, sun_up, clearsky_GHI) {
  if ((dim(nwp)[3]/dim(nwp)[2])%%1!=0 | (dim(nwp)[3]/dim(nwp)[2])==1) stop("Unknown handling for given number of issue and lead times")
  if (any(dim(GHI)!=dim(sun_up))) stop("Given incompatible number of forecasts and sun_up times")
  
  n <- ncol(GHI)
  nwp_ctrl <- matrix(aperm(nwp[,,1:(dim(nwp)[3]/dim(nwp)[2]),1], c(3,2,1)), nrow=nrow(nwp), byrow=T) # Reshape to [day x hour]
  # Fit a standard deviation for each hour of the day, based on all residuals from that hour
  daily_sd <- sapply(1:n, FUN=function(i) sd(nwp_ctrl[,i] - GHI[,i], na.rm=T))

  # Transform from matrices to vectors for simplicity of for loop
  sun_up_vector <- as.vector(t(sun_up))
  ctrl_vector <- as.vector(t(nwp_ctrl))
  clearsky_vector <- as.vector(t(clearsky_GHI))
  
  fc <- t(sapply(seq_along(sun_up_vector), FUN=function(i) {if (isTRUE(sun_up_vector[i])) qtruncnorm(p=percentiles, a=0, b=clearsky_vector[i], mean=ctrl_vector[i], sd=daily_sd[(i-1)%%n+1])
    else rep(0, times=length(percentiles))}))
  colnames(fc) <- percentiles
  return(fc)
}

#' Gaussian error dressing: intra-hour, using smart persistence errors from past nhours hours
#'
#' @param GHI A vector of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param ts_per_hour Time-steps per hour, e.g., 12 for a 5-minute resolution forecast
#' @param nhours Number of preceeding hours to collect training errors from
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
    std_dev <- sd(error_subset, na.rm=T)
    sapply(seq_len(ts_per_hour), FUN=function(j, std_dev) {
      if (isTRUE(sun_up[i+j-1])){
      return(qtruncnorm(p=percentiles, a=0, b=clearsky_GHI[i+j-1], mean=hourly_smart_persistence[i+j-1], sd=std_dev))
    } else return(rep(0, times=length(percentiles)))}, std_dev=std_dev, simplify="array")
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(GHI), length(percentiles)))
  
  colnames(fc) <- percentiles
  return(fc)
}
