
#' Climatology method using in-sample data
#' 
#' @param tel A [day x hour] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x hour] matrix of logicals, indicating whether the sun is up
forecast_climatology <- function(tel, percentiles, sun_up) {
  xseq <- stats::quantile(as.vector(tel[sun_up]), probs=percentiles, names=F, na.rm=T, type=1)
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
#' While training data is being collected at the beginning of the day, forecast starts as deterministic clear-sky forecast and evolves as errors appear. 
#' 
#' @param GHI A vector of telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param num_peen  Number of persistence ensemble members to take
forecast_PeEn_minute <- function(GHI, percentiles, sun_up, clearsky_GHI, num_peen) {
  
  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {
    if (isTRUE(sun_up[i])) {
      CSI <- GHI[max(c(i-num_peen,0)):(i-1)]/clearsky_GHI[max(c(i-num_peen,0)):(i-1)] # Get ensemble of clear-sky indices
      # If there is no training data available, either because at start of vector or start of day, assume clear sky. Ignore division by 0's. 
      if (!any(is.finite(CSI))) CSI <- 1
      return(stats::quantile(CSI[is.finite(CSI)]*clearsky_GHI[i], probs=percentiles, names=F, na.rm=T, type=1))
    } else return(rep(0, times=length(percentiles)))}, simplify="array"))
  
  colnames(fc) <- percentiles
  return(fc)
}

#' Complete-history persistence ensemble method
#' 
#' @param tel A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
forecast_Ch_PeEn <- function(tel, percentiles, sun_up) {
  n <- ncol(tel)
  daily_fc <- sapply(1:n, FUN=function(i)stats::quantile(tel[sun_up[,i],i], probs=percentiles, names=F, na.rm=T, type=1))
  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {if (isTRUE(as.vector(t(sun_up))[i]))  daily_fc[,(i-1)%%n+1]
    else rep(0, times=length(percentiles))}))
  colnames(fc) <- percentiles
  return(fc)
}

#' Gaussian error dressing: hourly, using hourly errors of the deterministic NWP forecast to fit truncated normal distributions
#' 
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble forecasts; First member is treated as the control member.
#' @param GHI A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
forecast_Gaussian_hourly <- function(nwp, GHI, percentiles, sun_up) {
  if ((dim(nwp)[3]/dim(nwp)[2])%%1!=0 | (dim(nwp)[3]/dim(nwp)[2])==1) stop("Unknown handling for given number of issue and lead times")
  if (any(dim(GHI)!=dim(sun_up))) stop("Given incompatible number of forecasts and sun_up times")
  
  n <- ncol(GHI)
  nwp_ctrl <- matrix(aperm(nwp[,,1:(dim(nwp)[3]/dim(nwp)[2]),1], c(3,2,1)), nrow=nrow(nwp), byrow=T) # Reshape to [day x hour]
  # Fit a standard deviation for each hour of the day, based on all residuals from that hour
  daily_sd <- sapply(1:n, FUN=function(i) sd(nwp_ctrl[,i] - GHI[,i], na.rm=T))

  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {if (isTRUE(as.vector(t(sun_up))[i])) qtruncnorm(p=percentiles, a=0, mean=as.vector(t(nwp_ctrl))[i], sd=daily_sd[(i-1)%%n+1])
    else rep(0, times=length(percentiles))}))
  colnames(fc) <- percentiles
  return(fc)
}

#' Gaussian error dressing: intra-hour, using smart persistence errors from the past hour
#'
#' @param GHI A vector of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A vector of logicals, indicating whether the sun is up
#' @param clearsky_GHI a vector of clear-sky irradiance estimates
#' @param n Time steps per hour, e.g., 60 for a 1-minute resolution forecast
forecast_Gaussian_minute <- function(GHI, percentiles, sun_up, clearsky_GHI, n) {
  
  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {
    if (isTRUE(sun_up[i])) {
      # Forecast is undefined in first time points when no error data is available
      if (i <= 2) {return(rep(NA, times=length(percentiles)))}
      else {
        
        # Get clear-sky indices from previous n time-steps, ending 2 steps before forecast
        CSI <- GHI[max(i-n-1, 0):max(i-2, 0)]/clearsky_GHI[max(i-n-1, 0):max(i-2, 0)] 
        # Use clear-sky indices to calculate smart persistence for previous (up to) n time-steps, ending 1 step before forecast
        smart_persistence <- clearsky_GHI[(i-length(CSI)):(i-1)]*CSI 
        errors <- smart_persistence - GHI[(i-length(CSI)):(i-1)]
        # Force deterministic forecast if errors are all NA/Inf in the training data (e.g., as the sun comes up)
        if (!any(is.finite(CSI))) errors <- rep(0, times=n) 
        
        CSI_forecast <- GHI[i-1]/clearsky_GHI[i-1]
        # Force deterministic smart persistence forecast when no training data is available
        smart_persistence_now <- clearsky_GHI[i]*ifelse(is.finite(CSI_forecast), CSI_forecast, 1)
        return(qtruncnorm(p=percentiles, a=0, mean=smart_persistence_now, sd=sd(errors[is.finite(errors)], na.rm=T)))  
      }
    } else return(rep(0, times=length(percentiles)))}, simplify="array"))
  
  colnames(fc) <- percentiles
  return(fc)
}
