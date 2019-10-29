
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

#' Persistence ensemble method
#' 
#' @param tel A [day x hour] matrix of the hourly average telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x hour] matrix of logicals, indicating whether the sun is up
#' @param num_peen  Number of persistence ensemble members to take
#' @param oos_tel A [day x hour] matrix of out-of-sample telemetry for the days leading up to the start of the sample period. Number of days must be >= num_peen
forecast_PeEn <- function(tel, percentiles, sun_up, num_peen, oos_tel) {
  if (dim(oos_tel)[1] < num_peen) stop(paste("Must have at least", num_peen, "days of data leading up to the sample people available. Given", dim(oos_tel)[1]))
  
  all_tel <- rbind(oos_tel, tel)
  
  fc <- sapply(seq(nrow(oos_tel)+1, nrow(all_tel)), FUN=function(day) {
    sapply(seq_len(ncol(tel)), FUN=function(hr) {if (isTRUE(sun_up[day-nrow(oos_tel), hr])) {
      stats::quantile(all_tel[(day-num_peen):(day-1),hr], probs=percentiles, names=F, na.rm=T, type=1)} else 
        rep(0, times=length(percentiles))})
  }, simplify="array")
  fc <- array(aperm(fc, c(2,3,1)), dim=c(length(tel), length(percentiles)))
  
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

#' NWP Gaussian error dressing, using hourly errors of the deterministic NWP forecast to fit truncated normal distributions.
#' 
#' @param nwp A [day x issue time x lead time x member] matrix of NWP ensemble forecasts; First member is treated as the control member.
#' @param tel A [day x time] matrix of the telemetry
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
#' @param sun_up A [day x time] matrix of logicals, indicating whether the sun is up
forecast_Gaussian <- function(nwp, tel, percentiles, sun_up) {
  if ((dim(nwp)[3]/dim(nwp)[2])%%1!=0 | (dim(nwp)[3]/dim(nwp)[2])==1) stop("Unknown handling for given number of issue and lead times")
  if (any(dim(tel)!=dim(sun_up))) stop("Given incompatible number of forecasts and sun_up times")
  
  n <- ncol(tel)
  nwp_ctrl <- matrix(aperm(nwp[,,1:(dim(nwp)[3]/dim(nwp)[2]),1], c(3,2,1)), nrow=nrow(nwp), byrow=T) # Reshape to [day x hour]
  # Fit a standard deviation for each hour of the day, based on all residuals from that hour
  daily_sd <- sapply(1:n, FUN=function(i) sd(nwp_ctrl[,i] - tel[,i], na.rm=T))
  
  fc <- t(sapply(seq_along(sun_up), FUN=function(i) {if (isTRUE(as.vector(t(sun_up))[i])) qtruncnorm(p=percentiles, a=0, mean=as.vector(t(nwp_ctrl))[i], sd=daily_sd[(i-1)%%n+1])
    else rep(0, times=length(percentiles))}))
  colnames(fc) <- percentiles
  return(fc)
}