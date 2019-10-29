

#' Get the quantile-weighted CRPS (continuous ranked probability score) for the forecast
#' Estimated by trapezoidal integration of its quantile decomposition
#' Quantile weighting functions as suggested in Gneiting & Ranjan 2012
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @return The weighted CRPS
qwCRPS <-function(fc, tel, sun_up, weighting="none"){
  qs <- QS(fc, tel, sun_up, as.numeric(colnames(fc)))

  wqs <- weight_QS(qs, as.numeric(colnames(fc)), weighting)
  return(trapz(as.numeric(colnames(fc)), wqs))
}
 
#' Get quantile score at each quantile
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param percentiles A vector (0,1) corresponding to the columns of fc
#' @return the Quantile scores at each quantile in the fc matrix
QS <- function(fc, tel, sun_up, percentiles) {
  if (nrow(fc) != length(tel) | length(sun_up)!=length(tel)) stop("Forecast, telemetry, and sun_up must have the same time range and resolution")
  
  indicator <- sapply(seq_len(ncol(fc)), FUN=function(i) {as.integer(tel <= fc[,i])}, simplify="array")
  qs <- sapply(seq_len(ncol(fc)), FUN=function(q) mean(2*(indicator[sun_up,q]-percentiles[q])*(fc[sun_up,q]-tel[sun_up]), na.rm = TRUE))
  return(qs)
}

#' Calculate a vector of weighted quantile scores, emphasizing one or both tails or center
#'
#' @param qs A vector of quantile scores
#' @param percentiles A vector of the percentiles in [0,1]
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @return A vector of the weighted scores
weight_QS <- function(qs, percentiles, weighting="none") {
  if (length(qs) != length(percentiles)) stop("Quantiles and quantile score must be the same length")
  
  weights <- switch(tolower(weighting),
                    "none" = 1,
                    "tails" = (2*percentiles-1)^2,
                    "right" = percentiles^2,
                    "left" = (1-percentiles)^2,
                    "center" = percentiles*(1-percentiles),
                    stop("Weighting method not recognized"))
  return(weights*qs)
}

#' Calculate average sharpness for central intervals from 10% to 90%. Negatively oriented
#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param intervals A vector of central intervals in (0,1) whose widths should be calculated. Should correspond to available percentile columns in fc.
#' @return a list of the intervals and their average widths
interval_width <- function(fc, sun_up, intervals=seq(0.1, 0.9, by=0.1)) {
  percentiles <- as.numeric(colnames(fc))
  
  widths <- sapply(intervals, FUN=function(i) mean(fc[sun_up,as.character(0.5+i/2)] - fc[sun_up,as.character(0.5-i/2)], na.rm=T))
  return(list(widths=widths, intervals=intervals))
}

#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param percentiles A vector of the percentiles corresponding to the desired forecast quantiles
reliability <- function(fc, tel, sun_up, percentiles) {
  percentiles <- as.numeric(colnames(fc))
  counts <- sapply(seq_along(percentiles), FUN= function(i) {sum(tel[sun_up] <= fc[sun_up,i], na.rm=T)})
  obs_proportion <- counts/sum(sun_up, na.rm=T)
  return(obs_proportion)
}

#' Plot a fan plot comparing the forecast quantiles to the telemetry
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param window A (min, max) vector of indices, corresponding to the rows of fc to plot
plot_fanplot <- function(fc, tel, window) {
  df <- rowid_to_column(data.frame(fc[window[1]:window[2],], check.names=F), var="x") 
  df <- gather(df, key="q", value="y", -x, convert=T)
  g <- ggplot(df,  aes(x=x,y=y,quantile=q)) + geom_fan() +
    geom_line(mapping=aes(x=x, y=y), data=data.frame(x=seq_len(diff(window)+1), y=tel[window[1]:window[2]]), col="chocolate3", inherit.aes = F) +
    theme_bw() + xlab("UTC Time (hours)") + ylab(expression(paste("Irradiance (W/m"^"2", ")")))  +
    scale_x_continuous(breaks=seq(0, diff(window)+1, by=6))
}
