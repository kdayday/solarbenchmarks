

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
    theme_bw() + xlab("Time (hours)") + ylab(expression(paste("Irradiance (W/m"^"2", ")")))  +
    scale_x_continuous(breaks=seq(0, diff(window)+1, by=6))
}

# Plot reliability of 1..99 quantiles, comparing among forecast methods for a single site
#' @param reliability Data frame of the quantile reliability results for each method
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as well
plot_reliability <- function(reliability_df, site, res, shapes, output_directory, R_graph_export) {
  df <- melt(reliability_df[site,,], varnames=c("Method", "levels"))
  
  g <- ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
    xlab("Nominal proportion") + ylab("Observed proportion") + 
    scale_shape_manual(values = shapes[1:nrow(reliability_df)]) + 
    geom_line(data = data.frame(x=c(0,1), y=c(0,1)), mapping=aes(x=x,y=y), col="black") + 
    theme(legend.justification=c(1,0), legend.position=c(0.98,0.02), text = ggplot2::element_text(size=14))
  ggsave(file.path(output_directory, paste(site, res, "reliability.pdf", sep="_")), height=4, width=8)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, "reliability_plot.R", sep="_")))
  }
}

# Plot average interval width of 10%, ..., 90% central intervals, comparing among forecast methods for a single site
#' @param interval_width_df Data frame of the interval width results for each method
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as well
plot_interval_width <- function(interval_width_df, site, res, shapes, output_directory, R_graph_export) {
  
  df <- melt(interval_width_df[site,,], varnames=c("Method", "levels"))
  
  g <- ggplot(df, aes(levels,value)) + geom_point(aes(colour = Method, shape=Method)) + geom_line(aes(colour = Method)) + 
    xlab("Central Interval (%)") + ylab(expression(paste("Average Width (W/m"^"2", ")"))) + 
    scale_shape_manual(values = shapes[1:nrow(interval_width_df)]) + 
    theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), text = ggplot2::element_text(size=14)) + 
    scale_x_continuous(breaks=intervals, labels=intervals*100)
  ggsave(file.path(output_directory, paste(site, res, "interval_width.pdf", sep="_")), height=4, width=8)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, "interval_width_plot.R", sep="_")))  
  }
}
