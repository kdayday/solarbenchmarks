
#' Calculate quantile-weighted CRPS
#'
#' Get the quantile-weighted CRPS (continuous ranked probability score) for the
#' forecast, estimated by trapezoidal integration of its quantile decomposition.
#' Can be weighted or unweighted (default).
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @return The weighted CRPS
#' @export
qwCRPS <-function(fc, tel, sun_up, weighting="none"){
  qs <- QS(fc, tel, sun_up, as.numeric(colnames(fc)))

  wqs <- weight_QS(qs, as.numeric(colnames(fc)), weighting)
  return(pracma::trapz(as.numeric(colnames(fc)), wqs))
}
 
#' Get quantile score at each quantile
#' 
#' Calculate quantile score at each quantile in the fc matrix
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile forecasts
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param percentiles A vector (0,1) corresponding to the columns of fc
#' @return the quantile scores at each quantile in the fc matrix
#' @export
QS <- function(fc, tel, sun_up, percentiles) {
  if (nrow(fc) != length(tel) | length(sun_up)!=length(tel)) stop("Forecast, telemetry, and sun_up must have the same time range and resolution")
  
  indicator <- sapply(seq_len(ncol(fc)), FUN=function(i) {as.integer(tel <= fc[,i])}, simplify="array")
  qs <- sapply(seq_len(ncol(fc)), FUN=function(q) mean(2*(indicator[sun_up,q]-percentiles[q])*(fc[sun_up,q]-tel[sun_up]), na.rm = TRUE))
  return(qs)
}

#' Calculate weighted quantile scores
#'
#' Calculate a vector of weighted quantile scores, emphasizing one or both tails
#' or center. Quantile weighting functions available as suggested in:
#'
#' Gneiting, Tilmann, and Roopesh Ranjan. "Comparing density forecasts using
#' threshold-and quantile-weighted scoring rules." Journal of Business &
#' Economic Statistics 29.3 (2011): 411-422.
#'
#' @param qs A vector of quantile scores
#' @param percentiles A vector of the percentiles in [0,1]
#' @param weighting One of "none" (default), "tails", "right", "left", "center"
#' @return A vector of the weighted scores
#' @export
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

#' Calculate average interval width
#'
#' Calculate average sharpness for central intervals from 10% to 90%. Negatively
#' oriented (smaller is better).
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param intervals A vector of central intervals in (0,1) whose widths should
#'   be calculated. Should correspond to available percentile columns in fc.
#' @return a list of the intervals and their average widths
#' @export
interval_width <- function(fc, sun_up, intervals=seq(0.1, 0.9, by=0.1)) {
  percentiles <- as.numeric(colnames(fc))
  
  widths <- sapply(intervals, FUN=function(i) mean(fc[sun_up,as.character(0.5+i/2)] - fc[sun_up,as.character(0.5-i/2)], na.rm=T))
  return(list(widths=widths, intervals=intervals))
}

#' Calculate reliability at each quantile
#'
#' Calculate proportion of observations that fell within a given quantile,
#' relative to each nominal proportion (percentile).
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @return A vector of observed proportions corresponding to each column in fc
#' @export
reliability <- function(fc, tel, sun_up) {
  percentiles <- as.numeric(colnames(fc))
  counts <- sapply(seq_along(percentiles), FUN= function(i) {sum(tel[sun_up] <= fc[sun_up,i], na.rm=T)})
  obs_proportion <- counts/sum(sun_up, na.rm=T)
  return(obs_proportion)
}

#' Plot probability integral transform histogram
#'
#' Plot PIT histogram for a single single and forecast method
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param nbins The number of histogram bins to use. Must be a factor of
#'   length(percentiles)+1 to make use of available fc format.
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param method Forecast method name, for naming output file
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as
#'   well
#' @export
plot_PIT_histogram  <- function(fc, tel, sun_up, nbins, site, res, method, output_directory, R_graph_export) {
  if ((ncol(fc)+1)%%nbins != 0 ) stop("nbins must be a factor of length(percentiles)+1 to make use of available fc format.")

  percentiles <- as.numeric(colnames(fc))
  breaks_indices <- 1:(nbins)*(length(percentiles)+1)/nbins
  
  valid <- is.finite(tel) & sun_up # Only assess metrics over times when sun is up and data is not missing
  bin_hits <- sapply(seq_along(breaks_indices), FUN = function(i) {
    if (i==1) {low_limit <- rep(-Inf, times=nrow(fc))} else {low_limit <- fc[,breaks_indices[i-1]]}
    if (i==length(breaks_indices)) {upper_limit <- rep(Inf, times=nrow(fc))} else {upper_limit <- fc[,breaks_indices[i]]}
    sum(tel[valid] > low_limit[valid] & tel[valid] <= upper_limit[valid], na.rm=T)
  })
  
  bin_width <- diff(percentiles[breaks_indices[1:2]])
  bin_means <- c(percentiles[breaks_indices[-nbins]],1) - bin_width/2
  
  g <- ggplot2::ggplot(data.frame(x=bin_means, y=bin_hits/sum(bin_hits)), ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_col(width=bin_width) +
    ggplot2::geom_line(data.frame(x=c(0, 1), y=c(1/nbins, 1/nbins)), mapping=ggplot2::aes(x=x, y=y), linetype="dashed") +
    ggplot2::xlab("Probability Integral Transform") +
    ggplot2::ylab("Relative Frequency")
  ggplot2::ggsave(file.path(output_directory, paste(site, res, method, "PIT_histogram.pdf", sep="_")), height=4, width=6)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, method, "PIT_histogram.R", sep="_")))
  }
  return(g)
}

#' Plot quantile score vs. quantile
#' 
#' Plot unweighted, left-, and right-weighted quantile scores of selected
#' quantiles
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param sun_up Logical vector of whether sun is up across the forecast times
#' @param percentiles A vector (0,1) corresponding to the columns of fc
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param method Forecast method name, for naming output file
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as
#'   well
#' @export
plot_quantile_score  <- function(fc, tel, sun_up, percentiles, site, res, method, output_directory, R_graph_export) {
  unweighted <- QS(fc, tel, sun_up, percentiles)
  l_weighted <- weight_QS(unweighted, percentiles, weighting="left")
  r_weighted <- weight_QS(unweighted, percentiles, weighting="right")

  df <- data.frame(a=unweighted, b=l_weighted, c=r_weighted, q=percentiles)
  colnames(df) <- c("Unweighted", "Left-weighted", "Right-weighted", "quantiles")
  df <- reshape2::melt(df, id.vars=c("quantiles"))
  
  g <- ggplot2::ggplot(df, ggplot2::aes(quantiles,value)) + 
    ggplot2::geom_point(ggplot2::aes(colour = variable, shape=variable)) + 
    ggplot2::geom_line(ggplot2::aes(colour = variable)) +
    ggplot2::xlab("Quantile") +
    ggplot2::ylab(expression(paste("Quantile Score (W/m"^"2", ")"))) + 
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.justification=c(0,1), 
                   legend.position=c(0.02,0.98), text = ggplot2::element_text(size=14))
  ggplot2::ggsave(file.path(output_directory, paste(site, res, method, "quantile_scores.pdf", sep="_")), 
         height=4, width=6)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, method, "quantile_scores.R", sep="_")))
  }
  return(g)
}

#' Plot fan plot of time-series forecasts
#'
#' Plot a fan plot comparing the forecast quantiles to the telemetry
#'
#' @param fc A [valid time x quantile] matrix of probabilistic quantile
#'   forecasts, with column names giving the [0,1] percentiles of the forecast
#' @param tel A vector of the telemetry values
#' @param window A (min, max) vector of time indices, corresponding to the rows
#'   of fc to plot
#' @param ts_per_hour Time-steps per hour
#' @param output_directory Directory to save graphs
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param method Forecast method name, for naming output file
#' @param R_graph_export Boolean, whether to save a .R object of the plot as
#'   well
#' @param ymax Y-axis maximum value (defaults to 1200 W/m^2)
#' @export
plot_fanplot <- function(fc, tel, window, ts_per_hour, output_directory, site, res, method, R_graph_export, ymax=1200) {
  df <- tibble::rowid_to_column(data.frame(fc[window[1]:window[2],], check.names=F), var="x") 
  df <- tidyr::gather(df, key="q", value="y", -x, convert=T)
  g <- ggplot2::ggplot(df,  ggplot2::aes(x=x,y=y,quantile=q)) + ggfan::geom_fan() +
    ggplot2::geom_line(mapping=ggplot2::aes(x=x, y=y), data=data.frame(x=seq_len(diff(window)+1), y=tel[window[1]:window[2]]), col="chocolate3", inherit.aes = F) +
    ggplot2::scale_fill_continuous(limits=c(0.02, 0.98), breaks=c(0.02, 0.250, 0.500, 0.750, 0.98), 
                          guide=ggplot2::guide_colourbar(nbin=98, draw.ulim=F, draw.llim=F)) +
    ggplot2::theme_bw() + ggplot2::xlab("Time (hours)") + ggplot2::ylab(expression(paste("Irradiance (W/m"^"2", ")"))) + 
    ggplot2::coord_cartesian(ylim = c(0,1200)) + # Don't truncate data if extends past window
    ggplot2::scale_y_continuous(breaks=seq(from=0, to=ymax, by=200)) +
    ggplot2::theme(text=ggplot2::element_text(size=14), legend.title = ggplot2::element_text(vjust=1)) + 
    # Add x scale with every 6 hours
    ggplot2::scale_x_continuous(breaks=seq(from=0, to=diff(window)+1, by=6*ts_per_hour), 
                       labels=seq(from=0, to=(diff(window)+1)/ts_per_hour, by=6))
  ggplot2::ggsave(file.path(output_directory, paste(site, res, method, "sample.pdf", sep="_")), height=3, width=6)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, method, "sample.R", sep="_")))  
  }
  return(g)
}

#' Plot average quantile reliability
#' 
#' Plot reliability of 1..99 quantiles, comparing among forecast methods for a
#' single site
#' @param reliability_df Data frame of the quantile reliability results for each
#'   method
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param shapes List of ggplot2 shapes, one for each method
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as
#'   well
#' @export
plot_reliability <- function(reliability_df, site, res, shapes, output_directory, R_graph_export) {
  df <- reshape2::melt(reliability_df[site,,], varnames=c("Method", "levels"))
  
  g <- ggplot2::ggplot(df, ggplot2::aes(levels,value)) + 
    ggplot2::geom_point(ggplot2::aes(colour = Method, shape=Method)) + 
    ggplot2::geom_line(ggplot2::aes(colour = Method)) + 
    ggplot2::xlab("Nominal proportion") + ggplot2::ylab("Observed proportion") + 
    ggplot2::scale_shape_manual(values = shapes[1:nrow(reliability_df)]) + 
    ggplot2::geom_line(data = data.frame(x=c(0,1), y=c(0,1)), mapping=ggplot2::aes(x=x,y=y), col="black") + 
    ggplot2::theme(legend.justification=c(1,0), legend.position=c(0.98,0.02), text = ggplot2::element_text(size=14))
  ggplot2::ggsave(file.path(output_directory, paste(site, res, "reliability.pdf", sep="_")), height=4, width=8)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, "reliability_plot.R", sep="_")))
  }
  return(g)
}

#' Plot average interval width
#' 
#' Plot average interval width of 10%, ..., 90% central intervals, comparing
#' among forecast methods for a single site
#' @param interval_width_df Data frame of the interval width results for each
#'   method
#' @param site Site name, for naming output file
#' @param res Temporal resolution, for naming output file
#' @param shapes List of ggplot2 shapes, one for each method
#' @param output_directory Directory to save graphs
#' @param R_graph_export Boolean, whether to save a .R object of the plot as
#'   well
#' @export
plot_interval_width <- function(interval_width_df, site, res, shapes, output_directory, R_graph_export) {
  
  df <- reshape2::melt(interval_width_df[site,,], varnames=c("Method", "levels"))
  
  g <- ggplot2::ggplot(df, ggplot2::aes(levels,value)) + 
    ggplot2::geom_point(ggplot2::aes(colour = Method, shape=Method)) + ggplot2::geom_line(ggplot2::aes(colour = Method)) + 
    ggplot2::xlab("Central Interval (%)") + ggplot2::ylab(expression(paste("Average Width (W/m"^"2", ")"))) + 
    ggplot2::scale_shape_manual(values = shapes[1:nrow(interval_width_df)]) + 
    ggplot2::theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), text = ggplot2::element_text(size=14)) + 
    ggplot2::scale_x_continuous(breaks=intervals, labels=intervals*100)
  ggplot2::ggsave(file.path(output_directory, paste(site, res, "interval_width.pdf", sep="_")), height=4, width=8)
  if (R_graph_export) {
    save(g, file=file.path(output_directory, paste(site, res, "interval_width_plot.R", sep="_")))  
  }
  return(g)
}
