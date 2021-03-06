% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation_functions.R
\name{weight_QS}
\alias{weight_QS}
\title{Calculate weighted quantile scores}
\usage{
weight_QS(qs, percentiles, weighting = "none")
}
\arguments{
\item{qs}{A vector of quantile scores}

\item{percentiles}{A vector of the percentiles in [0,1]}

\item{weighting}{One of "none" (default), "tails", "right", "left", "center"}
}
\value{
A vector of the weighted scores
}
\description{
Calculate a vector of weighted quantile scores, emphasizing one or both tails
or center. Quantile weighting functions available as suggested in:
}
\details{
Gneiting, Tilmann, and Roopesh Ranjan. "Comparing density forecasts using
threshold-and quantile-weighted scoring rules." Journal of Business &
Economic Statistics 29.3 (2011): 411-422.
}
