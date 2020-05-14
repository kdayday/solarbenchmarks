# This code reproduced courtesy of the authors, based on the method in: J.
# Munkhammar, J. Widén, D. W. van der Meer, Probabilistic forecasting of
# high-resolution clear-sky index time-series using a Markov-chain mixture
# distribution model, Solar Energy vol. 184, pp. 688-695, 2019. 
# Original code can be accessed in the mcmFunctions.R script at:
# https://github.com/SheperoMah/MCM-distribution-forecasting


#' @export
mcmFit <- function(timeSeries, numBins, numStepAhead=1){
  a <- min(timeSeries)
  b <- max(timeSeries)
  binWidth <- (b-a) / numBins
  
  state <- floor((timeSeries - a) / binWidth) + 1
  state[state > numBins] <- numBins # remove any state > numberOfStates
  
  p <- matrix(0, numBins, numBins)
  
  for (i in 2:length(timeSeries)){
    p[state[i-1], state[i]] <- p[state[i-1], state[i]] + 1
  } 
  
  #normalize p
  pNormalized <- p / rowSums(p)
  pNormalized[is.na(pNormalized)] <- 0.0
  
  pStepped <- matrixcalc::matrix.power(pNormalized, numStepAhead)
  return(pStepped)
}

#' @export
mcmForecast <- function(p, a, b, obs) {
  stopifnot(a < b) # a should be less than b
  stopifnot(diff(dim(p)) == 0) # p should be square matrix
  stopifnot(obs >= a) # observation should be larger than the samllest value in the range.
  
  numBins <- dim(p)[1]
  binWidth <- (b-a) / numBins
  
  binStartArray <- a + binWidth * (0:(numBins-1))
  obsBin <- max(which(obs >= binStartArray))
  
  returnList <- list(binStartingValues = binStartArray, 
                     transitionProbs = p[obsBin,])
  return(returnList)
}

#' @export
mcmRnd <- function(binStartVals, transProbs, numSamples){
  binWidth <- diff(binStartVals)[1]
  numBins <- length(binStartVals)
  
  empiricalCdf <- cumsum(transProbs)
  
  randWithinBin <- runif(numSamples, min=0, max=binWidth)
  randWithinECDF <- runif(numSamples, min=0, max=1.0)
  
  rndValue <- rep(0, numSamples)
  for (i in 1:numSamples){
    binIndex <- min(which(randWithinECDF[i] <= empiricalCdf))
    rndValue[i] <- binStartVals[binIndex] + randWithinBin[i]
  }
  return(rndValue)
}
