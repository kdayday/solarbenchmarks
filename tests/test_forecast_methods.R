library(testthat)
library(truncnorm)

# Currently set up to run from the Benchmarks_comparison folder
source("forecast_methods.R")

test_that("forecast_climatology calculation is correct", {
  percentiles <- seq(0.1, 0.9, by=0.1)
  tel <- c(15, 3, 6, 2, 7, 8, 1, 0, 0, 11, 9, 12, 19)
  clim <- c(0, 1, 2, 3, 6, 7, 8, 9, 11)
  fc <- rbind(rep(0, times=9), clim, clim, clim, clim, clim, clim, clim, clim, clim, clim, clim, rep(0, times=9))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_climatology(tel, percentiles=percentiles, sun_up = c(F, rep(T, 11), F)), fc)
})

test_that("forecast_nwp throws errors", {
  nwp <- array(1:24, dim=c(2, 2, 2, 3))
  expect_error(forecast_NWP(nwp, percentiles=NA, sun_up=NA), "Unknown handling*")
  nwp <- array(1:24, dim=c(2, 2, 3, 2))
  expect_error(forecast_NWP(nwp, percentiles=NA, sun_up=NA), "Unknown handling*")
  sun_up <- matrix(c(T,T,T,F), ncol=2)
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  expect_error(forecast_NWP(nwp, percentiles=NA, sun_up=sun_up), "Given incompatible*")
})

test_that("forecast_nwp calculation is correct", {
  percentiles <- c(0.25, 0.5, 0.75)
  sun_up <- matrix(c(rep(T,7), F), ncol=2)
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  fc <- rbind(c(1,17,33),
              c(5,21,37),
              c(3,19,35),
              c(7,23,39),
              c(2,18,34),
              c(6,22,38),
              c(4,20,36),
              c(0,0,0))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_NWP(nwp, percentiles, sun_up), fc)
})

test_that("forecast_PeEn_hourly throws errors", {
  lead_up_GHI <- matrix(c(0:2, 5:7), ncol=2)
  expect_error(forecast_PeEn_hourly(GHI=NA, percentiles=NA, sun_up=NA, num_peen=5, lead_up_GHI=lead_up_GHI), "days of data*")
})

# Includes test of longer lead_up_GHI than needed, sun down time, and forecasts both with and without lead_up_GHI data
test_that("forecast_PeEn_hourly calculation is correct", {
  tel <- matrix(c(3:5, 8:10), ncol=2)
  percentiles <- c(0.25, 0.5, 0.75)
  lead_up_GHI <- matrix(c(0:2, 5:7), ncol=2)
  sun_up <- matrix(c(rep(T, times=5), F), ncol=2)
  fc <- rbind(c(1,1,2),
              c(6,6,7),
              c(2,2,3),
              c(7,7,8),
              c(3,3,4),
              c(0,0,0))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_PeEn_hourly(tel, percentiles, sun_up, num_peen=2, lead_up_GHI=lead_up_GHI), fc)
})

# Include test of edge cases with no training data yet or no non-NA training data
test_that("forecast_PeEn_minute calculation is correct", {
  sun_up <- c(T, F, F, T, T, T, F)
  tel <- c(0, 0, 1, 9, 16, 10, 3)
  clearsky_GHI <- c(10, 0, 0, 12, 20, 40, 30)
  # CSI <- (NULL, NA, NA, 0.75, 0.8, 0.25, 0.1)
  percentiles <- c(0.25, 0.5, 0.75)
  fc <- rbind(c(10, 10, 10), # Deterministic forecast: start of vector
              c(0, 0, 0), # sundown
              c(0, 0, 0), # sundown
              c(12, 12, 12), # Deterministic forecast: Start of day
              c(15, 15, 15), # Deterministic forecast: only one error so far
              c(30, 30, 32), # Distributional forecast: two errors available
              c(0, 0, 0)) # sun down with non-zero cSI
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_PeEn_minute(tel, percentiles, sun_up, num_peen=2, clearsky_GHI = clearsky_GHI), fc)
})

test_that("forecast_Ch_PeEn calculation is correct", {
  tel <- matrix(c(1, 20, 300, 4, 50, 60, 7, 8, 900, 16), ncol=2)
  clearsky_GHI <- matrix(c(10, 100, 1000, 10, 100, 100, 10, 10, 1000, 10), ncol=2)
  percentiles <- c(0.25, 0.5, 0.75)
  fc1 <- c(0.2, 0.3, 0.4)
  fc2 <- c(0.6, 0.7, 0.8) # Last point is excluded because sun down
  fc <- rbind(fc1*10,
              fc2*100,
              fc1*100,
              fc2*10,
              fc1*1000,
              fc2*10,
              fc1*10,
              fc2*1000,
              fc1*100,
              rep(0, times=3))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_Ch_PeEn(tel, percentiles, cbind(rep(T, 5), c( T, T, T, T, F)), clearsky_GHI), fc)
})

test_that("forecast_Gaussian_hourly throws errors", {
  nwp <- array(1:24, dim=c(2, 3, 2, 2))
  expect_error(forecast_Gaussian_hourly(nwp=nwp, GHI=NA, percentiles=NA, sun_up=NA), "Unknown handling*")
  nwp <- array(1:24, dim=c(2, 2, 2, 3))
  expect_error(forecast_Gaussian_hourly(nwp=nwp, GHI=NA, percentiles=NA, sun_up=NA), "Unknown handling*")
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  sun_up <- matrix(c(T,T,T,F), ncol=2)
  tel <- matrix(1:6, ncol=3)
  expect_error(forecast_Gaussian_hourly(nwp=nwp, GHI=tel, percentiles=NA, sun_up=sun_up), "Given incompatible*")
})

test_that("forecast_Gaussian_hourly calculation is correct", {
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  tel <- matrix(rep(2, 8), ncol=4, byrow = T)
  percentiles <- c(0.25, 0.5, 0.75)
  sun_up <- matrix(c(rep(T, times=7), F), ncol=4)
  # residuals (-1 + 0), (3 + 4), (1 + 2), (5 + 6) = -1, 7, 3, 11
  fc <- rbind(rep(0, times=3),
              rep(12, times=3),
              rep(6, times=3),
              rep(18, times=3), 
              rep(1, times=3), 
              rep(13, times=3), 
              rep(7, times=3), 
              rep(0, times=3))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  with_mock(sd=function(...) return(sum(...)), 
            qtruncnorm=function(p, a, mean, sd) return(rep(mean+sd, times=length(p))),
            expect_equal(forecast_Gaussian_hourly(nwp, tel, percentiles, sun_up), fc))
})

# Include test of edge cases with no training data yet or no non-NA training data
test_that("forecast_Gaussian_minute calculation is correct", {
  sun_up <- c(T, T, T, T, F, F, F, T, T, T)
  GHI <- c(7, 8, 10, 16, 0, 0, 10, 9, 8, 16)
  clearsky_GHI <- c(10, 10, 10, 20, 0, 0, 0, 12, 16, 20)
  # csi (0.7, 0.8, 1, 0.8, NA, NA, Inf, 0.75, 0.5, 0.8)
  # smart persistence (NA, 7, 8, 20, NA, NA, NA, 12, 12, 10)
  # errors (NA, -1, -2, 4, NA, NA, NA, Inf, -4, -10) 
  # sum of errors (NA, NA, -1, -3, NA, NA, NA, NA, 4)
  percentiles <- c(0.25, 0.5, 0.75)
  fc <- rbind(rep(NA, times=3),
              rep(NA, times=3),
              rep(7, times=3),
              rep(17, times=3),
              rep(0, times=3),
              rep(0, times=3),
              rep(0, times=3),
              rep(12, times=3),
              rep(12, times=3),
              rep(14, times=3))
  
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  with_mock(sd=function(...) return(sum(...)),
            qtruncnorm=function(p, a, mean, sd) return(rep(mean+sd, times=length(p))),
            expect_equal(forecast_Gaussian_minute(GHI, percentiles, sun_up, clearsky_GHI, n=2), fc))
})
