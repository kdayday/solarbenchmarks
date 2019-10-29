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

test_that("forecast_PeEn throws errors", {
  oos_tel <- matrix(c(0:2, 5:7), ncol=2)
  expect_error(forecast_PeEn(tel=NA, percentiles=NA, sun_up=NA, num_peen=5, oos_tel=oos_tel), "days of data*")
})

# Includes test of longer oos_tel than needed, sun down time, and forecasts both with and without oos_tel data
test_that("forecast_PeEn calculation is correct", {
  tel <- matrix(c(3:5, 8:10), ncol=2)
  percentiles <- c(0.25, 0.5, 0.75)
  oos_tel <- matrix(c(0:2, 5:7), ncol=2)
  sun_up <- matrix(c(rep(T, times=5), F), ncol=2)
  fc <- rbind(c(1,1,2),
              c(6,6,7),
              c(2,2,3),
              c(7,7,8),
              c(3,3,4),
              c(0,0,0))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_PeEn(tel, percentiles, sun_up, num_peen=2, oos_tel=oos_tel), fc)
})

test_that("forecast_Ch_PeEn calculation is correct", {
  tel <- matrix(c(1:5, c(1:4, 15)), ncol=2)
  percentiles <- c(0.25, 0.5, 0.75)
  fc1 <- 2:4
  fc2 <- 1:3
  fc <- matrix(c(rep(c(fc1, fc2), times=4), fc1, rep(0, times=3)), ncol=3, byrow=T)
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  expect_equal(forecast_Ch_PeEn(tel, percentiles, cbind(rep(T, 5), c( T, T, T, T, F))), fc)
})

test_that("forecast_Gaussian throws errors", {
  nwp <- array(1:24, dim=c(2, 3, 2, 2))
  expect_error(forecast_Gaussian(nwp=nwp, tel=NA, percentiles=NA, sun_up=NA), "Unknown handling*")
  nwp <- array(1:24, dim=c(2, 2, 2, 3))
  expect_error(forecast_Gaussian(nwp=nwp, tel=NA, percentiles=NA, sun_up=NA), "Unknown handling*")
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  sun_up <- matrix(c(T,T,T,F), ncol=2)
  tel <- matrix(1:6, ncol=3)
  expect_error(forecast_Gaussian(nwp=nwp, tel=tel, percentiles=NA, sun_up=sun_up), "Given incompatible*")
})

test_that("forecast_Gaussian calculation is correct", {
  nwp <- array(1:48, dim=c(2, 2, 4, 3))
  tel <- matrix(rep(2, 8), ncol=4, byrow = T)
  percentiles <- c(0.25, 0.5, 0.75)
  sun_up <- matrix(c(rep(T, times=7), F), ncol=4)
  fc <- rbind(rep(0, times=3),
              rep(12, times=3),
              rep(7, times=3),
              rep(19, times=3), 
              rep(1, times=3), 
              rep(13, times=3), 
              rep(7, times=3), 
              rep(0, times=3))
  colnames(fc) <- percentiles
  rownames(fc) <- NULL
  with_mock(sd=function(...) return(sum(...)), 
            qtruncnorm=function(p, a, mean, sd) return(rep(mean+sd, times=length(p))),
            out <- forecast_Gaussian(nwp, tel, percentiles, sun_up))
            # expect_equal(forecast_Gaussian(nwp, tel, percentiles, sun_up), fc))
})



