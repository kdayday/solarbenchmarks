library(testthat)

# Currently set up to run from the Benchmarks_comparison folder
source("evaluation_functions.R")

test_that("QS function throws error", {
  expect_error(QS(matrix(1:6, ncol=2), tel=1:4, sun_up=c(T,T,T,T), percentiles=NA), "*range and resolution")
  expect_error(QS(matrix(1:8, ncol=2), tel=1:4, sun_up=c(T,T,T), percentiles=NA), "*range and resolution")
})

test_that("Quantile score is as expected", {
  fc <- rbind(c(20, 50, 70),
              c(20, 50, 70),
              c(10, 25, 35),
              c(10, 25, 35))
  expect_equal(QS(fc, tel=c(0, 40, 5, NA), sun_up=c(F,T,T,T), percentiles=c(0.2, 0.5, 0.7)), c(8, 15 , 18)) # 2*mean(0.2*20, 0.8*5)=8, (0.5*10, 0.5*20), (0.3*30, 0.3*30)
})

test_that("Weighted QS function throws error", {
  expect_error(weight_QS(qs=c(1,1), percentiles=c(0.2, 0.5, 0.7), weighting="none"), "*same length")
})

test_that("Weighted QS with no weight is correct", {
  expect_equal(weight_QS(qs=c(1,1,1), percentiles=c(0.2, 0.5, 0.7), weighting="none"), c(1,1,1))
})

test_that("Left-tail weighted quantile scores are correct", {
  expect_equal(weight_QS(qs=c(1,1,1), percentiles=c(0.2, 0.5, 0.7), weighting="left"), c(0.64, 0.25, 0.09))
})

test_that("Right-tail weighted quantile scores are correct", {
  expect_equal(weight_QS(qs=c(1,1,1), percentiles=c(0.2, 0.5, 0.7), weighting="right"), c(0.04, 0.25, 0.49))
})

test_that("Center weighted quantile scores are correct", {
  expect_equal(weight_QS(qs=c(1,1,1), percentiles=c(0.2, 0.5, 0.7), weighting="center"), c(0.16, 0.25, 0.21))
})

test_that("Tails weighted quantile scores are correct", {
  expect_equal(weight_QS(qs=c(1,1,1), percentiles=c(0.2, 0.5, 0.7), weighting="tails"), c(0.36, 0, 0.16))
})

test_that("Average interval width calculation is correct", {
  fc <- rbind(c(1, 3, 5, 6),
              c(1, 17, 36, 100),
              c(1, 3, 7, 8))
  colnames(fc) <- c(0.1, 0.4, 0.6, 0.9)
  sun_up <- c(T, F, T)
  out <- interval_width(fc, sun_up, intervals=c(0.2, 0.8))
  expect_equal(out$widths, c(3, 6))
  expect_equal(out$intervals, c(0.2, 0.8))
})

