context("Testing mutual information estimators")

data("iris")
n.bins <- 10

test_that("Testing ML estimator", {

  run.tests <- function(disc.method) {
    mim.fgmi <- get.mim.ML(iris[,1:4], disc.method, n.bins=n.bins)
    expect_identical(is.matrix(mim.fgmi), TRUE)
    expect_identical(isSymmetric(mim.fgmi), TRUE)
    expect_equal(nrow(mim.fgmi), 4)
    iris.disc <- infotheo::discretize(iris[,1:4], disc.method, nbins = n.bins)
    mim.it <- infotheo::mutinformation(iris.disc, method = "emp")
    expect_equivalent(mim.fgmi, mim.it)
  }
  
  lapply(c("equalfreq", "equalwidth"), run.tests)
  
})

test_that("Testing MM estimator", {
  
  run.tests <- function(disc.method) {
    mim.fgmi <- get.mim.MM(iris[,1:4], disc.method, n.bins=n.bins)
    expect_identical(is.matrix(mim.fgmi), TRUE)
    expect_identical(isSymmetric(mim.fgmi), TRUE)
    expect_equal(nrow(mim.fgmi), 4)
    iris.disc <- infotheo::discretize(iris[,1:4], disc.method, nbins = n.bins)
    mim.it <- infotheo::mutinformation(iris.disc, method = "mm")
    expect_equivalent(mim.fgmi, mim.it)
  }
  
  lapply(c("equalfreq", "equalwidth"), run.tests)
  
})

test_that("Testing CS estimator", {
  
  run.tests <- function(disc.method) {
    mim.fgmi <- get.mim.CS(iris[,1:4], disc.method, n.bins=n.bins)
    expect_identical(is.matrix(mim.fgmi), TRUE)
    expect_identical(isSymmetric(mim.fgmi), TRUE)
    expect_equal(nrow(mim.fgmi), 4)
  }
  
  lapply(c("equalfreq", "equalwidth"), run.tests)
  
})

test_that("Testing Shrinkage estimator", {
  
  run.tests <- function(disc.method) {
    mim.fgmi <- get.mim.shrink(iris[,1:4], disc.method, n.bins=n.bins)
    expect_identical(is.matrix(mim.fgmi), TRUE)
    expect_identical(isSymmetric(mim.fgmi), TRUE)
    expect_equal(nrow(mim.fgmi), 4)
    iris.disc <- infotheo::discretize(iris[,1:4], disc.method, nbins = n.bins)
    mim.it <- infotheo::mutinformation(iris.disc, method = "shrink")
    #expect_equivalent(mim.fgmi, mim.it)
  }
  
  lapply(c("equalfreq", "equalwidth"), run.tests)
  
})

test_that("Testing B-spline estiamtor", {
  mim.bspline <- get.mim.bspline(iris[,1:4], order = 1, n.bins = n.bins)
  mim.ml <- get.mim.ML(iris[,1:4], discretisation = "equalwidth", n.bins = n.bins)
  expect_identical(is.matrix(mim.bspline), TRUE)
  expect_identical(isSymmetric(mim.bspline), TRUE)
  expect_equal(nrow(mim.bspline), 4)
  expect_equal(mim.bspline, mim.ml, 1e-2)
})