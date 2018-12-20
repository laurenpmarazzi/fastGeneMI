context("Testing network inference")

test_that("Non-symmetric mim", {
  expect_error(infer.net(matrix(runif(10*10),10,10), "clr"))
})

test_that("Non-square mim test", {
  expect_error(infer.net(matrix(runif(5*10),5,10), "clr"))
})

test_that("Non-symmetric gold standard", {
  expect_error(infer.net(diag(10), "clr", matrix(sample(0:1, 10*10, replace=T), 10, 10)))
})

test_that("Mismatching mim and gold standard", {
  expect_error(infer.net(diag(10), "clr", matrix(sample(0:1, 11*11, replace=T), 11, 11)))
})

test_that("Non-binary gold standard", {
  gs.net <- matrix(sample(0:2, 10*10, replace=T), 10, 10)
  expect_error(infer.net(diag(10), "clr", gs.net + t(gs.net)))
})