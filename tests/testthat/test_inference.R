context("Testing network inference")

test_that("Non-symmetric mim error", {
  expect_error(infer.net(matrix(runif(10*10),10,10), "clr"))
})

test_that("Non-square mim test error", {
  expect_error(infer.net(matrix(runif(5*10),5,10), "clr"))
})

test_that("Non-symmetric gold standard error", {
  expect_error(infer.net(diag(10), "clr", matrix(sample(0:1, 10*10, replace=T), 10, 10)))
})

test_that("Mismatching mim and gold standard error", {
  expect_error(infer.net(diag(10), "clr", matrix(sample(0:1, 11*11, replace=T), 11, 11)))
})

test_that("Non-binary gold standard error", {
  gs.net <- matrix(sample(0:2, 10*10, replace=T), 10, 10)
  expect_error(infer.net(diag(10), "clr", gs.net + t(gs.net)))
})

test_that("Inference size error and regulator warning", {
  mim <- get.mim.ML(matrix(runif(100), 10, 10))
  gs.net <- matrix(sample(0:1, 10*10, replace=T), 10, 10)
  gs.net[upper.tri(gs.net)] <- 0
  gs.net <- gs.net + t(gs.net)
  diag(gs.net) <- 0
  out <- infer.net(mim, gs.net=gs.net, n.reg=4, plot.pr = FALSE)
  expect_equal(nrow(out$network), ncol(out$network))
  expect_equal(nrow(out$network), nrow(mim))
  expect_equal(ncol(out$pr), 3)
  expect_warning(infer.net(mim, gs.net=gs.net))
})