test_that("Model w intercept works.", {
  set.seed(123)
  x <- matrix(1:250, 50)
  storage.mode(x) <- "double"
  y <- 1:50 + rnorm(50)

  res <- reg.sgl(x, y, gindex = c(1, 1, 1, 2, 2), nlambda = 3)
  res$set.seed <- NULL

  expect_snapshot(res)
})

test_that("Model w/o intercept works.", {
  set.seed(123)
  x <- matrix(1:250, 50)
  y <- 1:50 + rnorm(50)

  res <- reg.sgl(x, y, gindex = c(1, 1, 1, 2, 2), nlambda = 3, intercept = F)
  res$set.seed <- NULL

  expect_snapshot(res)
})
