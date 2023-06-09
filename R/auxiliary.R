chkvars2 <- function (x) {
  apply(x, 2, function (y) !all(y == y[1]))
}

standard2 <- function (x, isd, intr) {
  if (intr) {
    xmean <- colMeans(x)
  } else {
    xmean <- FALSE
  }

  if (isd) {
    xnorm <- apply(x, 2, function (y) sqrt(mean((y - mean(y))^2)))
  } else {
    xnorm <- FALSE
  }

  x <- scale(x, center = xmean, scale = xnorm)
  maj <- colMeans(x^2)

  return(list(x = x, maj = maj))
}