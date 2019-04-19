library(sobolnp)

ishigami.fun <- function(X) {
  A <- 7
  B <- 0.1
  sin(X[, 1]) + A * sin(X[, 2]) ^ 2 + B * X[, 3] ^ 4 * sin(X[, 1])
}

X <- matrix(runif(3 * 10, -pi, pi), ncol = 3)
Y <- ishigami.fun(X)

#colnames(X) <- c("1", "2", "3")
ss <- sobolnp::sobolnp(
  Y = Y,
  X = X,
  nboot = 10,
  mc.cores = 1
)

print(ss)

plot(ss)
