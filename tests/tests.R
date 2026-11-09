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
  nboot = 2,
  mc.cores = 1
)

print(ss)

plot(ss)

# Two input variables (previously failed in the residual bootstrap)
ss2 <- sobolnp::sobolnp(Y = Y, X = X[, 1:2], nboot = 2, mc.cores = 1)
stopifnot(length(ss2$Sboot) == 2)

# Without bootstrap the result is still a sobolnp object
ss_cv <- sobolnp::sobolnp(Y = Y, X = X, bootstrap = FALSE)
stopifnot(inherits(ss_cv, "sobolnp"), is.null(ss_cv$Sboot),
          identical(names(ss_cv$S), c("X1", "X2", "X3")))
print(ss_cv)
plot(ss_cv)

# One fixed bandwidth per variable
ss_bw <- sobolnp::sobolnp(Y = Y, X = X, bandwidth = c(1, 1, 1),
                          bandwidth.compute = FALSE, bootstrap = FALSE)
stopifnot(all(ss_bw$bws == 1))

# Invalid input is rejected
stopifnot(inherits(try(sobolnp::sobolnp(Y = Y[-1], X = X), silent = TRUE),
                   "try-error"))
stopifnot(inherits(try(sobolnp::sobolnp(Y = Y, X = X, bandwidth = -1),
                       silent = TRUE), "try-error"))
