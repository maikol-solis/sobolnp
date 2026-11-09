ishigami <- function(X) {
  sin(X[, 1]) + 7 * sin(X[, 2])^2 + 0.1 * X[, 3]^4 * sin(X[, 1])
}

set.seed(123)
X <- matrix(runif(3 * 20, -pi, pi), ncol = 3)
Y <- ishigami(X)

fit <- sobolnp(Y = Y, X = X, nboot = 2, verbose = FALSE)

test_that("sobolnp returns a well-formed object", {
  expect_s3_class(fit, "sobolnp")
  expect_equal(fit$num.obs, 20)
  expect_equal(fit$num.var, 3)
  for (el in c("S", "bws", "Sboot", "bwsboot")) {
    expect_length(fit[[el]], 3)
    expect_named(fit[[el]], c("X1", "X2", "X3"))
  }
  expect_true(all(fit$S >= 0 & fit$Sboot >= 0))
  expect_true(all(fit$bws > 0 & fit$bwsboot > 0))
  expect_true(all(fit$bwsboot <= fit$bws * (1 + 1e-8)))
})

test_that("bootstrap results are reproducible with set.seed", {
  set.seed(1)
  a <- sobolnp(Y = Y, X = X, nboot = 2, verbose = FALSE)
  set.seed(1)
  b <- sobolnp(Y = Y, X = X, nboot = 2, verbose = FALSE)
  expect_identical(a$Sboot, b$Sboot)
})

test_that("two input variables work (residual bootstrap indexing)", {
  fit2 <- sobolnp(Y = Y, X = X[, 1:2], nboot = 2, verbose = FALSE)
  expect_length(fit2$Sboot, 2)
})

test_that("bootstrap = FALSE returns a sobolnp object", {
  fit_cv <- sobolnp(Y = Y, X = X, bootstrap = FALSE, verbose = FALSE)
  expect_s3_class(fit_cv, "sobolnp")
  expect_null(fit_cv$Sboot)
  expect_equal(fit_cv$S, fit$S)
})

test_that("column names and data frames are respected", {
  df <- as.data.frame(X)
  names(df) <- c("a", "b", "c")
  fit_df <- sobolnp(Y = Y, X = df, bootstrap = FALSE, verbose = FALSE)
  expect_named(fit_df$S, c("a", "b", "c"))
})

test_that("a fixed bandwidth per variable is used as given", {
  fit_bw <- sobolnp(Y = Y, X = X, bandwidth = c(1, 2, 3),
                    bandwidth.compute = FALSE, bootstrap = FALSE,
                    verbose = FALSE)
  expect_equal(unname(fit_bw$bws), c(1, 2, 3))
})

test_that("invalid input is rejected", {
  expect_error(sobolnp(Y = Y[-1], X = X), "one value per row")
  expect_error(sobolnp(Y = Y, X = X, bandwidth = -1), "positive")
  expect_error(sobolnp(Y = Y, X = X, bandwidth = c(1, 1)), "length")
  expect_error(sobolnp(Y = Y, X = X, bandwidth.compute = FALSE),
               "must be supplied")
  expect_error(sobolnp(Y = Y, X = X, ckerorder = 3), "ckerorder")
  expect_error(sobolnp(Y = Y, X = X, nboot = 0), "nboot")
  X_na <- X
  X_na[1, 1] <- NA
  expect_error(sobolnp(Y = Y, X = X_na), "missing")
})

test_that("verbose = FALSE is silent", {
  expect_silent(sobolnp(Y = Y, X = X, bootstrap = FALSE, verbose = FALSE))
})

test_that("print and plot return the object invisibly", {
  expect_output(out <- withVisible(print(fit)), "First order indices")
  expect_false(out$visible)
  expect_identical(out$value, fit)

  pdf(NULL)
  on.exit(dev.off())
  out <- withVisible(plot(fit, xlab = "custom", main = "title"))
  expect_false(out$visible)
  expect_identical(out$value, fit)
  expect_silent(plot(sobolnp(Y = Y, X = X[, 1, drop = FALSE],
                             bootstrap = FALSE, verbose = FALSE)))
})

test_that("sobolnp does not mask base generics", {
  expect_identical(print, base::print)
  expect_identical(plot, graphics::plot)
})
