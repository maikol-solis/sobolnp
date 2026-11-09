#' Print method for objects \code{sobolnp}
#'
#' @param x an object of class \code{sobolnp}.
#' @param digits minimal number of significant digits, see
#'   \code{\link{print.default}}.
#' @param ... further arguments passed to \code{\link{print.default}}.
#'
#' @return The object \code{x}, invisibly. Called for its side effect of
#'   printing the estimated indices and bandwidths.
#' @export
#'
#' @examples
#'
#' ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }
#'
#' X <- matrix(runif(3*100, -pi, pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- sobolnp(Y = Y, X = X, nboot = 5)
#'
#' print(estimation)
#'
print.sobolnp <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x[["call"]]), collapse = "\n"), "\n", sep = "")
  cat("\nNumber of observations:", x[["num.obs"]], "\n")
  cat("Number of variables:", x[["num.var"]], "\n")
  cat("\nFirst order indices\n")

  m <- cbind(x[["S"]], x[["bws"]], x[["Sboot"]], x[["bwsboot"]])
  colnames(m) <- c("Si CV", "bw CV", "Si Bootstrap", "bw Bootstrap")[seq_len(ncol(m))]
  rownames(m) <- names(x[["S"]])
  print.default(m, digits = digits, ...)
  invisible(x)
}
