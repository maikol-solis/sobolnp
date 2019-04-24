#' Print method for objects \code{sobolnp}
#'
#' @param snp an object of class \code{sobolnp}
#' @param ... further arguments passed to the \code{print} function
#'
#' @return
#' A formatted table with the results of the \code{\link{sobolnp}}
#' function.
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

print <- function(snp, ...) {
  UseMethod("print", snp)
}

#' @export
#' @rdname print
print.sobolnp <- function(snp, ...) {
  cat("\nCall:\n", deparse(snp[["call"]]), "\n", sep = "")
  cat("\nNumber of observations:", snp[["num.obs"]], "\n")
  cat("\nNumber of variables:", snp[["num.var"]], "\n")
  cat("\nFirst order indices\n")

  m <- cbind(snp[["S"]], snp[["bws"]], snp[["Sboot"]], snp[["bwsboot"]])
  colnames(m) <- c("Si CV", "bw CV", "Si Bootstrap", "bw Bootstrap")
  rownames(m) <- names(snp[["S"]])
  print(m, ...)
}
