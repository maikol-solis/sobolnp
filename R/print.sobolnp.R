#' Print method for objects \code{sobolnp}
#'
#' @param x an object of class \code{sobolnp}
#' @param ... further arguments passed to the \code{print} function
#'
#' @return
#' A formatted table with the results of the \code{\link{sobolnp}}
#' function.
#' @export
#'
#' @examples
#'
#' \dontrun{
#' ishigami.fun <- function(X) {
#' A <- 7
#' B <- 0.1
#' sin(X[, 1]) + A * sin(X[, 2])^2 + B * X[, 3]^4 * sin(X[, 1])
#' }
#'
#' X<- matrix(runif(3*10, -pi,pi), ncol = 3)
#' Y<- ishigami.fun(X)
#'
#' estimation <- sobolnp::sobolnp(Y = Y,X = X, mc.cores = 1)
#'
#' print(estimation)
#' }
#'

print <- function(x,...){
  UseMethod("print", x)
}

#' @export
#' @rdname print
print.sobolnp <- function(x, ...) {
  cat("\nCall:\n", deparse(x[["call"]]), "\n", sep = "")
  cat("\nNumber of observations:", x[["num.obs"]], "\n")
  cat("\nNumber of variables:", x[["num.var"]], "\n")
  cat("\nFirst order indices\n")

  m <- cbind(x[["S"]], x[["bws"]], x[["Sboot"]], x[["bwsboot"]])
  colnames(m) <- c("Si CV", "bw CV", "Si Bootstrap", "bw Bootstrap")
  rownames(m) <- names(x[["S"]])
  print(m, ...)
}
