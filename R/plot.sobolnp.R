#' Plot method for objects \code{sobolnp}
#'
#' Plot the Sobol indices based in a non-parametric regression
#' with cross-validation and bootstrap bandwidth
#'
#' @param x an object of class \code{sobolnp}
#' @param ... further arguments passed to the \code{plot} function
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
#' plot(estimation)
#' }
#'
#' @importFrom  graphics axis legend segments

plot <- function(x, ...) {
  UseMethod("plot", x)
}

#' @export
#' @rdname plot
plot.sobolnp <- function(x, ...) {

  yrang <- c(0, max(c(x[["Sboot"]], x[["S"]])) + 0.1)
  xrang <- range(1:x[["num.var"]])
  dx <- diff(xrang) / (3 * x[["num.var"]])

  xrang[1] <- xrang[1] - dx
  xrang[2] <- xrang[2] + dx

  plot(
    x = 1:x[["num.var"]],
    xlim = xrang,
    ylim = yrang,
    xaxt = "n",
    ylab = "Sobol  indices",
    xlab = "Variables",
    ...
  )

  segments(
    x0 = 1:x[["num.var"]] - dx,
    x1 = 1:x[["num.var"]] + dx,
    y0 = x[["S"]],
    lwd = 3,
    col = "blue"
  )

  segments(
    x0 = 1:x[["num.var"]] - dx,
    x1 = 1:x[["num.var"]] + dx,
    y0 = x[["Sboot"]],
    lwd = 3,
    col = "red"
  )

  axis(side = 1,
       at = 1:x[["num.var"]],
       labels = names(x[["S"]]))

  legend(
    "topright",
    legend = c("CV", "Bootstrap"),
    lty = 1,
    lwd = 3,
    col = c("blue", "red"),
    horiz = FALSE,
    bty = "n"
  )

}
