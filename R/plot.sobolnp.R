#' Plot method for objects \code{sobolnp}
#'
#' Plot the Sobol indices based in a non-parametric regression
#' with cross-validation and bootstrap bandwidth
#'
#' @param snp an object of class \code{sobolnp}
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
#' X <- matrix(runif(3*1000, -pi,pi), ncol = 3)
#' Y <- ishigami.fun(X)
#'
#' estimation <- sobolnp::sobolnp(Y = Y,X = X, mc.cores = 1)
#'
#' plot(estimation)
#' }
#'
#' @importFrom  graphics plot.default axis legend segments

plot <- function(snp, ...) {
  UseMethod("plot", snp)
}

#' @export
#' @rdname plot
plot.sobolnp <- function(snp, ...) {

  yrang <- c(0, max(c(snp[["Sboot"]], snp[["S"]])) + 0.1)
  xrang <- range(1:snp[["num.var"]])
  dx <- diff(xrang) / (3 * snp[["num.var"]])

  xrang[1] <- xrang[1] - dx
  xrang[2] <- xrang[2] + dx

  graphics::plot.default(
    x = 1:snp[["num.var"]],
    xlim = xrang,
    ylim = yrang,
    xaxt = "n",
    ylab = "Sobol  indices",
    xlab = "Variables",
    ...
  )

  segments(
    x0 = (1:snp[["num.var"]]) - dx,
    x1 = (1:snp[["num.var"]]) + dx,
    y0 = snp[["S"]],
    lwd = 3,
    col = "blue"
  )

  segments(
    x0 = (1:snp[["num.var"]]) - dx,
    x1 = (1:snp[["num.var"]]) + dx,
    y0 = snp[["Sboot"]],
    lwd = 3,
    col = "red"
  )

  axis(side = 1,
       at = 1:snp[["num.var"]],
       labels = names(snp[["S"]]))

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
