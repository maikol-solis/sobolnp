#' Plot method for objects \code{sobolnp}
#'
#' Bar plot of the first order Sobol indices estimated by
#' \code{\link{sobolnp}}, with the cross-validation and the bootstrap
#' bandwidth side by side for each variable.
#'
#' @param x an object of class \code{sobolnp}.
#' @param ylim limits of the y axis. By default from 0 to the larger of 1 and
#'   the largest estimated index.
#' @param xlab,ylab axis labels.
#' @param col colors for the cross-validation and bootstrap bars.
#' @param legend.pos position of the legend, passed to
#'   \code{\link[graphics]{legend}}. Use \code{NULL} to omit it.
#' @param ... further graphical arguments passed to
#'   \code{\link[graphics]{barplot}}.
#'
#' @return The object \code{x}, invisibly. Called for its side effect of
#'   drawing a plot.
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
#' plot(estimation)
#'
#' @importFrom graphics barplot legend
plot.sobolnp <- function(x,
                         ylim = NULL,
                         xlab = "Variables",
                         ylab = "First order Sobol indices",
                         col = c("#0072B2", "#E69F00"),
                         legend.pos = "topright",
                         ...) {
  indices <- rbind(CV = x[["S"]], Bootstrap = x[["Sboot"]])
  if (is.null(colnames(indices))) {
    colnames(indices) <- paste0("X", seq_len(ncol(indices)))
  }
  col <- rep_len(col, nrow(indices))

  if (is.null(ylim)) {
    ylim <- c(0, 1.1 * max(1, indices, na.rm = TRUE))
  }

  barplot(
    indices,
    beside = TRUE,
    col = col,
    ylim = ylim,
    xlab = xlab,
    ylab = ylab,
    ...
  )

  if (!is.null(legend.pos)) {
    legend(
      legend.pos,
      legend = rownames(indices),
      fill = col,
      horiz = TRUE,
      bty = "n"
    )
  }

  invisible(x)
}
