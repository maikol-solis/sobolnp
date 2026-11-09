#' Nonparametric Sobol Estimator with Bootstrap Bandwidth
#'
#' Algorithm to estimate the Sobol indices using a non-parametric
#' fit of the regression curve. The bandwidth is estimated using
#' bootstrap to reduce the finite-sample bias.
#'
#' @param Y Response continuous variable
#' @param X Numeric matrix or data frame of independent variables, with one
#'   row per element of \code{Y}. Unnamed columns are labelled
#'   \code{X1, ..., Xp}.
#' @param bandwidth \code{NULL} or a positive numeric vector of length 1 or
#'   \code{ncol(X)} (one bandwidth per variable). If
#'   \code{bandwidth.compute = TRUE}, it sets the starting value of the
#'   least-squares cross-validation search; if \code{NULL} np chooses the
#'   starting values. If \code{bandwidth.compute = FALSE}, the bandwidths are
#'   used as given and must be supplied. Defaults to \code{NULL}.
#' @param bandwidth.compute Logical value. Indicates if the bandwidth should be estimated or not. Defaults to \code{TRUE}.
#' @param bootstrap Logical value. Indicates if the estimation should be with bootstrap or not. Defaults to \code{TRUE}.
#' @param nboot Number of bootstrap samples taken for the method. Ignored if `bootstrap = FALSE`. Defaults to \code{100}.
#' @param ckerorder Numeric value specifying kernel order (should be one of
#' \code{(2,4,6,8)}). Defaults to \code{2}.
#' @param mc.cores Number of cores used. Defaults to \code{1}.
#' @param verbose Logical value. If \code{TRUE} (default), report progress
#'   messages and progress bars.
#'
#' @details
#' For each input \eqn{X_i}, the regression curve \eqn{E[Y | X_i]} is
#' estimated by kernel regression with a least-squares cross-validation (CV)
#' bandwidth, and \code{S} is the variance of the fitted curve divided by the
#' variance of \code{Y}.
#'
#' With \code{bootstrap = TRUE}, the residuals of the CV fit are standardized
#' with a kernel estimate of their scale and resampled to build \code{nboot}
#' bootstrap responses. For a bandwidth \eqn{h}, the bootstrap curves give a
#' bias-corrected curve \eqn{2 \hat m_{CV} - \bar m^*_h} and a bootstrap
#' index. The bootstrap bandwidth \code{bwsboot} minimizes the mean squared
#' distance between the bootstrap indices and \code{S}, searching between
#' \code{bws / 1000} and \code{bws}, and \code{Sboot} is the mean bootstrap
#' index at that bandwidth.
#'
#' @return A list of class \code{sobolnp} with the following elements:
#' \describe{
#'   \item{\strong{call}}{The matched call.}
#'   \item{\strong{num.obs}, \strong{num.var}}{Number of observations and of
#'   input variables.}
#'   \item{\strong{S}}{First order Sobol indices estimated with nonparametric
#'   regression and a cross-validation bandwidth}
#'   \item{\strong{bws}}{Bandwidth estimated with cross-validation}
#'   \item{\strong{Sboot}}{First order Sobol indices estimated with
#'   nonparametric regression and a bootstrap bandwidth}
#'   \item{\strong{bwsboot}}{Bandwidth estimated with bootstrap}
#' }
#' \code{Sboot} and \code{bwsboot} are \code{NULL} when
#' \code{bootstrap = FALSE}.
#' @references
#'
#' Solís, M. (2021). Non-parametric estimation of the first-order Sobol indices with bootstrap bandwidth. \emph{Communications in Statistics - Simulation and Computation}, 50(9), 2497-2512. \doi{10.1080/03610918.2019.1655575}
#'
#' @export
#' @examples
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
#' @importFrom np npreg npregbw
#' @importFrom pbmcapply pbmclapply
#' @importFrom stats optimize residuals var
#' @importFrom parallel mclapply

sobolnp <- function(Y,
                    X,
                    bandwidth = NULL,
                    bandwidth.compute = TRUE,
                    bootstrap = TRUE,
                    nboot = 100,
                    ckerorder = 2,
                    mc.cores = 1,
                    verbose = TRUE){
  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }
  if (is.null(dim(X))) {
    X <- matrix(X, ncol = 1)
  }
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("'X' and 'Y' must be numeric.", call. = FALSE)
  }
  Y <- as.vector(Y)
  n <- length(Y)
  p <- ncol(X)
  if (nrow(X) != n) {
    stop("'Y' must have one value per row of 'X'.", call. = FALSE)
  }
  if (anyNA(X) || anyNA(Y)) {
    stop("'X' and 'Y' must not contain missing values.", call. = FALSE)
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(p))
  }
  if (length(ckerorder) != 1 || !ckerorder %in% c(2, 4, 6, 8)) {
    stop("'ckerorder' must be one of 2, 4, 6 or 8.", call. = FALSE)
  }
  if (bootstrap &&
      (length(nboot) != 1 || !is.numeric(nboot) || nboot < 1 ||
       nboot != round(nboot))) {
    stop("'nboot' must be a positive integer.", call. = FALSE)
  }
  if (!is.null(bandwidth)) {
    if (!length(bandwidth) %in% c(1, p) || anyNA(bandwidth) ||
        any(bandwidth <= 0)) {
      stop("'bandwidth' must be NULL or positive, of length 1 or ncol(X).",
           call. = FALSE)
    }
    bandwidth <- rep_len(bandwidth, p)
  } else if (!bandwidth.compute) {
    stop("'bandwidth' must be supplied when 'bandwidth.compute = FALSE'.",
         call. = FALSE)
  }

  old_options <- options(np.messages = verbose)
  on.exit(options(old_options), add = TRUE)
  lapply_cores <- if (verbose) pbmclapply else mclapply
  say <- function(...) if (verbose) message(...)

  output <- list(call = match.call(),
                 num.obs = n,
                 num.var = p)

  say("Estimating non-parametric Sobol indices with cross-validation bandwidth")
  SList <- lapply_cores(
    X = 1:p,
    FUN = function(k) {
      compute.sobol.indices(
        xdat = X[, k],
        ydat = Y,
        bws = bandwidth[k],
        bandwidth.compute = bandwidth.compute,
        ckerorder = ckerorder
      )
    },
    mc.cores = mc.cores
  )

  S <- sapply(SList, function(l)
    l$S)
  bws_initial <- sapply(SList, function(l)
    l$bws)

  names(S) <- names(bws_initial) <- colnames(X)
  output[["S"]] <- S
  output[["bws"]] <- bws_initial

  if (bootstrap) {
    # Bootstrap row indices are drawn once, in the main process, so results
    # are reproducible with set.seed() for any number of cores.
    boot_idx <- replicate(nboot, sample.int(n, n, replace = TRUE))
    if (!is.matrix(boot_idx)) {
      boot_idx <- matrix(boot_idx, nrow = n)
    }

    say("Estimating non-parametric Sobol indices with bootstrap bandwidth")
    Sfitboot <- lapply_cores(
      X = 1:p,
      FUN = function(k) {
        bootstrap.sobol.index(
          xdat = X[, k],
          ydat = Y,
          fit = SList[[k]],
          boot_idx = boot_idx,
          ckerorder = ckerorder
        )
      },
      mc.cores = mc.cores
    )

    Sboot <- sapply(Sfitboot, function(x) x$Sboot)
    bwsboot <- sapply(Sfitboot, function(x) x$bwsboot)
    names(Sboot) <- names(bwsboot) <- colnames(X)
    output[["Sboot"]] <- Sboot
    output[["bwsboot"]] <- bwsboot
  }

  class(output) <- "sobolnp"
  output
}

compute.sobol.indices <-
  function(xdat,
           ydat,
           bws,
           bandwidth.compute,
           ckerorder,
           only.mean = FALSE) {
    ghat <- npreg(
      npregbw(
        bws = bws,
        xdat = xdat,
        ydat = ydat,
        bandwidth.compute = bandwidth.compute,
        ckertype = "epanechnikov",
        ckerorder = ckerorder,
        bwmethod = "cv.ls"
      ),
      residuals = TRUE
    )
    g <- ghat$mean
    if (only.mean) {
      return(g)
    } else {
      return(list(
        S = var(g) / var(ydat),
        mean = g,
        bws = ghat$bw,
        residuals = residuals(ghat)
      ))
    }
  }

bootstrap.sobol.index <-
  function(xdat,
           ydat,
           fit,
           boot_idx,
           ckerorder) {
    # Scale of the residuals: standard error of a kernel regression of the
    # residuals on xdat (standard errors are opt-in since np 0.80-1).
    res <- fit$residuals
    scale_fit <- npreg(
      npregbw(
        xdat = xdat,
        ydat = res,
        ckertype = "epanechnikov",
        bwmethod = "cv.ls"
      ),
      se = TRUE
    )
    sderr <- scale_fit$merr
    error <- (res - mean(res)) / sderr
    error[sderr <= 1e-5] <- 0

    # Bootstrap responses: CV curve plus rescaled, resampled residuals.
    y_boot <- apply(boot_idx, 2, function(i) fit$mean + sderr * error[i])
    var_y_boot <- apply(y_boot, 2, var)

    # Bootstrap indices from the bias-corrected curve 2 * m_CV - mean(m*_h).
    index_boot <- function(h) {
      m_boot <- apply(y_boot, 2, function(yb) {
        compute.sobol.indices(
          xdat = xdat,
          ydat = yb,
          bws = h,
          bandwidth.compute = FALSE,
          ckerorder = ckerorder,
          only.mean = TRUE
        )
      })
      m_corrected <- 2 * fit$mean - rowMeans(m_boot)
      var(m_corrected) / var_y_boot
    }

    criterion <- function(h) {
      value <- mean((index_boot(h) - fit$S)^2)
      if (is.finite(value)) value else .Machine$double.xmax
    }

    h_cv <- fit$bws
    h_boot <- optimize(criterion, lower = h_cv / 1000, upper = h_cv)$minimum

    list(Sboot = mean(index_boot(h_boot)), bwsboot = h_boot)
  }
