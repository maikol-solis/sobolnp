#' Nonparametric Sobol Estimator with Bootstrap Bandwidth
#'
#' Algorithm to estimate the Sobol indices using a non-parametric
#' fit of the regression curve. The bandwidth is estimated using
#' bootstrap to reduce the finite-sample bias.
#'
#' @param Y Response continuous variable
#' @param X Matrix of independent variables
#' @param bandwidth If \code{bandwidth.compute = TRUE}, it sets the starting bandwidth to find the bootstrap bandwidth. If \code{NULL} the least-square cross-validation bandwidth is used.
#' If \code{bandwidth.compute = FALSE}, it will use the values provided fixed over all the simulation. Defaults to \code{NULL}.
#' @param bandwidth.compute Logical value. Indicates if the bandwidth should be estimated or not. Defaults to \code{TRUE}.
#' @param bootstrap Logical value. Indicates if the estimation should be with bootstrap or not. Defaults to \code{TRUE}.
#' @param nboot Number of bootstrap samples taken for the method. Ignored if `bootstrap = FALSE`. Defaults to \code{100}.
#' @param ckerorder Numeric value specifying kernel order (should be one of
#' \code{(2,4,6,8)}). Defaults to \code{2}.
#' @param mc.cores Number of cores used. Defaults to \code{1}.
#'
#' @return A list of class \code{sobolnp} with the following elements:
#' \describe{
#'   \item{\strong{S}}{First order Sobol indices estimated with nonparametric
#'   regression and a cross-validation bandwidth}
#'   \item{\strong{bws}}{Bandwidth estimated with cross-validation}
#'   \item{\strong{Sboot}}{First order Sobol indices estimated with
#'   nonparametric regression and a bootstrap bandwidth}
#'   \item{\strong{bwsboot}}{Bandwidth estimated with bootstrap}
#' }
#' @references
#'
#'Solís, Maikol. "Nonparametric estimation of the first order Sobol indices with bootstrap bandwidth." \emph{arXiv preprint arXiv:1803.03333} (2018).
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
#' estimation <- sobolnp(Y = Y, X = X)
#'
#' @import np pbmcapply minqa
#' @importFrom stats residuals var

sobolnp <- function(Y,
                    X,
                    bandwidth = NULL,
                    # bandwidth.total = NULL,
                    bandwidth.compute = TRUE,
                    bootstrap = TRUE,
                    nboot = 100,
                    ckerorder = 2,
                    mc.cores = 1){
  n <- length(Y)
  p <- ncol(X)
  output <- list(call = match.call(),
                 num.obs = n,
                 num.var = p)

  if (any(bandwidth < 0)) {
    return(list(S = rep(NA, p),
                #   ST = rep(NA, p),
                bws = rep(NA, p)
                #   bwsT = rep(NA, p)
                ))
  }

  if (is.null(bandwidth)) {
    bandwidth <- NULL
  }

  # if (is.null(bandwidth.total)) {
  #   bandwidth.total <- rep(NULL, p - 1)
  # }

  message("Estimating non-parametric Sobol indices with cross-validation bandwidth")
  SList <- pbmclapply(
    X = 1:p,
    FUN = function(k) {
      compute.sobol.indices(
        xdat = X[, k],
        ydat = Y,
        bws = bandwidth,
        bandwidth.compute = bandwidth.compute,
        ckerorder = ckerorder
      )
    },
    mc.cores = mc.cores
  )

  # STotalList <- pbmclapply(
  #   X = 1:p,
  #   FUN = function(k) {
  #     x <- compute.sobol.indices(
  #       xdat = X[, -k],
  #       ydat = Y,
  #       bws = bandwidth.total,
  #       bandwidth.compute = bandwidth.compute,
  #       ckerorder = ckerorder
  #     )
  #     x$S <- 1 - x$S
  #     x$bws <-
  #       append(x$bws , NA, after = k - 1)
  #     return(x)
  #   },
  #   mc.cores = mc.cores
  # )



  S <- vapply(SList, function(l)
    l$S)
  bws_initial <- vapply(SList, function(l)
    l$bws)
  # ST <- sapply(STotalList, function(l)
  #   l$S)
  # bwsT_initial <- t(sapply(STotalList, function(l)
  #   l$bws))

  if (bootstrap == FALSE) {
    return(
      list(
        S = S,
        # ST = ST,
        bws = bws_initial
        # bwsT = bwsT_initial
        ))
  } else {
    meanS <- vapply(SList, function(l)
      l$mean)
    sderr <- vapply(SList, function(l)
      l$sderr)
    # meanres <- sapply(SList, function(l)
    #   l$meanres)
    error <- vapply(SList, function(l)
      l$error)
    # meanT <- sapply(STotalList, function(l)
    #   l$mean)
    # sderrT <- sapply(STotalList, function(l)
    #   l$sderr)
    # errorT <- sapply(STotalList, function(l)
    #   l$error)

    Y_boot <- mapply(
      FUN = function(n) {
        #browser()
        idx <-
          matrix(sample(1:n, size = n * p, replace = TRUE), nrow = n)
        return(meanS + sderr * error[idx])
      },
      rep(n, nboot),
      SIMPLIFY = FALSE
    )

    # Y_bootT <- mapply(
    #   FUN = function(n) {
    #     #browser()
    #     idx <-
    #       matrix(sample(1:n, size = n * p, replace = TRUE), nrow = n)
    #     return(meanT + sderrT * errorT[idx])
    #   },
    #   rep(n, nboot),
    #   SIMPLIFY = FALSE
    # )

    # Sboot <- numeric(p)
    # STboot <- numeric(p)
    # bwsboot <- numeric(p)
    # bwsSTboot <- matrix(nrow = p, ncol = p)
    #
    message("Estimating non-parametric Sobol indices with bootstrap bandwidth")
    Sfitboot <-   pbmclapply(
      X = 1:p,
      FUN = function(k) {
        cv_opt <-
          minqa::bobyqa(
            par = bws_initial[k],
            fn = BLS,
            xdat = X[, k],
            ydat = Y,
            ydat_boot = Y_boot,
            var.index = k,
            ckerorder = ckerorder,
            #method = "L-BFGS-B"#,
            lower = 0,
            upper = 10 * bws_initial[k],
            # method = "L-BFGS-B",
            control = list(iprint = 0)
          )


        gS <- compute.sobol.indices(
          xdat = X[, k],
          ydat = Y,
          bws = cv_opt$par,
          bandwidth.compute = FALSE,
          ckerorder = ckerorder,
          only.mean = TRUE
        )
        # Sboot[k] <- var(gS) / var(Y)
        # bwsboot[k] <-  cv_opt$par

        return(list(Sboot = var(gS) / var(Y),
                    bwsboot = cv_opt$par))
      },
      mc.cores = mc.cores
    )

    Sboot <- vapply(Sfitboot, function(x) x$Sboot)
    bwsboot <- vapply(Sfitboot, function(x) x$bwsboot)
  # for (k in 1:p) {
  #     cv_opt <-
  #       minqa::bobyqa(
  #         par = bws_initial[k],
  #         fn = BLS,
  #         xdat = X[, k],
  #         ydat = Y,
  #         ydat_boot = Y_boot,
  #         var.index = k,
  #         ckerorder = ckerorder,
  #         #method = "L-BFGS-B"#,
  #         lower = 0,
  #         upper = 10 * bws_initial[k],
  #         # method = "L-BFGS-B",
  #         control = list(iprint = 2)
  #       )
  #
  #
  #     gS <- compute.sobol.indices(
  #       xdat = X[, k],
  #       ydat = Y ,
  #       bws = cv_opt$par,
  #       bandwidth.compute = FALSE,
  #       ckerorder = ckerorder ,
  #       only.mean = TRUE
  #     )
  #     Sboot[k] <- var(gS) / var(Y)
  #     bwsboot[k] <-  cv_opt$par
  #
  #
  #     #   cv_opt_T <-
  #     #     bobyqa(
  #     #       par = bwsT_initial[k, -k],
  #     #       fn = BLS,
  #     #       xdat = X[, -k],
  #     #       ydat = Y,
  #     #       ydat_boot = Y_bootT,
  #     #       var.index = k,
  #     #       ckerorder = ckerorder,
  #     #       #method = "L-BFGS-B"#,
  #     #       lower = 0,
  #     #       upper = 10 * bwsT_initial[k, -k],
  #     #       control = list(iprint = 2)
  #     #     )
  #     #
  #     #   gT <- compute.sobol.indices(
  #     #     xdat = X[, -k],
  #     #     ydat = Y ,
  #     #     bws = cv_opt_T$par,
  #     #     bandwidth.compute = FALSE,
  #     #     ckerorder = ckerorder ,
  #     #     only.mean = TRUE
  #     #   )
  #     #   STboot[k] <- 1 - var(gT) / var(Y)
  #     #
  #     #
  #     #   bwsSTboot[k, ] <-
  #     #     append(cv_opt_T$par, NA, after = k - 1)
  #   }
    names(S) <- names(Sboot) <- colnames(X)
    output[["S"]] <- S
    output[["bws"]] <- bws_initial
    output[["Sboot"]] <- Sboot
    output[["bwsboot"]] <- bwsboot

    class(output) <- "sobolnp"

    return(output)
  }#end-else-bootstrap
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
      return(mean = g)
    } else{
      bws <- ghat$bw

      res <- residuals(ghat)

      sderr <- ghat$merr
      sderr.idx <- sderr <= 1e-5
      #meanres <- mean(res)
      error  <- (res - mean(res)) / sderr
      error[sderr.idx] <- 0

      return(list(
        S = var(g) / var(ydat),
        mean = g,
        bws = bws,
        sderr = sderr,
        #meanres = meanres,
        error = error
      ))
    }
  }


BLS <-
  function(bws,
           xdat,
           ydat,
           ydat_boot,
           var.index,
           ckerorder) {
    #message(bws)

    if (any(bws < 0) | any(is.na(bws))) {
      return(10e10)
    }

    mean_boot <-  sapply(ydat_boot, function(yboot,
                                             xdat,
                                             bws,
                                             ckerorder) {
      g <- compute.sobol.indices(
        xdat = xdat,
        ydat = yboot[, var.index],
        bws = bws,
        bandwidth.compute = FALSE,
        ckerorder = ckerorder,
        only.mean = TRUE
      )
      return(g)
    },
    xdat = xdat,
    bws = bws,
    ckerorder)

    return(mean(abs(ydat - rowMeans(mean_boot))))

    # delta <- 10
    # yerr <- ydat - rowMeans(mean_boot)
    # loss <-
    #   ifelse(abs(yerr) <= delta, 0.5 * yerr ^ 2, delta * (abs(yerr) - delta *                                                           0.5))
    # print(mean(loss))
    #     return(mean(loss))

    # return(mean(log(cosh(
    #
    # ))))
    # return(mean((ydat - rowMeans(mean_boot)) ^ 2))
  }
