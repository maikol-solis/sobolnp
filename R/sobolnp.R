# Author: Maikol Solís
#To create a robust function, I took as a template the script of Filippo Monari,
#'sobolSmthSpl' from the package 'sensitivity'. However, the core estimations
#are my original work.

sobolnp <- function(Y,
                    X,
                    nboot = 100,
                    bootstrap = FALSE,
                    ckerorder = 2) {
  #Arguments:
  #Y: matrix of model outputs (only one column)
  #X: matrix model parameters
  #nboot: number of bootstrap samples
  #bootstrap: logical if apply bootstrap or not
  #ckerorder: Kernel order for the nonparametric estimator
  #MAIN

  if (!requireNamespace("np", quietly = TRUE)) {
    stop("Please install the package np: install.packages('np')")
  }#end-require-np

  ANS = list()
  ANS[['call']] = match.call()
  ANS[['X']] = X
  ANS[['Y']] = Y
  ANS[['bootstrap']] <- bootstrap
  if (bootstrap) {
    ANS[['nboot']] = nboot
  }
  ANS[['ckerorder']] <- ckerorder
  par.names = colnames(X)	#gets parameters names
  if (is.null(colnames(X)))
    par.names = paste0('X', 1:ncol(X))
  # X = normalize(X)    		#normalize inpiuts between [0,1]
  #Y = Y - mean(Y)     		#center model responses
  Y = sapply(1:ncol(X), function(i)
    return(Y[order(X[, i])]))    #order Y before X (or create a new variable for X)
  X = sapply(1:ncol(X), function(i)
    return(X[order(X[, i]), i]))
  NP = optNP(Y,
             X,
             nboot = nboot,
             bootstrap = bootstrap,
             ckerorder = ckerorder)
  SA.tab = t(sapply(NP, est.Si))
  colnames(SA.tab) = c('Si', 'bw' , 'se', 'q0.05')
  rownames(SA.tab) = par.names
  ANS[['S']] = SA.tab
  class(ANS) = 'sobolnp'
  return(ANS)
}

est.Si <- function (NP) {
  gi <- NP[['mean']]
  yi <- NP[['y']]
  bw <- NP[['bws']]
  Si = var(gi) / var(yi)
  #calculates the standard error of the main effect estimates
  yi.sc = (yi - mean(yi)) / sd(yi)            	#scaled yi
  u = (yi - gi) / (sd(yi) * abs(1 - Si) ** 0.5)   #scales residuals
  Si.se = abs(1 - Si) * sd(yi.sc ** 2 - u ** 2) / length(yi) ** 0.5
  q0.05 = qnorm(0.05, Si, Si.se)
  return(c(Si, bw, Si.se, q0.05))
}

optNP <- function (Y, X, nboot, bootstrap, ckerorder) {
  #DOC
  #Estimates nonparametrically the regression curve across all the column of the
  #matrix Y.
  #Arguments:
  #Y: matrix containing the response variable
  #X: matrix of inputs
  #nboot: number of bootstrap samples
  #bootstrap: logical if apply bootstrap or not
  #ckerorder: Kernel order for the nonparametric estimator

  #CONTAINS

  ghat_boot <- function(Y_boot, Xi, h, ckerorder = 2) {
    ghat_boot <- lapply(
      X = Y_boot,
      FUN = function(ydat, xdat, h, ckerorder) {
        ghat.bw <- np::npregbw(
          xdat = xdat,
          ydat = ydat,
          bws = h,
          bandwidth.compute = FALSE,
          ckertype = "epanechnikov",
          ckerorder = ckerorder,
          bwmethod = "cv.ls"
        )
        ghat <- np::npreg(ghat.bw)
        return(mean = ghat$mean)
      },
      #end-function-ghat_boot
      xdat = Xi,
      h = h,
      ckerorder  = ckerorder
    )
    ghat_boot_mean <- apply(simplify2array(ghat_boot), 1, mean)
    return(ghat_boot_mean)
  }

  objective <- function(h,
                        Y_boot,
                        Y,
                        Xi,
                        ckerorder = 2) {
    ghat_boot_mean <- ghat_boot(
      Y_boot = Y_boot,
      Xi = Xi,
      h = h,
      ckerorder = ckerorder
    )
    return(mean((ghat_boot_mean - Y) ^ 2))
  } # end-function-objetive


  estimation_mean_and_bws <-
    function(i,
             Ydat,
             Xdat,
             nboot,
             bootstrap,
             ckerorder) {
      Y <- Ydat
      X <- Xdat
      ghat.bw <- np::npregbw(
        xdat = X[, i],
        ydat = Y[, i],
        ckertype = "epanechnikov",
        bwmethod = "cv.ls"
      )
      ghat <- np::npreg(ghat.bw, residuals = TRUE)
      if (!bootstrap) {
        return(list(
          mean = ghat$mean,
          bws = ghat$bw,
          x = X[, i],
          y = Y[, i]
        ))
      } #end-if-!bootstrap
      else if (bootstrap) {
        resid.np <- np::npreg(np::npregbw(
          xdat = X[, i],
          ydat = residuals(ghat),
          ckertype = "epanechnikov"
        ))
        g <- ghat$mean
        sderr <- resid.np$merr
        sderr.idx <- sderr <= 1e-5
        error  <-
          (residuals(ghat) - mean(residuals(ghat))) / sderr
        error[sderr.idx] <- 0

        n <- length(g)
        Y_boot <- mapply(
          FUN = function(n)
          {
            idx <- sample(1:n, replace = TRUE)
            return(g +  sderr * error[idx])
          },
          # end-function-Y_boot
          rep(n, nboot),
          SIMPLIFY = FALSE
        )

        cv.opt <- optimize(
          f = objective,
          interval = c(0,  2 * ghat$bw),
          Y_boot = Y_boot,
          Y = Y[, i],
          Xi = X[, i],
          ckerorder = ckerorder
        )

        h_boot_min <- cv.opt$minimum

        g_boot_min_value <-
          ghat_boot(
            Y_boot = Y_boot,
            Xi = X[, i],
            h = h_boot_min,
            ckerorder = ckerorder
          )

        return(list(
          mean = g_boot_min_value,
          bws = h_boot_min,
          x = X[, i],
          y = Y[, i]
        ))

      } #end-if-bootstrap
    } #end-function-estimation_mean_and_bws

  nc <-  min(ncol(Y), parallel::detectCores())

  ANS <- try(parallel::mclapply(
    X = 1:ncol(Y),
    FUN = estimation_mean_and_bws,
    Ydat = Y,
    Xdat = X,
    nboot = nboot,
    bootstrap = bootstrap,
    ckerorder = ckerorder,
    mc.cores = nc
  ),
  silent = T)

  if (is(ANS, 'try-error')) {
    #Windows does not support mclapply...
    NP <- lapply(
      X = 1:ncol(Y),
      FUN = estimation_mean_and_bws,
      Ydat = Y,
      Xdat = X,
      nboot = nboot,
      bootstrap = bootstrap,
      ckerorder = ckerorder
    )
  } else {
    NP <<- ANS
  }

  return(NP)
} #end-optNP

normalize <- function (X,
                       MAXS = NULL,
                       MINS = NULL,
                       inv = F) {
  #DOC
  #Normalizes a vector or a matrix according to the given 'MAXS' and 'MINS'.
  #If 'MAXS' and 'MINS' are not provided 'X' is scaled so that each column is within [0,1].
  #ARGUMENTS
  #X: matrix or vector to normalize
  #MAXS, MINS: maxima and minima. NULL or NA values are set to max(X[,i]) and min(X[,i]) respectively.
  #inv: if T performs the inverse transformation
  #MAIN
  X = as.matrix(X)
  if (inv) {
    return(scale(
      scale(X, center = F, scale = (MAXS - MINS) ** (-1)),
      center = -MINS,
      scale = F
    ))
  } else {
    if (is.null(MAXS))
      MAXS = apply(X, 2, max)
    if (is.null(MINS))
      MINS = apply(X, 2, min)
    for (i in 1:ncol(X)) {
      if (is.na(MAXS[i]))
        MAXS[i] = max(X[, i])
      if (is.na(MINS[i]))
        MINS[i] = min(X[, i])
    }
    return(scale(X, center = MINS, scale = MAXS - MINS))
  }
}

plot.sobolnp <- function(x, ...) {
  yrng = range(c(1 + x[['S']][, 'se'], x[['S']][, 'q0.05']))
  #plot estimates
  plot(
    x = 1:nrow(x[['S']]),
    y = x[['S']][, 'Si'],
    pch = 19,
    ylim = yrng,
    xaxt = 'n',
    xlab = 'parameter',
    ylab = 'Si and se',
    ...
  )
  #plot q 0.05
  points(
    x = 1:nrow(x[['S']]),
    y = x[['S']][, 'q0.05'],
    pch = 19,
    col = 2
  )
  #plot se
  arrows(
    x0 = 1:nrow(x[['S']]),
    y0 = x[['S']][, 'Si'],
    y1 = x[['S']][, 'Si'] - x[['S']][, 'se'],
    length = 0
  )	#lower
  arrows(
    x0 = 1:nrow(x[['S']]),
    y0 = x[['S']][, 'Si'],
    y1 = x[['S']][, 'Si'] + x[['S']][, 'se'],
    length = 0
  )	#upper
  #plot 0
  abline(h = 0, lty = 2)
  #x axis
  axis(side = 1,
       at = 1:nrow(x[['S']]),
       labels = row.names(x[['S']]))
  #legend
  legend(
    'topright',
    legend = c('Si', 'se', 'q0.05'),
    pch = c(19, NA, 19),
    lty = c(NA, 1, NA),
    col = c(1, 1, 2),
    horiz = T,
    bty = 'n'
  )
}

print.sobolnp <- function(x, ...) {
  cat("\nCall:\n", deparse(x[['call']]), "\n", sep = "")
  cat("\nNumber of observations:", length(x[['Y']]), "\n")
  if (x[['bootstrap']]) {
    cat("\nBootstrap Samples:", x[['nboot']], "\n")
  }
  cat("\nFirst order indices:\n")
  print(x[['S']])
}
