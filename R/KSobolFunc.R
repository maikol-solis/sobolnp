# library(np)
# library(matlab)
# library(parallel)
# library(boot)
# library(pbmcapply)

sobolNP <-
  function(data,
           bandwidth,
           ncores = 2,
           bandwidth.compute = TRUE,
           bootstrap,
           ckerorder = 2) {

    n <- dim(data)[1]
    p <- dim(data)[2] - 1
    
    #if no bandwith is shown, the bandwith is found by bootstrapping
    if (missing(bandwidth)) {
      if (bootstrap == TRUE) {
        #browser()
        bandwidth <-
          bootsobol(data,
                    nboot = 100,
                    ncores = ncores,
                    ckerorder = ckerorder)
        return(bandwidth)
      } else if (bootstrap == FALSE) {
        bandwidth <- rep(NULL, p)
      } else {
        message("The parameter 'bootstrap' has no default. Possible values: TRUE or FALSE. ")
      }
      
    } else {
      for (i in 1:p) {
        if (bandwidth[i] <= 0) {
          return(list(idx = rep(NaN, p),
                      bandwidth = rep(NaN, p)))
        }
      }
    }
    
    
    
    X <- data[, 2:(p + 1)]
    Y <- data[, 1]
    
    g <- matrix(nrow = n, ncol = p) #matrix of X
    resid <- matrix(nrow = n, ncol = p) #matrix of X
    error <- matrix(nrow = n, ncol = p) #matrix of X
    sderr.resid <- matrix(nrow = n, ncol = p) #matrix of X
    bws <- matrix(nrow = 1, ncol = p) #matrix of Y
    
    for (k in 1:p) {
      ghat.bw <- npregbw(
        xdat = X[, k],
        ydat = Y,
        bws = bandwidth[k],
        bandwidth.compute = bandwidth.compute,
        ckertype = "epanechnikov",
        ckerorder = ckerorder,
        bwmethod = "cv.ls"
      )
      # ghat.bias <-
      #   npplot(
      #     ghat.bw,
      #     ydat = Y,
      #     xdat = X[, k],
      #     neval = n,
      #     plot.errors.method = "bootstrap",
      #     plot.errors.center = "bias-corrected",
      #     plot.behavior = "data",
      #     plot.errors.boot.num = 100
      #   )
     
      # boot.out <-
      #   boot(
      #     data = data.frame(xdat = X[, k], ydat = Y),
      #     statistic = function(d, indices) {
      #       npreg(txdat = d$xdat[indices],
      #             tydat = d$ydat[indices],
      #             exdat = d$xdat[indices],
      #             ckertype = "epanechnikov",
      #             bws =  ghat.bw)$mean
      #     },
      #     R = 200
      #   )
      #browser()
      bws[, k] <- ghat.bw$bw
      ghat <- npreg(ghat.bw, residuals = TRUE)
      # ghat <-
      #   npreg(ghat.bw, exdat = ghat.bias$r1$eval[, 1], residuals = TRUE)
 
      resid[, k] <- residuals(ghat)
      
      resid.np <- npreg(npregbw(
        xdat = X[, k],
        ydat = resid[,k],
        ckertype = "epanechnikov"
      ) )
      
      g[, k] <- ghat$mean
      sderr.resid[, k] <- resid.np$merr
      error[, k]  <- (resid[, k] - mean(resid[, k])) / sderr.resid[, k]
        
        
    }
    #browser()
    sigma.Y <- var(Y)
    var.g <- apply(g, 2, var)  #apply(g, 2, var)
    sobol.idx <- var.g / sigma.Y
    #message("Sobol final: ", sobol.idx)
    return(
      list(
        idx = sobol.idx,
        bandwidth = bws,
        resid = resid,
        error = error,
        mean = g, 
        sderr = sderr.resid
      )
    )
  }

cv_sobol <-
  function(h,
           Sobol_original,
           Y_boot,
           X,
           nboot,
           nindex,
           ckerorder,
           minimizer = TRUE,
           ...)
  {
    if (h <= 0) {
      return(NaN)
    }
    sigma <- Sobol_original$idx[nindex]
    sigma.boot <- mclapply(
      X = Y_boot,
      FUN = function(ydat, xdat, h = 0.5, nindex, ...) {
        # Ksobol(
        #   data = data,
        #   bandwidth = rep(const, times = dim(data_boot[[1]])[2] - 1),
        #   ckerorder = ckerorder
        # )$idx[nindex]
        #browser()
        # print(xdat)
        # message(nindex)
        ghat.bw <- npregbw(
          xdat = xdat[, nindex],
          ydat = ydat[, nindex],
          bws = h,
          bandwidth.compute = FALSE,
          ckertype = "epanechnikov",
          ckerorder = ckerorder,
          bwmethod = "cv.ls"
        )
        ghat <- npreg(ghat.bw)
        g <- ghat$mean 
        # sigma.Y <- var(ydat[, nindex])
        # var.g <- mean(g ^ 2) - mean(ydat[, nindex]) ^ 2 #var(g)
        # sobol.idx <- var.g / sigma.Y
        # return(sobol.idx)
        return(list(
          meang = g,
          varY = var(ydat[, nindex]),
          meanY = mean(ydat[, nindex])
        ))
      },
      xdat = X,
      h = h,
      nindex = nindex,
      SIMPLIFY = TRUE,
      mc.cores = 20
    )
    #message("Sobol indice #", nindex, " = ",  mean(simplify(sigma.boot)))
    Sboot.list <- simplify2array(sigma.boot)
    # browser()
    gboot.corrected.bias <-
      apply(simplify2array(Sboot.list["meang", ]), 1, mean)
    
    g.corrected <-
      2 * Sobol_original$mean[, nindex] - gboot.corrected.bias
    
    Sobol.boot <- var(g.corrected) / simplify(Sboot.list["varY",])
    
    # if (mean(Sobol.boot) < 0)
    #   return(NA)
    
    cv_fun <- mean((Sobol.boot - sigma) ^ 2)
    # message("CV: ", cv_fun, " bw: ", h)
    # message("Sobol base: ", sigma)
    # message("Sobol estimado: ", mean(simplify(sigma.boot)))
    if (minimizer == TRUE) {
      if (is.finite(cv_fun)) {
        return(cv_fun)
      } else{
        return(NA)
      }
    } else {
      return(Sobol.boot)
    }
    
  }

bootsobol <- function (data,
                       nboot = 100,
                       ncores = ncores,
                       ckerorder) {
  n <- dim(data)[1]
  m <- dim(data)[2] - 1
  
  # data_boot <- mcmapply(
  #   FUN = function(n) {
  #     idx <- sample(1:n, replace = TRUE)
  #     return(data[idx,])
  #   },
  #   rep(n, nboot),
  #   SIMPLIFY = FALSE,
  #   mc.cores = ncores
  # )
  
  #estimate original value of bandwidth
  # bw <- NULL
  # for(i in 1:m){
  #   bw <- rbind(bw,cv_bws_npreg(y=data[,1],x=data[,m])$best.bw)
  # }
  #
  Sobol_original <-
    sobolNP(
      data = data,
      # bandwidth = bw,
      ncores = ncores,
      bandwidth.compute = TRUE,
      bootstrap = FALSE,
      ckerorder = ckerorder
    )
  #browser()
  Y_boot <- pbmcmapply(
    FUN = function(n) {
      #browser()
      idx <- sample(1:n, replace = TRUE)
      # epsilon <- .75 * (1 - runif(n = 1, min = -1, max = 1) ^ 2)
      # return(data[idx, ] + epsilon * Sobol_original$bandwidth)
      p <- length(Sobol_original$idx)
      return(Sobol_original$mean +  Sobol_original$sderr * Sobol_original$error[idx,])
    },
    rep(n, nboot),
    SIMPLIFY = FALSE,
    mc.cores = ncores
  )
  
  cv.min <- NULL
  cv.value <- NULL
  for (i in 1:m) {
    # cv.opt <- optimx(
    #   par = Sobol_original$bandwidth[i],
    #   fn = cv_sobol,
    #   Sobol_original = Sobol_original,
    #   Y_boot = Y_boot,
    #   X = data[, -1],
    #   nboot = nboot,
    #   nindex = i,
    #   ckerorder = ckerorder,
    #   hessian = TRUE
    #   # control = list(all.methods = TRUE)
    #   # method = "Brent",
    #   # lower = 0,
    #   # upper = 10
    # )
    #browser()
    cv.opt <- nlm(
      f = cv_sobol,
      p = Sobol_original$bandwidth[i],
      Sobol_original = Sobol_original,
      Y_boot = Y_boot,
      X = data[, -1],
      nboot = nboot,
      nindex = i,
      ckerorder = ckerorder,
      minimizer = TRUE,
      print.level = 0
    )
    cv.min <- cbind(cv.min, cv.opt$estimate)
    # cv.opt <-
    #   GenSA(
    #     par = Sobol_original$bandwidth[i],
    #     fn = cv_sobol,
    #     lower = 0,
    #     upper = Inf,
    #     Sobol_original = Sobol_original,
    #     data_boot = data_boot,
    #     nboot = nboot,
    #     nindex = i,
    #     control = list(verbose = TRUE, maxit = 3)
    #   )
    
    
    cv.min.sobol <- cv_sobol(
      cv.opt$estimate,
      Sobol_original = Sobol_original,
      Y_boot = Y_boot,
      X = data[, -1],
      nindex = i,
      nboot = nboot,
      ckerorder = ckerorder,
      minimizer = FALSE
    )
    cv.value <- cbind(cv.value,  mean(cv.min.sobol))
    # browser()
    # cv.boot.mean  <-
    #   BCa(x = cv.min.sobol,
    #       delta = 0.01,
    #       theta = mean)
    # cv.value <-
    #   cbind(cv.value,  cv.boot.mean[3] )
    # cv.boot.mean <- boot(
    #   data = cv.min.sobol,
    #   statistic = function(data, indices) {
    #     return(mean(data[indices]))
    #   },
    #   R = 200,
    #   type = "bca",
    #   parallel = "multicore",
    #   ncpus = ncores
    # )
    # message("Sobol biased: ", mean(cv.min.sobol))
    # message("Bias: ", mean(cv.boot.mean$t) - mean(cv.min.sobol))
    # message("Sobol corrected-bias: ",
    #         2 * mean(cv.min.sobol) - mean(cv.boot.mean$t))
    # cv.value <-
    #   cbind(cv.value, 2*mean(cv.min.sobol) - mean(cv.boot.mean$t) )
  }
  
  # message("bw: ", paste(cv.min, collapse = "\n"))
  # message("Sobol base: ", paste(Sobol_original$idx, collapse = "\n"))
  # message("Sobol estimado: ", paste(cv.value, collapse = "\n"))
  
  return(list(idx = cv.value, bandwidth = cv.min))
}

# cv_bws_npreg <-
#   function(x,
#            y,
#            bandwidths = seq(1, 2, length.out = 100),
#            training.fraction = 0.9,
#            folds = 10) {
#     require(np)
#     n = length(x)
#     # Sanity-check inputs
#     stopifnot(n > 1, length(y) == n)
#     stopifnot(length(bandwidths) > 1) # CV pointless otherwise
#     stopifnot(training.fraction > 0, training.fraction < 1)
#     stopifnot(folds > 0, folds == trunc(folds))
#
#     # Make a matrix to store MSEs for each fold/bandwidth combination
#     fold_MSEs = matrix(0, nrow = folds, ncol = length(bandwidths))
#     # Name the columns after bandwidths for easy reference later
#     # coerces the numerical bandwidths into character strings
#     colnames(fold_MSEs) = bandwidths
#
#     n.train = max(1, floor(training.fraction * n))
#     for (i in 1:folds) {
#       train.rows = sample(1:n, size = n.train, replace = FALSE)
#       x.train = x[train.rows]
#       y.train = y[train.rows]
#       x.test = x[-train.rows]
#       y.test = y[-train.rows]
#       for (bw in bandwidths) {
#         fold_MSEs[i, paste(bw)] <-
#           npreg(
#             txdat = x.train,
#             tydat = y.train,
#             exdat = x.test,
#             eydat = y.test,
#             bws = bw
#           )$MSE
#         # see help(npreg) for arguments to that function
#         # paste(bw): turns numerical bandwidth to type character, so R knows it's
#         # the name of a column, not a column index
#       }
#     }
#     CV_MSEs = colMeans(fold_MSEs)
#     best.bw = bandwidths[which.min(CV_MSEs)]
#     return(list(
#       best.bw = best.bw,
#       CV_MSEs = CV_MSEs,
#       fold_MSEs = fold_MSEs
#     ))
#   }

# plotErrorFunction <-
#   function(lBound, uBound, steps, data, nboot, nindex = 1) {
#     n <- dim(data)[1]
#     m <- dim(data)[2] - 1
#     c.vec <- seq(from = lBound,
#                  to = uBound,
#                  length.out = steps)
#
#     Sobol_original <- Ksobol(
#       data = data,
#       bandwidth = rep(0.5, dim(data)[2] - 1),
#       ncores = 1,
#       bandwidth.compute = TRUE
#     )
#     data_boot <- mcmapply(
#       FUN = function(n) {
#         idx <- sample(1:n, replace = TRUE)
#         return(data[idx,])
#       },
#       rep(n, nboot),
#       SIMPLIFY = FALSE,
#       mc.cores = 3
#     )
#
#     cv.vec <- NULL
#     for (k in 1:length(c.vec)) {
#       print(k)
#       s <- cv_sobol(
#         c.vec[k],
#         Sobol_original = Sobol_original,
#         data_boot = data_boot,
#         nboot = nboot,
#         nindex = nindex
#       )
#       cv.vec <- rbind(cv.vec, s)
#     }
#     cv.num <- as.numeric(cv.vec)
#
#     plot.new()
#     plot(c.vec, cv.num, type = "l")
#   }
