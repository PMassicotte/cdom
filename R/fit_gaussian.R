.removeleverageexp <- function(x,y){

  xx <- x
  yy <- y

  fit <- cdom_fit_exponential(wl = xx,
                              absorbance = yy,
                              wl0 = 350,
                              startwl = min(xx),
                              endwl = max(xx))

  #   fit <- tryCatch(nlsLM(formula = yy ~ a0 * exp(-S*(xx - 350)) + K, control = list(maxiter = 1024, maxfev = 2000), start = list(S = 0.05, K = 0, a0 = 0.1))
  if(is.null(fit)){return(NULL)}

  flag <- TRUE

  while(flag == TRUE){

    #plot(xx, yy)
    #lines(x, predict(fit, newdata = list(xx = x)), col = "red")

    toremove <- which(fit$data$.resid > 1 * mean(abs(fit$data$.resid)))

    ## Do not remove points at the very low wavelength
    toremove <- toremove[which(toremove - seq(1:length(toremove)) > 0)]


    #points(xx[toremove], yy[toremove], pch = 21, col = "red")

    if(length(toremove) > 0) {

      xx <- xx[-toremove]
      yy <- yy[-toremove]

      fit <- cdom_fit_exponential(wl = xx,
                                  absorbance = yy,
                                  wl0 = 350,
                                  startwl = min(xx),
                                  endwl = max(xx))

      if(is.null(fit)) {

        return(NULL)

      }

      ## Do not attempt to model data if not "enough" observations.
      if(length(xx) < 50) {

        return(list(x = xx, y = yy, model = fit$model))

      }


    } else {

      flag = FALSE

    }

  }

  return(list(x = xx, y = yy, model = fit$model))
}


func1 <- function(par, x, y, functooptimise){

  fo <- formula(sub(".*~", "~", deparse(functooptimise)))
  func <- gsubfn::fn$identity(fo)
  y. <- do.call(func, c(list(x = x), par))

  sum((y - y.)^2)

}


fitCDOMcomponents <- function(x,y, stopmethod = c("BIC", "meanfiterror")){

  expfit <- .removeleverageexp(x = x, y = y)

  if(is.null(expfit)){return(NULL)}

  ff <- predict(expfit$model, newdata = list(x = x))

  res <- y - ff

  res <- signal::sgolayfilt(res, p = 3, n = 21)

  xx <- x[which(x < 500)]
  res <- res[which(x < 500)]

  ## Temporary hack...
  assign("x", xx, envir = .GlobalEnv)
  assign("y", res, envir = .GlobalEnv)

  ## Get initial estimations for peaks only. This script is based on a matlab toolbox.
  bm <- getbestmodel(x = xx, y = res, nmin = 1, nmax = 5, NumTrials = 1, stopmethod = stopmethod)

  ## Sometimes no peaks are found
  if(is.null(bm)){return(NULL)}

  ## Create a data frame with peak estimations.
  peaks <- data.frame(peak = 1:ncol(bm$components), p1 = unlist(bm$params$p1), p0 = unlist(bm$params$p0), p2 = unlist(bm$params$p2))

  ## Construct a fit including the exponential part + gauss components
  C <- "(p0%s*exp(-((x-p1%s)^2/(2*p2%s^2))))"
  n <- nrow(peaks)

  ## Create the equation with the appropriate number of Gaussian components
  C <- rep(C, n)
  C <- paste(C, collapse = "+")

  C <- do.call(sprintf, as.list(c(C, rep(letters[1:n], each = 3, times = 1))))

  exppart <- "y ~ a0 * exp(-S * (x - 350)) + K +"
  f <- paste(exppart, C)
  f <- as.formula(f)

  ## Set the starting values for the parameters
  mystart  <- c(coefficients(expfit$model)[c(3,1,2)], as.vector(t(as.matrix(peaks[, c(3,2,4)]))))
  names(mystart)  <- all.vars(f)[-c(1,4)]

  bestguesses <- optim(par = mystart,
                       fn = func1 ,
                       x = x,
                       y = y,
                       functooptimise = f,
                       control = list(maxit = 5000),
                       lower = rep(0, length(mystart)),
                       method = "L-BFGS-B")

  print(bestguesses$par)

  # bestguesses2 <- nmkb(par = mystart,
  #                      fn = func1 ,
  #                      x = x,
  #                      y = y,
  #                      functooptimise = f,
  #                      control = list(maxfeval = 5000),
  #                      lower = rep(0, length(mystart)))
  #
  # names(bestguesses2$par) = names(mystart)

  fit <- tryCatch(nlsLM(formula = f,
                        start = bestguesses$par,
                        data = data.frame(y = y, x = x),
                        control = list(maxiter = 1024, maxfev = 60000),
                        upper = bestguesses$par,
                        lower = bestguesses$par),
                  error=function(e) NULL,  warning=function(e) NULL)

  if(is.null(fit)){

    return(NULL)

  }else{

    # write.csv(data.frame(mystart = mystart, estimated = coefficients(fit)), "tmp/diagnosticNLS.csv")
    #
    # write.csv(x = cbind(sprintf("%2.2f", bestguesses$par), sprintf("%2.2f", coef(fit))), file = "tmp/compare_optim_nls.csv", row.names = FALSE)
    #
    # write.csv(data.frame(x = xx, y = res), "tmp/expfit.csv")

    res <- list(fit = fit, xraw = x, yraw = y)

    return(res)
  }

}
