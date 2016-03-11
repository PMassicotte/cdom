#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## Gaussian probability density function
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
gaussian <- function(x, pos, wid){
  g <- exp(-((x - pos) / (0.6005615 * wid))^2)

  return(g)
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
##
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
fitgaussian <- function(lambda){

  x <- get("x", envir = .GlobalEnv)
  y <- get("y", envir = .GlobalEnv)

  numpeaks <- round(length(lambda) / 2)

  A <- optimbase::zeros(length(x), numpeaks)

  for(j in 1:numpeaks) {

    A[, j] <- gaussian(x, lambda[2 * j - 1], lambda[2 * j])
  }

  A <- cbind(pracma::ones(n = length(y), 1), A)

  PEAKHEIGHTS <- abs(pracma::mldivide(A,y))

  z <- A %*% PEAKHEIGHTS

  err <- pracma::Norm(z - y)

  assign("PEAKHEIGHTS", PEAKHEIGHTS, envir = .GlobalEnv)

  return(err)
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## Main method to find the desired numbers of Gaussian components
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
findgaussian <- function(NumPeaks, x, y, NumTrials){

  n <- max(x) - min(x)
  MINWIDTH <- x[2] - x[1]
  startpos <- seq(from = n / (NumPeaks + 1),
                  to = n - (n / (NumPeaks + 1)),
                  by = n / (NumPeaks + 1)) + min(x)

  start <- NULL
  for(marker in 1:NumPeaks) {

    markx <- startpos[marker]

    start <- c(start, markx, n/(3*NumPeaks))

  }

  myoptions <- neldermead::optimset(MaxIter = 1000,
                        MaxFunEvals = 1000,
                        TolX = 0.001,
                        Display = FALSE)
  newstart <- start
  LowestError <- 10000000

  for(k in 1:NumTrials) {

    ## Set the optimisation options and proceed
    opt <- neldermead::fminsearch(x0 = newstart,
                                  fun = fitgaussian,
                                  options = myoptions)

    bc <- buildcomponents(opt = opt,
                          x = x,
                          NumPeaks = NumPeaks,
                          MINWIDTH = MINWIDTH)

    model <- rowSums(bc$components)

    MeanFitError <- 100 * pracma::Norm(y - bc$model) / (sqrt(length(y)) * max(y))

    ## We found better parameter estimations
    if(MeanFitError < LowestError){
      LowestError <- MeanFitError
      BestStart <- newstart
      BestComponents <- bc$components
      BestParams <- bc$params
    }

    ## Randomize parameters to setup a new trial
    delta <- 1

    for(parameter in seq(from = 1, by = 2, to = 2*NumPeaks)) {

      newstart[parameter] <- newstart[parameter] * (1 + delta * (pracma::rand(1, 1) - 0.5) / 500)

      newstart[parameter + 1] <- newstart[parameter + 1] * (1 + delta * (pracma::rand(1, 1) - 0.5) / 100)

    }
  }

  ## Here we are done with NumTrials
  myBIC <- calculateBIC(rowSums(BestComponents),
                        y,
                        ncol(BestComponents) * 3 + 1)

  res <- list(params = BestParams,
              components = BestComponents,
              BIC = myBIC,
              MeanFitError = MeanFitError)

  return(res)
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## Method to build the y vectors from the Gaussian estimated parameters.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
buildcomponents <- function(opt, x, NumPeaks, MINWIDTH){

  components <- optimbase::zeros(length(x), NumPeaks)
  A <- optimbase::zeros(length(x), NumPeaks)
  params <- list()

  for(i in 1:NumPeaks){

    if(opt$optbase$xopt[2*i] < MINWIDTH){
      opt$optbase$xopt[2*i] <- abs(opt$optbase$xopt[2*i])
    }


    p0est <- PEAKHEIGHTS[i+1]
    p1est <- opt$optbase$xopt[2*i-1]
    p2est <- opt$optbase$xopt[2*i]/(2*sqrt(2*log(2)))

    components[,i] <- p0est*exp(-((x-p1est)^2/(2*p2est^2)))
    A[, i] <- gaussian(x, opt$optbase$xopt[2*i-1], opt$optbase$xopt[2*i])

    params$p0[i] <- PEAKHEIGHTS[i+1]
    params$p1[i] <- opt$optbase$xopt[2*i-1]
    params$p2[i] <- opt$optbase$xopt[2*i]/(2*sqrt(2*log(2))) # This is the "gaussian" width

  }

  baseline <- PEAKHEIGHTS[1]
  Heights <- PEAKHEIGHTS[2:(1+NumPeaks)]
  model <- (Heights %*% t(A) +baseline)
  res <- list(components = components, params = params, model = model)
  return(res)
}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## Method to the loglik value.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
loglik <- function(yhat, y){

  N <- length(y)
  w <- rep(1, N)
  zw <- rep(FALSE, N)

  myresiduals <- y - yhat

  val <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + log(sum(w * myresiduals^2)))/2

  return(val)

}

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## Methods to calculate the BIC metric.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
calculateBIC <- function(yhat, y, k){

  ll <- loglik(yhat = yhat, y = y)

  BIC <- -2 * ll + k * log(length(y))

  return(BIC)

}

getBIC <- function(x){
  return(x$BIC)
}

getmeanfiterror <- function(x){
  return(x$MeanFitError)
}


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
## This method loops from nmin to nmax model components and return the best one
## based on the BIC values.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
getbestmodel <- function(x, y, nmin = 1, nmax = 5, NumTrials = 1, stopmethod){

  n <- nmax - nmin + 1

  res <- vector("list", n)

  for(i in nmin:nmax){
    res[[i]] <- findgaussian(NumPeaks = i, x = x, y = y, NumTrials = NumTrials)
  }

  if(stopmethod == "BIC"){
    ## Extract the BIC for all models
    BICs <- unlist(lapply(res, getBIC))

    ## Find the one with the lowest BIC and return it
    index <- which.min(BICs)
  }else if(stopmethod == "meanfiterror") {
    ## Extract the mean fit error for all models
    meanfiterror <- unlist(lapply(res, getmeanfiterror))

    ## Find the one with the lowest BIC and return it
    index <- which.min(meanfiterror)
  }


  return(res[[index]])

}



#================================================================================
## Get the Gauss component based on components found un a nls model.
#================================================================================
.getGaussComponent <- function(x, fit){

  ncomp <- length(grep("p",names(coef(fit))))/3

  res <- data.frame(matrix(nrow = length(x), ncol = ncomp))
  colnames(res) <- paste("gauss", 1:ncomp, sep = "")

  pos <- 4
  for(i in 1:ncomp){

    p0 <- coef(fit)[pos]
    p1 <- coef(fit)[pos + 1]
    p2 <- coef(fit)[pos + 2]

    pos <- pos + 3

    res[, i] <- p0*exp(-0.5*((x-p1)/p2)^2)

  }

  return(res)
}

#================================================================================
#================================================================================
.sp <- function(x, Y){

  ncomp <- ncol(Y)

  cols <- brewer.pal(name = "Dark2", 8)

  m <- lapply(Y, function(x) {which(x > 1e-10)})
  mymin <- x[min(unlist(m))]
  mymax <- x[max(unlist(m))]

  plot(x, rowSums(Y), type = "l", axes = FALSE, xlab = "Wavelength (nm)", ylab = "abs.", xlim = c(mymin, mymax), col = "gray", lwd = 3)

  axis(1, tck = -0.02)
  axis(2, tck = -0.02, las = 1)
  box(bty = "l")

  for(i in 1:ncomp){
    lines(x, Y[,i], type = "l", col = cols[i], lty = 2)
  }
}

#================================================================================
## Plot all xxx.
#================================================================================
.plotGlobalFit <- function(fitCDOMobj){

  if(is.null(fitCDOMobj)){
    stop("Can not plot fitCDOMobj")
    return
  }


  plot(fitCDOMobj$xraw, fitCDOMobj$yraw, col = "gray75", axes = FALSE, xlab = "wavelength (nm)", ylab = expression(paste("Absorption (", m^{-1}, ")")))
  axis(1, tck = 0.02)
  axis(2, tck = 0.02, las = 2)
  box()

  lines(fitCDOMobj$xraw, predict(fitCDOMobj$fit), col = "red")

  components <- .getGaussComponent(x = fitCDOMobj$xraw, fit = fitCDOMobj$fit)

  xpos <- par()$usr[2]*0.85
  ypos <- par()$usr[4]*0.75
  Hmisc::subplot(.sp(fitCDOMobj$xraw, components), xpos, ypos, size = c(2,1), pars = list(cex = 0.75))

}
