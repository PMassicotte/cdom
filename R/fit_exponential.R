#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FILE:         fit_exponential.R
#
# AUTHOR:       Philippe Massicotte
#
# DESCRIPTION:  Fit an exponential curve to CDOM data.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#' Fit an exponential model to CDOM data.
#'
#' @details \deqn{y = a0 + e^{(-S(x - \lambda_0))} + K}
#'
#'
#' @param wl The wavelength vector.
#' @param spectra The spectra vector.
#' @param wl0 The reference wavelength (ex.: 350).
#' @param startwl The starting wavelength (ex.: 240).
#' @param endwl The ending wavelength (ex.: 600).
#'
#' @return A list contianing the \code{nls} object and the \code{R2}. NULL if
#' the model did not converged.
#' @export
#'
#' @import minpack.lm
#'
#' @examples
#' # Fit an exponential model using the reference wavelength 350 between 190 and 900 nm.
#'
#' data(spectra)
#'
#' fit <- fit_exponential(spectra$wavelength, spectra$absorbance, 350, 190, 900)
#' summary(fit)
#'
#' plot(spectra$wavelength, spectra$absorbance)
#' lines(spectra$wavelength, predict(fit), col = "red")

fit_exponential <- function(wl, spectra, wl0, startwl, endwl){

  if(length(wl) != length(spectra)){
    stop("wl and spectra are not of the same length.")
  }

  if(!is.numeric(wl) | !is.numeric(spectra)){
    stop("wl and spectra need to be numeric.")
  }

  #--------------------------------------------
  # Get a0 value.
  #--------------------------------------------
  sf <- splinefun(wl, spectra)
  a0 <- sf(wl0)

  #--------------------------------------------
  # Extract CDOM data based on user inputs.
  #--------------------------------------------
  xx <- wl[which(wl >= startwl & wl <= endwl)]
  yy <- spectra[which(wl >= startwl & wl <= endwl)]

  #--------------------------------------------
  # Fit the data.
  #--------------------------------------------
  control <- list(minFactor = 1e-10,
                  warnOnly = FALSE,
                  maxiter = 1024,
                  maxfev = 600)

  out <- tryCatch(
    {
      fit <- nlsLM(yy ~ a0 * exp(-S*(xx - wl0)) + K,
                   start = c(S = 0.02, K = 0.01, a0 = a0),
                   lower = c(S = 0, K = -2, a0 = 0),
                   upper = c(S = 1, K = 2, a0 = max(yy)),
                   control = control)

      fit$R2 <- 1 - sum((yy - predict(fit))^2) / (length(yy) * var(yy))

      return(fit)

    },error = function(cond) {
      message("Error in fit_exponential() when trying to fit. Check your data.")
      message(cond)
      # Choose a return value in case of error
      return(NULL)

    },warning = function(cond) {
      message(cond)
      message("Warning in fit_exponential().")
      # Choose a return value in case of warning
      return(NULL)

    },finally={

    }
  )








}
