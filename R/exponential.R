#' Fit an exponential model to CDOM data.
#'
#' @details \deqn{y = a0 * e^{(-S(x - \lambda_0))} + K}
#'
#' @param wl The wavelength vector.
#' @param absorbance The absorbance vector.
#' @param wl0 The reference wavelength (ex.: 350).
#' @param startwl The starting wavelength (ex.: 240).
#' @param endwl The ending wavelength (ex.: 600).
#'
#' @return A list containing:
#' \describe{
#'   \item{params}{A data frame with values of fitted parameters.}
#'   \item{r2}{R2 of the nls model.}
#'   \item{data}{A data frame with fitted (predicted) values of the model.}
#' }
#'
#' The function will return \code{NULL} if the model did not converged.
#' @export
#' @import minpack.lm
#' @importFrom broom tidy augment
#' @importFrom stats splinefun predict var coef lm na.omit
#' @importFrom purrr safely
#'
#' @examples
#' # Fit an exponential model using the reference wavelength 350 between 190 and 900 nm.
#'
#' data(spectra)
#'
#' fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#'
#' plot(spectra$wavelength, spectra$spc1)
#' lines(spectra$wavelength, fit$data$.fitted, col = "red")

cdom_exponential <- function(wl, absorbance, wl0 = 350, startwl, endwl){

  stopifnot(length(wl) == length(absorbance),
            is.numeric(absorbance),
            is.numeric(wl),
            is.vector(wl),
            is.vector(absorbance),
            is.numeric(wl0),
            is.numeric(startwl),
            is.numeric(endwl),
            wl0 > min(wl) & wl0 < max(wl))

  if (missing(startwl)) {startwl = min(wl)}
  if (missing(endwl)) {endwl = max(wl)}

  #--------------------------------------------
  # Get a0 value.
  #--------------------------------------------
  sf <- splinefun(wl, absorbance)
  a0 <- sf(wl0)

  #--------------------------------------------
  # Extract CDOM data based on user inputs.
  #--------------------------------------------
  x <- wl[which(wl >= startwl & wl <= endwl)]
  y <- absorbance[which(wl >= startwl & wl <= endwl)]

  df <- data.frame(x = x, y = y)

  #--------------------------------------------
  # Fit the data.
  #--------------------------------------------
  control <- list(minFactor = 1e-10,
                  warnOnly = FALSE,
                  maxiter = 1024,
                  maxfev = 600)

  safe_nls <- purrr::safely(minpack.lm::nlsLM)

  fit <- safe_nls(y ~ a0 * exp(-S * (x - wl0)) + K,
                 start = c(S = 0.02, K = 0.01, a0 = a0),
                 lower = c(S = 0, K = -Inf, a0 = 0),
                 upper = c(S = 1, K = Inf, a0 = max(y)),
                 control = control,
                 data = df)

  if (is.null(fit$error)) {

    r2 <- 1 - sum((y - predict(fit$result))^2) / (length(y) * var(y))

    res <- list(model = fit$result, r2 = r2, x = x, y = y)
    class(res) <- "exponential_fit"

    return(res)

  } else {
    return(NULL)
  }
}

#' Predict method for CDOM exponential fit
#'
#' @param object An object returned by \code{cdom_exponential}.
#' @param ... other arguments.
#'
#' @return A numerical vector with predicted values.
#' @export
#'
#' @examples
#' data(spectra)
#'
#' fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#' predict(fit)
predict.exponential_fit <- function(object, ...) {

  res <- predict(object$model)

  return(res)

}

#' Extract Model Coefficients from a CDOM exponential fit
#'
#' @param object An object returned by \code{cdom_exponential}.
#' @param ... other arguments.
#'
#' @return A numerical vector with estimated coefficients.
#' @export
#'
#' @examples
#' data(spectra)
#'
#' fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#' coef(fit)
coef.exponential_fit <- function(object, ...) {

  res <- coef(object$model)

  return(res)

}

#' Plot a Fitted CDOM Exponential Curve
#'
#' @param x An object returned by \code{cdom_exponential}.
#' @param ... other arguments.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' library(ggplot2)
#'
#' data(spectra)
#'
#' fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#' p <- plot(fit)
#' p
#' p + ggtitle("My super fit")
plot.exponential_fit <- function(x, ...) {

  df <- data.frame(x = x$x, y = x$y, yy = predict(x))

  p <- ggplot(df, aes(x = x)) +
    geom_point(aes_string(y = "y")) +
    geom_line(aes_string(y = "yy"), col = "red")

  invisible(p)
}

#' Summary of a CDOM exponential fit
#'
#' @param object An object returned by \code{cdom_exponential}.
#' @param ... other arguments.
#'
#' @return A numerical vector with estimated coefficients.
#' @export
#'
#' @examples
#' data(spectra)
#'
#' fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
#' summary(fit)
summary.exponential_fit <- function(object, ...) {

  print(summary(object$model))
  cat("r2 = ", object$r2)

}
