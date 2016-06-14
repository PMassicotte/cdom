#' @importFrom stats residuals
extract_outlier <- function(x, y){

  xx <- x
  yy <- y

  fit <- cdom_exponential(wl = xx,
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

    toremove <- which(residuals(fit$model) > 1 * mean(abs(residuals(fit$model))))

    ## Do not remove points at the very low wavelength
    toremove <- toremove[which(toremove - seq(1:length(toremove)) > 0)]


    #points(xx[toremove], yy[toremove], pch = 21, col = "red")

    if(length(toremove) > 0) {

      xx <- xx[-toremove]
      yy <- yy[-toremove]

      fit <- cdom_exponential(wl = xx,
                                  absorbance = yy,
                                  wl0 = 350,
                                  startwl = min(xx),
                                  endwl = max(xx))

      if(is.null(fit)) {
        return(NULL)
      }

      ## Do not attempt to model data if not "enough" observations.
      if(length(xx) < 50) {

        outlier <- y - predict(fit$model, newdata = list(x = x))

        exponential_coef <- setNames(as.vector(coef(fit$model)),
                                     names(coef(fit$model)))

        df <- list(outlier = outlier , exponential_coef = exponential_coef)

        return(df)
        # return(list(x = xx, y = yy, model = fit$model))
      }

    } else {
      flag = FALSE
    }

  }

  outlier <- y - predict(fit$model, newdata = list(x = x))

  exponential_coef <- setNames(as.vector(coef(fit$model)),
                               names(coef(fit$model)))

  df <- list(outlier = outlier , exponential_coef = exponential_coef)

  return(df)

}


func1 <- function(par, x, y, functooptimise){

  fo <- formula(sub(".*~", "~", deparse(functooptimise)))
  func <- gsubfn::fn$identity(fo)
  y. <- do.call(func, c(list(x = x), par))

  sum((y - y.)^2)

}


#' Build Gaussian Model
#'
#' @param ngaussian Numerical. Number of Gaussian component to consider.
#'
#' @return A \code{formula} object.
#' @export
#'
#' @examples
#' myfunction <- build_model(2)
build_model <- function(ngaussian = 1) {

  exponential_part <- "y ~ a0 * exp(-S * (x - 350)) + K +"

  gaussian_part <- "(p0%s*exp(-((x-p1%s)^2/(2*p2%s^2))))"

  ## Create the equation with the appropriate number of Gaussian components
  gaussian_part <- rep(gaussian_part, ngaussian)
  gaussian_part <- paste(gaussian_part, collapse = "+")

  gaussian_part <-
    do.call(sprintf, as.list(c(gaussian_part, rep(letters[1:ngaussian],
                                         each = 3,
                                         times = 1))))

  myfunc <- as.formula(paste(exponential_part, gaussian_part))

  return(myfunc)
}

gaussian <- function(df, x) {

  index_p0 <- which(df$param == "p0")
  index_p1 <- which(df$param == "p1")
  index_p2 <- which(df$param == "p2")

  res <- df$estimated[index_p0] * exp(-((x - df$estimated[index_p1])^2 / (2 * df$estimated[index_p2]^2)))

  return(res)
}



find_segment <- function(df, merge = TRUE, min_distance) {

  spl <- smooth.spline(df$y_gaussian ~ df$x)

  pred1 <- predict(spl, x = df$x, deriv = 1)
  deriv1 <- pred1$y

  df$rising <- c(ifelse(diff(deriv1) > 0, TRUE, FALSE), NA)

  myrle <- rle(df$rising)


  df$segment <- rep(1:length(myrle$lengths), times = myrle$lengths)
  df$deriv1 <- deriv1

  if(merge){
    df <- merge_segment(df, min_distance = min_distance)
  }


  res <- filter_(df, ~rising == FALSE)
  res <- group_by_(res, "segment")
  res <- nest_(res, "data")
  res <- mutate_(res, integral = ~purrr::map(data, ~pracma::trapz(.$x, .$y)))
  res <- mutate_(res, start_pos = ~purrr::map(data, ~min(.$x)))
  res <- mutate_(res, end_pos = ~purrr::map(data, ~max(.$x)))
  res <- mutate_(res, has_peak = ~purrr::map(data, ~has_peak(.)))
  res <- mutate_(res, y_max = ~purrr::map(data, ~max(.$y_gaussian)))
  res <- unnest_(res, c("start_pos", "end_pos", "integral"))
  res <- arrange_(res, ~desc(integral))
  res <- mutate_(res, segment = ~1:nrow(res)) # make sur this way of sorting data works
  res <- filter_(res, ~integral != 0) # drop segments with juste 1 point (integral = 0)

  return(res)
}

merge_segment <- function(segment, min_distance) {

  res <- group_by_(segment, c("rising", "segment"))
  res <- nest_(res, "data")
  res <- mutate_(res, end_pos = ~purrr::map(data, ~max(.$x)))
  res <- unnest_(res, "end_pos")

  d <- which((diff(res$end_pos[res$rising == FALSE]) < min_distance) == TRUE)

  # Nothing to merge
  if(length(d) == 0) {return(segment)}

  index <- res[res$rising == TRUE, ]$segment[d]
  res$rising[index] <- FALSE

  myrle <- rle(res$rising)

  res$segment <- rep(1:length(myrle$lengths), times = myrle$lengths)

  res2 <- unnest_(res)
  res2 <- select_(res2, "-end_pos")

  return(res2)
}

find_initial_gaussian <- function(spectra, segment) {

  ngaussian <- nrow(segment)

  # *************************************************************************
  # Build a named vector with starting parameters.
  #
  # If there is a "true" peak in the segment, use it to estimate initial
  # values.
  # *************************************************************************

  params <- lapply(segment$data, p)

  starting_gaussian <- as.vector(unlist(params))

  gaussian_name <- paste0(rep(c("p0", "p1", "p2"), time = ngaussian),
                          rep(letters[1:ngaussian], each = 3))

  starting_gaussian <- setNames(starting_gaussian, gaussian_name)

  return(starting_gaussian)
}

p <- function(seg) {

  # Is the segment has a peak, then use it to estimate starting parameters
  if (has_peak(seg)) {

    params <- list(
      p0 = max(seg$y_gaussian),
      p1 = seg$x[which.max(seg$y_gaussian)],
      p2 = (max(seg$x) - min(seg$x)) / 2
      )

    return(params)

  }

  # If the segment has no peak, then use the point in the middle of the
  # segment to estimate parameters
  mid_point <-  min(seg$x) + (max(seg$x) - min(seg$x)) / 2

  sf <- splinefun(seg$x, seg$y_gaussian)

  params <- list(p0 = sf(mid_point),
                 p1 = mid_point,
                 p2 = (max(seg$x) - min(seg$x)) / 2)

  return(params)

}

#' @importFrom gridExtra grid.arrange
plot_segment <- function(spectra, segment) {

  segment <- unnest(segment, data)

  p1 <- ggplot() +
    geom_point(data = spectra, aes_string(x = "x", y = "y"), col = "black") +
    geom_point(data = segment, aes_string(x = "x", y = "y", col = factor(segment)))

  p2 <- ggplot() +
    geom_point(data = spectra, aes_string(x = "x", y = "y_gaussian"), col = "black") +
    geom_point(data = segment, aes_string(x = "x", y = "y_gaussian", col = factor(segment)))

  p <- gridExtra::grid.arrange(p1, p2, ncol = 1)

  return(p)

}

#' Gaussian decomposition on CDOM spectra
#'
#' @param x Numerical vector of wavelengths.
#' @param y Numerical vector of absorption/absorbance values.
#' @param filter Logical. Should the spectra filtered? Default is TRUE.
#' @param min_distance Minimum distance in nm allowed between possible Gaussian
#'   components. Default is 50 nm.
#' @param merge Logical. Determines if initially estimated Gaussian components
#'   should be merged.
#' @param min_height Numerical. Determines the minimum height (per meter) to
#'   model Gaussian components. This is usefull to set a threshold to avoid
#'   modeling noise at higher wavelenghts. Default is 1.
#'
#' @return A data frame with estimated parameters.
#' @export
#'
#' @importFrom  signal sgolayfilt
#' @importFrom stats as.formula formula setNames smooth.spline
#' @importFrom utils data
#' @import dplyr
#' @import minpack.lm
#' @importFrom  purrr map
#' @import proto
#'
#' @examples
#' \dontrun{cdom_gaussian(x = wl, y = spc, min_distance = 50)}
cdom_gaussian <- function(x,
                          y,
                          filter = TRUE,
                          merge = TRUE,
                          min_distance = 50,
                          min_height = 1) {

  spectra <- dplyr::data_frame(x = x, y = y)
  spectra$y <- signal::sgolayfilt(y, p = 3, n = 21)

  out <- extract_outlier(x, y)

  spectra$y_gaussian <- out$outlier

  # segment <- find_segment(spectra, merge = FALSE, min_distance = min_distance)
  # plot_segment(spectra, segment)

  segment <- find_segment(spectra, merge = merge, min_distance = min_distance)
  # plot_segment(spectra, segment)

  # *************************************************************************
  # Find the "right" number of gaussian component and estimate starting
  # parameters for each of them.
  # *************************************************************************

  segment <- guess_ngaussian(segment, min_height = min_height)
  ngaussian <- nrow(segment)
  message(sprintf("Estimated number of components: %d\n", ngaussian))

  starting_values <- find_initial_gaussian(spectra, segment)
  starting_values <- c(out$exponential_coef, starting_values)

  bound_values <- set_bounds(starting_values)

  starting_values <- data_frame(param = names(starting_values),
                                start = starting_values,
                                lower = bound_values$lower_values,
                                upper = bound_values$upper_values)

  # *************************************************************************
  # Build model equation.
  # *************************************************************************

  myfunc <- build_model(ngaussian = nrow(segment))

  # print(starting_values)
  # print(myfunc)

  # *************************************************************************
  # Start the optimisation process.
  # *************************************************************************

  maxit <- 1000

  # bestguesses <- optim(par = starting_values$start,
  #                      fn = func1 ,
  #                      x = spectra$x,
  #                      y = spectra$y,
  #                      functooptimise = myfunc,
  #                      control = list(maxit = maxit),
  #                      lower = starting_values$lower,
  #                      #upper = starting_values * 2,
  #                      method = c("L-BFGS-B"))
  #
  # starting_values$estimated <- bestguesses$par

  safe_nls <- purrr::safely(minpack.lm::nlsLM)

  fit <- safe_nls(formula = myfunc,
                  start = starting_values$start,
                  lower = starting_values$lower,
                  # upper = starting_values$upper,
                  data = data.frame(y = spectra$y, x = spectra$x),
                  control = list(maxiter = maxit, maxfev = 1000))

  # A valide nls model as been returned
  if (is.null(fit$error)) {

    r2 <- 1 - sum((spectra$y - predict(fit$result))^2) / (length(spectra$y) * var(spectra$y))

    res <- list(model = fit$result, x = spectra$x, y = spectra$y, r2 = r2)
    class(res) <- "gaussian_fit"

    return(res)

  } else {
    return(NULL)
  }

  # fit <- tryCatch(
  #   minpack.lm::nlsLM(formula = myfunc,
  #         start = starting_values$start,
  #         lower = starting_values$lower,
  #         # upper = starting_values$upper,
  #         data = data.frame(y = spectra$y, x = spectra$x),
  #         control = list(maxiter = maxit, maxfev = 1000)),
  #   error = function(e) NULL,  warning = function(e) NULL)

  # starting_values$estimated <- coef(fit)

  # bestguesses <- optimx::optimx(par = starting_values$start,
  #                               fn = func1 ,
  #                               x = spectra$x,
  #                               y = spectra$y,
  #                               functooptimise = myfunc,
  #                               method = "nmkb")
  #
  # starting_values$estimated <- unlist(bestguesses[1:nrow(starting_values)])

  # res <- list(params = starting_values)
  #
  # class(res) <- "cdom_gaussian"
  #
  # return(res)

}

# *************************************************************************
# Determine if a segment has a peak in it.
# *************************************************************************

has_peak <- function(x) {
  peak <- pracma::findpeaks(x$y_gaussian)
  return(ifelse(is.null(peak), FALSE, TRUE))
}

# *************************************************************************
# Keep only peaks with a minimum height.
# *************************************************************************

guess_ngaussian <- function(segment, min_height = 1) {

  segment <- filter_(segment, ~y_max >= min_height)
  return(segment)

}

# *************************************************************************
# Set lower bounds to starting guesses.
# *************************************************************************

set_bounds <- function(starting_values) {

  # *************************************************************************
  # Lower values
  # *************************************************************************

  lower_values <- starting_values

  # Minimum values for the exponential parameters
  index <- grepl("S", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.5

  index <- grepl("K", names(starting_values))
  lower_values[index] <- starting_values[index] / 5

  index <- grepl("a0", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.5

  # Minimum values for Gaussian parameters
  index <- grepl("p0", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.5

  index <- grepl("p1", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.5

  index <- grepl("p2", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.15

  # *************************************************************************
  # Upper values
  # *************************************************************************

  upper_values <- starting_values

  # Maximum values for the exponential parameters
  index <- grepl("S", names(starting_values))
  upper_values[index] <- starting_values[index] * 1.5

  index <- grepl("K", names(starting_values))
  upper_values[index] <- starting_values[index] * 5

  index <- grepl("a0", names(starting_values))
  upper_values[index] <- starting_values[index] * 1.5

  # Minimum values for Gaussian parameters
  index <- grepl("p0", names(starting_values))
  upper_values[index] <- starting_values[index] * 1.5

  index <- grepl("p1", names(starting_values))
  upper_values[index] <- starting_values[index] * 1.5

  index <- grepl("p2", names(starting_values))
  upper_values[index] <- starting_values[index] * 1.15

  df <- data_frame(lower_values, upper_values)

  return(df)

}

#' Gaussian Model Predictions
#'
#' @param object a \code{cdom_gaussian} model object for which prediction is
#'   desired.
#' @param ...	other arguments.
#'
#' @return A data frame containing
#'  \describe{
#'   \item{exponential_part}{The fitted exponential part of the model}
#'   \item{component1 ... componentn}{The individual Gaussian components}
#'   \item{spectra}{The fitted spectra}
#' }
#' @export
#'
#' @examples
#' data(spectra)
#' myfit <- cdom_gaussian(spectra$wavelength,
#'                        spectra$spc21,
#'                        min_distance = 50,
#'                        min_height = 10)
#'
#' predict(myfit)
predict.gaussian_fit <- function(object, ...) {

  res <- extract_components(object)

  return(res)

}

#' Extract Model Coefficients
#'
#' @param object An object returned by \code{cdom_gaussian()}.
#' @param ...	other arguments.
#'
#' @return Coefficients extracted from the model object \code{object}.
#' @export
#'
#' @examples
#' data(spectra)
#' myfit <- cdom_gaussian(spectra$wavelength,
#'                        spectra$spc21,
#'                        min_distance = 50,
#'                        min_height = 10)
#'
#' coef(myfit)
coef.gaussian_fit <- function(object, ...) {

  coef(object$model)

}

#' Plot CDOM Gaussian components
#'
#' @param x An object returned by the \code{cdom_gaussian()} function.
#' @param ...	other arguments.
#'
#' @return A ggplot2 plot.
#' @export
#' @importFrom gridExtra grid.arrange
#'
#' @examples
#' data(spectra)
#' myfit <- cdom_gaussian(spectra$wavelength,
#'                        spectra$spc21,
#'                        min_distance = 50,
#'                        min_height = 10)
#' plot(myfit)

plot.gaussian_fit <- function(x, ...) {

  df <- data.frame(x = x$x, y = x$y)

  components <- predict(x)
  components$x <- x$x

  # Plot showing raw data and the fitted curve
  p1 <- ggplot(df, aes_string(x = "x", y = "y")) +
    geom_point() +
    geom_line(data = components, aes_string(x = "x", y = "spectra"), col = "red")

  # Plot showing the Gaussian components
  n <- names(components)[grepl("component", names(components))]
  df <- gather_(components, "component", "value", n)

  p2 <- ggplot(df, aes_string(x = "x", y = "value")) +
    geom_line(aes_string(color = "component"))

  p <- gridExtra::grid.arrange(p1, p2)

  invisible(p)

}

extract_components <- function(x) {

  params <- coef(x)

  # Extract the exponential part
  a0 <- params[names(params) == "a0"]
  S <- params[names(params) == "S"]
  K <- params[names(params) == "K"]

  exponential_part <- a0 * exp(-S * (x$x - 350)) + K
  exponential_part <- as.data.frame(exponential_part)

  # Extract the Gaussiance components
  ngaussian <- (length(params) - 3) / 3

  df <- data.frame(param = names(coef(x)[grepl("p\\d", names(coef(x)))]),
                   estimated = coef(x)[grepl("p\\d", names(coef(x)))])
  df <- tidyr::separate_(df, "param", c("param", "component"), 2)
  df <- group_by_(df, "component")
  df <- nest(df)

  components <- purrr::map(df$data, ~ gaussian(., x = x$x))

  components <- do.call(cbind, components)
  components <- as.data.frame(components)

  res <- bind_cols(exponential_part, components)
  res$spectra <- rowSums(res)

  names(res) <- c("exponential_part", paste0("component", 1:ngaussian), "spectra")

  return(res)
}
