extract_outlier <- function(x, y){

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
        return(y - predict(fit$model, newdata = list(x = x)))
        # return(list(x = xx, y = yy, model = fit$model))
      }

    } else {
      flag = FALSE
    }

  }

  return(y - predict(fit$model, newdata = list(x = x)))

  # return(list(x = xx, y = yy, model = fit$model))
}


func1 <- function(par, x, y, functooptimise){

  fo <- formula(sub(".*~", "~", deparse(functooptimise)))
  func <- gsubfn::fn$identity(fo)
  y. <- do.call(func, c(list(x = x), par))

  sum((y - y.)^2)

}

#' @export
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

find_segment <- function(df, merge = TRUE, min_distance) {

  spl <- smooth.spline(df$y ~ df$x)

  pred1 <- predict(spl, x = df$x, deriv = 1)
  deriv1 <- pred1$y

  rising <- c(ifelse(diff(deriv1) > 0, TRUE, FALSE), NA)

  myrle <- rle(rising)

  df$rising <- rep(myrle$values, times = myrle$lengths)
  df$segment <- rep(1:length(myrle$lengths), times = myrle$lengths)
  df$deriv1 <- deriv1

  if(merge){
    df <- merge_segment(df, min_distance = min_distance)
  }

  res <- df %>%
    filter(rising == FALSE) %>%
    group_by(segment) %>%
    nest() %>%
    mutate(integral = purrr::map(data, ~pracma::trapz(.$x, .$y))) %>%
    mutate(start_pos = purrr::map(data, ~min(.$x))) %>%
    mutate(end_pos = purrr::map(data, ~max(.$x))) %>%
    mutate(has_peak = purrr::map(data, ~has_peak(.)))

  res <- unnest(res, start_pos, end_pos, integral) %>%
    arrange(desc(integral)) %>%
    mutate(segment = 1:nrow(.)) # make sur this way of sorting data works

  return(res)
}

merge_segment <- function(segment, min_distance) {

  res <- segment %>%
    group_by(rising, segment) %>%
    nest() %>%
    mutate(end_pos = purrr::map(data, ~max(.$x))) %>%
    unnest(end_pos)

  d <- which((diff(res$end_pos[res$rising == FALSE]) < min_distance) == TRUE)

  # Nothing to merge
  if(length(d) == 0) {return(segment)}

  index <- purrr::map2(d[1:length(d) - 1], d[2:length(d)], function(x, y) x + y) %>%
    unlist()

  res2 <- res

  res2$rising[index] <- FALSE

  myrle <- rle(res2$rising)

  res2$segment <- rep(1:length(myrle$lengths), times = myrle$lengths)

  res3 <- unnest(res2) %>%
    select(-end_pos)

  return(res3)
}

find_initial <- function(spectra, segment) {

  ngaussian <- nrow(segment)

  # *************************************************************************
  # Standrad exponential parameters
  # *************************************************************************

  sf <- splinefun(spectra$x, spectra$y)

  a0 <- sf(350)           # estimate a0 at 350 nm
  K <- mean(sf(600:700))  # mean value between 600-700 nm
  S <- 0.02               # standard value of S

  starting_exp <- c("a0" = a0, "K" = K, "S" = S)

  # *************************************************************************
  # Build a named vector with starting parameters.
  #
  # If there is a "true" peak in the segment, use it to estimate initial
  # values.
  # *************************************************************************

  segment <- segment %>%
    mutate(p1 = purrr::map(data, ~.$x[which.min(abs(.$deriv1))])) %>%
    mutate(p0 = purrr::map(data, ~.$y_gaussian[which.min(abs(.$deriv1))])) %>%
    mutate(p2 = purrr::map(data, ~max(.$x) - min(.$x))) %>%
    unnest(p1, p0, p2) %>%
    filter(segment %in% 1:ngaussian)

  p0 <- segment$p0
  p1 <- segment$p1
  p2 <- segment$p2 / 2 # Find a better way

  starting_gaussian <- as.vector(rbind(p0, p1, p2))
  gaussian_name <- paste0(rep(c("p0", "p1", "p2"), time = ngaussian),
                          rep(letters[1:ngaussian], each = 3))

  starting_gaussian <- setNames(starting_gaussian, gaussian_name)

  # *************************************************************************
  # Bind everything together
  # *************************************************************************

  starting_value <- c(starting_exp, starting_gaussian)

  return(starting_value)
}

#' @importFrom gridExtra grid.arrange
plot_segment <- function(spectra, segment) {

  segment <- unnest(segment, data)

  p1 <- ggplot() +
    geom_point(data = spectra, aes(x = x, y = y), col = "black") +
    geom_point(data = segment, aes(x = x, y = y, col = factor(segment)))

  p2 <- ggplot() +
    geom_point(data = spectra, aes(x = x, y = y_gaussian), col = "black") +
    geom_point(data = segment, aes(x = x, y = y_gaussian, col = factor(segment)))

  p <- gridExtra::grid.arrange(p1, p2, ncol = 1)

  return(p)

}

#' Gaussian decomposition on CDOM spectra
#'
#' @param x Numerical vector of wavelengths.
#' @param y Numerical vector of absorption/absorbance values.
#' @param filter Logical. Should the spectra filtered? Default is TRUE.
#' @param min_distance Minimum distance in nm allowed between possible Gaussian
#'   components.
#'
#' @return A data frame with estimated parameters.
#' @export
#'
#' @importFrom  signal sgolayfilt
#' @import dplyr
#' @import minpack.lm
#' @importFrom  purrr map
#' @import proto
#'
#' @examples
#' \dontrun{fitCDOMcomponents(x = wl, y = spc, min_distance = 50)}
fitCDOMcomponents <- function(x, y, filter = TRUE, min_distance) {

  spectra <- dplyr::data_frame(x = x, y = y) %>%
    dplyr::mutate(y = signal::sgolayfilt(y, p = 3, n = 21)) %>%
    dplyr::mutate(y_gaussian = extract_outlier(x, y))

  # segment <- find_segment(spectra, merge = FALSE, min_distance = min_distance)
  # plot_segment(spectra, segment)

  segment <- find_segment(spectra, merge = TRUE, min_distance = min_distance)
  plot_segment(spectra, segment)

  # *************************************************************************
  # Find the "right" number of gaussian component and estimate starting
  # parameters for each of them.
  # *************************************************************************

  segment <- guess_ngaussian(segment)
  ngaussian <- nrow(segment)
  message(sprintf("Estimated number of components: %d\n", ngaussian))

  starting_values <- find_initial(spectra, segment)
  lower_values <- set_lower(starting_values)

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

  # bestguesses2 <- BB::spg(par = starting_values,
  #            fn = func1,
  #                        x = spectra$x,
  #                        y = spectra$y,
  #                        functooptimise = myfunc,
  #                        control = list(maxit = maxit),
  #                        lower = starting_values - starting_values,
  #                         upper = starting_values * 2,
  #                        method = 1,
  #                        quiet = TRUE)
  #
  # bestguesses <- optim(par = starting_values,
  #                      fn = func1 ,
  #                      x = spectra$x,
  #                      y = spectra$y,
  #                      functooptimise = myfunc,
  #                      control = list(maxit = maxit),
  #                      lower = starting_values - starting_values,
  #                      upper = starting_values * 2,
  #                      method = c("L-BFGS-B"))

  fit <- tryCatch(
    minpack.lm::nlsLM(formula = myfunc,
          start = starting_values,
          lower = lower_values,
          #upper = starting_values * 1.5,
          data = data.frame(y = spectra$y, x = spectra$x),
          control = list(maxiter = maxit, maxfev = 1000)),
    error = function(e) NULL,  warning = function(e) NULL)

  data_frame(params = names(starting_values), coef(fit), starting_values, lower_values)

  # fo <- formula(sub(".*~", "~", deparse(myfunc)))
  # func <- gsubfn::fn$identity(fo)
  # y. <- do.call(func, c(list(x = spectra$x), bestguesses$par))
  # y2. <- do.call(func, c(list(x = spectra$x), bestguesses2$par))

  plot(spectra$x, spectra$y)
  lines(spectra$x, predict(fit, newdata = list(x = spectra$x)), col = "red")
  #lines(spectra$x, y2., col = "green")

  res <- data_frame(
    parameter = names(starting_values),
    guess = starting_values,
    nls = coef(fit))

  return(res)

}

# *************************************************************************
# Determine if a segment has a peak in it.
# *************************************************************************

has_peak <- function(x) {
  peak <- pracma::findpeaks(x$y_gaussian)
  return(ifelse(is.null(peak), FALSE, TRUE))
}

# *************************************************************************
# All integral greater than the average should be modeled?
# *************************************************************************

guess_ngaussian <- function(segment) {

  segment <- filter(segment, integral >= mean(integral))
  return(segment)
  # ngaussian <- which(segment$integral > mean(segment$integral))
}

# *************************************************************************
# Set lower bounds to starting guesses.
# *************************************************************************
set_lower <- function(starting_values, segment) {

  lower_values <- starting_values

  # Minimum values for the exponential parameters
  index <- grepl("S", names(starting_values))
  lower_values[index] <- 0.001

  index <- grepl("K", names(starting_values))
  lower_values[index] <- -starting_values[index]

  index <- grepl("a0", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.05

  # Minimum values for Gaussian parameters
  index <- grepl("p0", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.95

  index <- grepl("p1", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.95

  index <- grepl("p2", names(starting_values))
  lower_values[index] <- starting_values[index] * 0.15

  return(lower_values)

}
