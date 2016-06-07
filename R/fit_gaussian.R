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


  res <- df %>%
    filter(rising == FALSE) %>%
    group_by(segment) %>%
    nest() %>%
    mutate(integral = purrr::map(data, ~pracma::trapz(.$x, .$y))) %>%
    mutate(start_pos = purrr::map(data, ~min(.$x))) %>%
    mutate(end_pos = purrr::map(data, ~max(.$x))) %>%
    mutate(has_peak = purrr::map(data, ~has_peak(.))) %>%
    mutate(y_max = purrr::map(data, ~max(.$y_gaussian))) %>%
    unnest(start_pos, end_pos, integral) %>%
    arrange(desc(integral)) %>%
    mutate(segment = 1:nrow(.)) %>% # make sur this way of sorting data works
    filter(integral != 0) # drop segments with juste 1 point (integral = 0)

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

  index <- res[res$rising == TRUE, ]$segment[d]
  res$rising[index] <- FALSE

  myrle <- rle(res$rising)

  res$segment <- rep(1:length(myrle$lengths), times = myrle$lengths)

  res2 <- unnest(res) %>%
    select(-end_pos)

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

  params <- list(
    p0 = sf(mid_point),
    p1 = mid_point,
    p2 = (max(seg$x) - min(seg$x)) / 2
    )

  return(params)

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
fitCDOMcomponents <- function(x,
                              y,
                              filter = TRUE,
                              min_distance = 50,
                              min_height = 1) {

  spectra <- dplyr::data_frame(x = x, y = y)
  spectra$y <- signal::sgolayfilt(y, p = 3, n = 21)

  out <- extract_outlier(x, y)

  spectra$y_gaussian <- out$outlier

  # segment <- find_segment(spectra, merge = FALSE, min_distance = min_distance)
  # plot_segment(spectra, segment)

  segment <- find_segment(spectra, merge = TRUE, min_distance = min_distance)
  plot_segment(spectra, segment)

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
          start = starting_values$start,
          lower = starting_values$lower,
          #upper = starting_values$upper,
          data = data.frame(y = spectra$y, x = spectra$x),
          control = list(maxiter = maxit, maxfev = 1000)),
    error = function(e) NULL,  warning = function(e) NULL)

  starting_values$estimated <- coef(fit)

  # fo <- formula(sub(".*~", "~", deparse(myfunc)))
  # func <- gsubfn::fn$identity(fo)
  # y. <- do.call(func, c(list(x = spectra$x), bestguesses$par))
  # y2. <- do.call(func, c(list(x = spectra$x), bestguesses2$par))

  plot(spectra$x, spectra$y)
  lines(spectra$x, predict(fit, newdata = list(x = spectra$x)), col = "red")
  #lines(spectra$x, y2., col = "green")


  return(starting_values)

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

  segment <- filter(segment, y_max >= min_height)
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
