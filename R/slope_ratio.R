#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# FILE:         slope_ratio.R
#
# AUTHOR:       Philippe Massicotte
#
# DESCRIPTION:  Function to calculate the slope ratio (SR) proposed by Helms.
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

#' Calculate the slope ratio (SR) from an absorption spectra.
#'
#' @details Calculate the slope ratio (SR) as defined by Helms et al. (2008).
#'
#' \deqn{SR = \frac{S_{275-295}}{S_{350-400}}}
#'
#' @references \url{http://www.aslo.org/lo/toc/vol_53/issue_3/0955.html}
#'
#' @param wl The wavelength vector.
#' @param spectra The spectra vector.
#'
#' @return The value of the slope ratio.
#' @export
#'
#' @examples
#' data("spectra")
#' slope_ratio(spectra$wavelength, spectra$absorbance)
#'

slope_ratio <- function(wl, spectra) {

  if(length(wl) != length(spectra)){
    stop("wl and spectra are not of the same length.")
  }

  if(!is.numeric(wl) | !is.numeric(spectra)){
    stop("wl and spectra need to be numeric.")
  }

  #--------------------------------------------
  # Get data
  #--------------------------------------------
  sf <- splinefun(wl, spectra)

  wl_275_295 <- seq(from = 275, to = 295, length.out = 25)
  wl_350_400 <- seq(from = 350, to = 400, length.out = 25)

  data_275_295 <- sf(wl_275_295)
  data_350_400 <- sf(wl_350_400)

  #--------------------------------------------
  # Calculate the ratio.
  #--------------------------------------------
  slope_275_295 <- coef(lm(log(data_275_295) ~ wl_275_295))[2]
  slope_350_400 <- coef(lm(log(data_350_400) ~ wl_350_400))[2]

  #--------------------------------------------
  # Return the result.
  #--------------------------------------------
  sr <- slope_275_295 / slope_350_400

  return(sr)
}
