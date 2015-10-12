#' CDOM absorbance data.
#'
#' Simple absorbance spectra used to test package's functions.
#'
#' \itemize{
#'   \item wavelength.  Wavelengths used for measurements (190-900 nm.)
#'   \item Absorbance.
#' }
#'
#' @import ggplot2
#' @import tidyr
#' @docType data
#' @keywords datasets
#' @name spectra
#' @usage data(spectra)
#' @format A data frame with 711 rows and 26 variables
#' @examples
#' library(ggplot2)
#' library(tidyr)
#' data("spectra")
#' spectra <- gather(spectra, sample, absorbance, -wavelength)
#' ggplot(spectra, aes(x = wavelength, y = absorbance, group = sample)) +
#'  geom_line(size = 0.1)
NULL
