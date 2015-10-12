\name{spectra}
\alias{spectra}
\docType{data}
\title{
CDOM absorbance data.
}
\description{
Simple absorbance spectra used to test package's functions.
}
\usage{data("spectra")}
\format{
  A data frame with 711 observations on the following 2 variables.
  \describe{
    \item{\code{wavelength}}{a numeric vector (190-900 nm.)}
    \item{\code{absorbance}}{a numeric vector}
  }
}

\examples{
data(spectra)
plot(spectra$wavelength, spectra$absorbance)
}
\keyword{datasets}
