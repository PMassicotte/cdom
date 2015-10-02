\name{spectra}
\alias{spectra}
\docType{data}
\title{
CDOM absorbance data.
}
\description{
Simple absorbance used to test package's functions.
}
\usage{data("spectra")}
\format{
  A data frame with 1421 observations on the following 2 variables.
  \describe{
    \item{\code{wavelength}}{a numeric vector}
    \item{\code{absorbance}}{a numeric vector}
  }
}

\examples{
data(spectra)
plot(spectra$wavelength, spectra$absorbance)
}
\keyword{datasets}
