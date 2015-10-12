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
library(ggplot2)
library(tidyr)

data("spectra")

spectra <- gather(spectra, sample, absorbance, -wavelength)

ggplot(spectra, aes(x = wavelength, y = absorbance, group = sample)) +
  geom_line(size = 0.1)
}
\keyword{datasets}
