% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exponential.R
\name{coef.exponential_fit}
\alias{coef.exponential_fit}
\title{Extract Model Coefficients from a CDOM exponential fit}
\usage{
\method{coef}{exponential_fit}(object, ...)
}
\arguments{
\item{object}{An object returned by \code{cdom_exponential}.}

\item{...}{other arguments.}
}
\value{
A numerical vector with estimated coefficients.
}
\description{
Extract Model Coefficients from a CDOM exponential fit
}
\examples{
data(spectra)

fit <- cdom_exponential(spectra$wavelength, spectra$spc1, 350, 190, 900)
coef(fit)
}
