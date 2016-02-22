
cdom [![Travis-CI Build Status](https://api.travis-ci.org/PMassicotte/cdom.svg?branch=master)](https://travis-ci.org/PMassicotte/cdom) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/PMassicotte/cdom?branch=master&svg=true)](https://ci.appveyor.com/project/PMassicotte/cdom) [![Package-License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

<!-- [![CRAN](http://www.r-pkg.org/badges/version/cdom)](http://cran.rstudio.com/package=cdom) [![Downloads](http://cranlogs.r-pkg.org/badges/cdom?color=brightgreen)](http://www.r-pkg.org/pkg/cdom) -->
The **cdom** package implements various functions used to model and calculate metrics from absorption spectra of chromophotic dissolved organic matter (CDOM).

This package provides:

1.  Simple wrappers to calculate common metrics found in the literature.
    -   The **spectral curve** (Loiselle et al. 2009).
    -   The **slope ratio (Sr)** (Helms et al. 2008).
    -   The **spectral slope (S)** (Jerlov 1968; Lundgren 1976; Bricaud, Morel, and Prieur 1981).

2.  The function to use the **Gaussian decomposition approach** proposed in Massicotte and Markager, (2015).

The package can be installed using the following command.

``` r
devtools::install_github("PMassicotte/cdom")
```

Please note that this is a developing version of the package for testing only. Please fill an issue when you find bugs.

All functions from the package start with the `cdom_` prefix.

``` r
library(cdom)
ls("package:cdom")
## [1] "cdom_fit_exponential" "cdom_slope_ratio"     "cdom_spectral_curve" 
## [4] "spectra"
```

Examples
========

The spectral slope (S)
----------------------

The `cdom_fit_exponential()` function fits an exponential curve to CDOM data using the simple model proposed by Jerlov (1968), Lundgren (1976), Bricaud, Morel, and Prieur (1981).

``` tex
a(\lambda) = a(\lambda0)e^{-S(\lambda - \lambda0)} + K
```

``` r
library(ggplot2)
library(cdom)
data("spectra")

fit <- cdom_fit_exponential(wl = spectra$wavelength,
                       absorbance = spectra$spc3,
                       wl0 = 350,
                       startwl = 190,
                       endwl = 900)

ggplot(spectra, aes(x = wavelength, y = spc3)) +
  geom_point() +
  geom_line(aes(y = fit$data$.fitted), col = "red") +
  xlab("Wavelength (nm)") +
  ylab(expression(paste("Absorption (", m ^ {-1}, ")")))
```

![](inst/images/README-exponential-1.png)<!-- -->

The slope ratio (SR)
--------------------

The `cdom_slope_ratio()` function calculates the slope ratio (S<sub>R</sub>) which is defined as: S<sub>275-295</sub>/S<sub>350-400</sub>. See Helms et al. (2008) for detailed information.

``` r
library(cdom)
data("spectra")

cdom_slope_ratio(spectra$wavelength, spectra$spc1)
## [1] 1.325082
```

The spectral curve
------------------

The `cdom_spectral_curve()` function generates the spectral curve using the slope of the linear regression between the natural log absorption spectrum and wavelengths over a sliding window of 21 nm interval (default) at 1 nm resolution. See Loiselle et al. (2009) for detailed information.

``` r
library(cdom)
data("spectra")

res <-  cdom_spectral_curve(wl = spectra$wavelength,
                       absorbance = spectra$spc10,
                       interval = 21,
                       r2threshold = 0.98) # Maybe to restrictive...

ggplot(res, aes(x = wl, y = s)) +
  geom_point() +
  geom_line() +
  xlab("Wavelength (nm)") +
  ylab(expression(paste("Spectral slope (", nm ^ {-1}, ")")))
```

![](inst/images/README-spectral_curve-1.png)<!-- -->

Data
====

A total 25 absorption spectra are provided in the package.

``` r
library(ggplot2)
library(tidyr)
data("spectra")

spectra <- gather(spectra, sample, absorption, -wavelength)

ggplot(spectra, aes(x = wavelength, y = absorption, group = sample)) +
  geom_line(size = 0.1) +
  xlab("Wavelength (nm)") +
  ylab(expression(paste("Absorption (", m ^ {-1}, ")")))
```

![](inst/images/README-data-1.png)<!-- -->

How to cite the package
=======================

``` r
citation("cdom")
## 
## To cite cdom in publications use:
## 
##   Massicotte, P., and Markager, S. (2016). Using a Gaussian
##   decomposition approach to model absorption spectra of
##   chromophoric dissolved organic matter. Mar. Chem. 180, 24-32.
##   doi:10.1016/j.marchem.2016.01.008.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Using a Gaussian decomposition approach to model absorption spectra of chromophoric dissolved organic matter},
##     author = {Philippe Massicotte and Stiig Markager},
##     journal = {Marine Chemistry},
##     year = {2016},
##     volume = {180},
##     pages = {24--32},
##     url = {http://linkinghub.elsevier.com/retrieve/pii/S0304420316300081},
##   }
```

References
==========

Bricaud, Annick, André Morel, and Louis Prieur. 1981. “Absorption by dissolved organic matter of the sea (yellow substance) in the UV and visible domains.” *Limnology and Oceanography* 26 (1): 43–53. doi:[10.4319/lo.1981.26.1.0043](https://doi.org/10.4319/lo.1981.26.1.0043).

Helms, John R., Aron Stubbins, Jason D. Ritchie, Elizabeth C. Minor, David J. Kieber, and Kenneth Mopper. 2008. “Absorption spectral slopes and slope ratios as indicators of molecular weight, source, and photobleaching of chromophoric dissolved organic matter.” *Limnology and Oceanography* 53 (3): 955–69. doi:[10.4319/lo.2008.53.3.0955](https://doi.org/10.4319/lo.2008.53.3.0955).

Jerlov, N.G. 1968. *Optical oceanography*. New York: Elsevier Publishing Company.

Loiselle, Steven A., Luca Bracchini, Arduino M. Dattilo, Maso Ricci, Antonio Tognazzi, Andres Cézar, and Claudio Rossi. 2009. “The optical characterization of chromophoric dissolved organic matter using wavelength distribution of absorption spectral slopes.” *Limnology and Oceanography* 54 (2): 590–97. doi:[10.4319/lo.2009.54.2.0590](https://doi.org/10.4319/lo.2009.54.2.0590).

Lundgren, Bo. 1976. “Spectral transmittance measurements in the Baltic.” Copenhagen: Institute Physical Oceanography University of Copenhagen.
