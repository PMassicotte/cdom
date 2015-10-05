CDOM
====

The **cdom** package implements various functions used to model and calculate metrics from absorption spectra of chromophotic dissolved organic matter (CDOM).

There are two main goals:

1.  Provides simple wrappers to calculate commonly metrics found in the literature.
    -   The **spectral slope** (Loiselle et al., 2009).
    -   The **slope ratio** (Helms et al., 2008).
    -   The well-known spectral slope (**S**) proposed by Bricaud et al., 1981.

2.  The **Gaussian decomposition approach** proposed in MAssicotte and Markager, 2015.

Please note that this is a developing version of the package for testing only. Please fill an issue if you find bugs.

Examples
========

Calculate the *standard* spectral slope.
----------------------------------------

The `fit_exponential()` function fit an exponential curve to CDOM data using the simple model proposed by Bricaud et al. 1981.

\[
a(\lambda) = a(\lambda0)e^{-S(\lambda - \lambda0)} + K
\]

``` r
library(cdom)
data("spectra")

plot(spectra$wavelength, spectra$absorbance)

fit <- fit_exponential(wl = spectra$wavelength,
                       spectra = spectra$absorbance,
                       wl0 = 350,
                       startwl = 190,
                       endwl = 900)

lines(spectra$wavelength, predict(fit), col = "red")
```

![](README-unnamed-chunk-2-1.png)

### Calculate the slope ratio (SR)

``` r
library(cdom)
data("spectra")

slope_ratio(spectra$wavelength, spectra$absorbance)
#> wl_275_295 
#>   6.838766
```

References
----------

Helms, John R., Aron Stubbins, Jason D. Ritchie, Elizabeth C. Minor, David J. Kieber, and Kenneth Mopper. 2008. “Absorption Spectral Slopes and Slope Ratios as Indicators of Molecular Weight, Source, and Photobleaching of Chromophoric Dissolved Organic Matter.” Limnology and Oceanography 53 (3): 955–69. <doi:10.4319/lo.2008.53.3.0955>.

Loiselle, Steven A., Luca Bracchini, Arduino M. Dattilo, Maso Ricci, Antonio Tognazzi, Andres Cézar, and Claudio Rossi. 2009. “The Optical Characterization of Chromophoric Dissolved Organic Matter Using Wavelength Distribution of Absorption Spectral Slopes.” Limnology and Oceanography 54 (2): 590–97. <doi:10.4319/lo.2009.54.2.0590>.

Massicotte, Philippe, and Stiig Markager. 2015. “Using a Gaussian Decomposition Approach to Model Absorption Spectra of Chromophoric Dissolved Organic Matter.”
