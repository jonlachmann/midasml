[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/midasml)](https://cran.r-project.org/package=midasml)
[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/midasml)](https://cran.rstudio.com/web/packages/midasml/index.html) 
[![Downloads](http://cranlogs.r-pkg.org/badges/midasml)](http://www.r-pkg.org/pkg/midasml)
# midasml

midasml - Estimation and Prediction Methods for High-Dimensional Mixed Frequency Time Series and Panel Data

## About

The *midasml* package implements *estimation* and *prediction* methods for high-dimensional mixed-frequency (MIDAS) time-series and panel data regression models. The regularized MIDAS models are estimated using orthogonal (e.g. Legendre) polynomials and sparse-group LASSO estimator. For more information on the *midasml* approach see [^1][^2][^3]. 

The package is equipped with the fast implementation of the sparse-group LASSO estimator by means of proximal block coordinate descent. High-dimensional mixed frequency time-series data can also be easily manipulated with functions provided in the package.

## Software in other languages

- Julia implmentation of the midasml method is available [here](https://github.com/ababii/Pythia.jl).
- MATLAB implmentation of the midasml method is available [here](https://github.com/jstriaukas/midasml_mat).
- Python implmentation of the midasml method is being developed at [here](https://github.com/jstriaukas/midasmlpy).

## Run to install the package

```{r }
# CRAN version - 0.1.10
install.packages("midasml") 

# Development version - 0.1.10
# install.packages("devtools")
library(devtools)
install_github("jstriaukas/midasml")
```
## Acknowledgements

Jonas Striaukas acknowledges that this material is based upon work supported by the Fund for Scientific Research-FNRS (Belgian National Fund for Scientific Research) under Grant #FC21388.

## References

[^1]: Babii, A., Ghysels, E., & Striaukas, J. Machine learning time series regressions with an application to nowcasting, (2022) *Journal of Business & Economic Statistics*, Volume 40, Issue 3, 1094-1106. https://doi.org/10.1080/07350015.2021.1899933. 

[^2]: Babii, A., Ghysels, E., & Striaukas, J. High-dimensional Granger causality tests with an application to VIX and news, (2022) *Journal of Financial Econometrics*, Forthcoming.

[^3]: Babii, A., R. Ball, Ghysels, E., & Striaukas, J. Machine learning panel data regressions with heavy-tailed dependent data: Theory and application, (2022) *Journal of Econometrics*, Forthcoming.
