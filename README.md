<!-- badges: start -->
[![R-CMD-check](https://github.com/CeMOS-Mannheim/MALDIcellassay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CeMOS-Mannheim/MALDIcellassay/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/github/CeMOS-Mannheim/MALDIcellassay/graph/badge.svg?token=URVX29WTDX)](https://codecov.io/github/CeMOS-Mannheim/MALDIcellassay)
[![CRAN status](https://www.r-pkg.org/badges/version/MALDIcellassay)](https://CRAN.R-project.org/package=MALDIcellassay)
<!-- badges: end -->

# MALDIcellassay
Detects high variance signals and generates dose-response curves to further investigate candidate signals from MALDI cell based assays.
More information in [Unger et. Al. 2021](https://doi.org/10.1038/s41596-021-00624-z)

## Main functionality
The main function in this package is `fitCurve()` which will not only do logistic regression and variance filtering but also handle all preprocessing necessary:
- Spectral alignment ("single-point re-calibration").
- Normalization by an internal standard or known endogenous signal

All these function can also be called on there own.

This package makes heavy use of the `MALDIquant` framework.

## GUI
If you want to use `MALDIcellassay` interactivly to explore your data you might be interested in [`M²ara`](https://github.com/CeMOS-Mannheim/M2ara) which offers a Shiny based interface.
It also comes as a [docker container](https://hub.docker.com/repository/docker/thomasenzlein/m2ara) and as a stand-alone [Windows installer](https://github.com/CeMOS-Mannheim/m2ara/releases/latest/).
Also check the [preprint](https://chemrxiv.org/engage/chemrxiv/article-details/663a1d0f418a5379b0aa286b) on M²ara. 
