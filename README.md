<!-- badges: start -->

[![R-CMD-check](https://github.com/CeMOS-Mannheim/MALDIcellassay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CeMOS-Mannheim/MALDIcellassay/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/github/CeMOS-Mannheim/MALDIcellassay/graph/badge.svg?token=URVX29WTDX)](https://codecov.io/github/CeMOS-Mannheim/MALDIcellassay) [![CRAN status](https://www.r-pkg.org/badges/version/MALDIcellassay)](https://CRAN.R-project.org/package=MALDIcellassay)

<!-- badges: end -->

<img src="inst/figures/MALDIcellassay_sticker.png" height="150"/>

# MALDIcellassay

Detects high variance signals and generates dose-response curves to further investigate candidate signals from MALDI cell based assays. More information in Unger et. Al. 2021

## Main functionality

The main function in this package is `fitCurve()` which will not only do logistic regression and variance filtering but also handle all preprocessing necessary: - Spectral alignment ("single-point re-calibration"). - Normalization by an internal standard or known endogenous signal

All these function can also be called on there own.

This package makes heavy use of the `MALDIquant` framework.
