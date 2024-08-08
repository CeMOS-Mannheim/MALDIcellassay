# MALDIcellassay
Detects high variance signals and generates dose-response curves to further investigate candidate signals from MALDI cell based assays.
More information in Unger et. Al. 2021

## Main functionality
The main function in this package is ```fitCurve()``` which will not only do logistic regression and variance filtering but also handle all preprocessing nessesary:
- Spectral alignment ("single-point recalibration").
- Normalization by an internal standard or known endogenous signal

All these function can also be called on there own.

This package makes heavy use of the ```MALDIquant``` framework.
