# Apply normalization factors to spectra

Apply normalization factors to spectra

## Usage

``` r
normalizeByFactor(spec, factors)
```

## Arguments

- spec:

  List of MALDIquant::MassSpectrum or MALDIquant::MassPeaks

- factors:

  Numeric vector of normalization factors. See getNormFactors().

## Value

            List of normalized Spectra or Peaks

## Examples

``` r
#' data(Blank2022peaks)
normFactors <-  getNormFactors(Blank2022peaks, 
                               targetMz = 760.585, 
                               tol = 0.1, 
                               tolppm = FALSE)
normPeaks <- normalizeByFactor(Blank2022peaks, 
                               normFactors$norm_factor)
```
